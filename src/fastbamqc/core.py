"""
fastbamqc.core
==============
Mapping and depth metrics for ONE BAM file with the two available analysis (whole-genome or amplicon).
"""
from __future__ import annotations

import math
import multiprocessing
from array import array
from collections import defaultdict

import numpy as np
import pysam

CHUNK_SIZE = 5_000_000

FLAG_UNMAPPED = 0x4

FLAG_PAIRED = 0x1
FLAG_PROPER_PAIR = 0x2
FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100
FLAG_QCFAIL = 0x200
FLAG_DUPLICATE = 0x400
FLAG_SUPPLEMENTARY = 0x800
FLAG_NOT_PRIMARY = FLAG_SECONDARY | FLAG_SUPPLEMENTARY


SHAPE_COLUMNS = ["median_coverage", "cv_coverage", "pct_above_0.2_mean", "fold_80"]
READ_STATS_COLUMNS = ["primary_mapped_reads", "pct_properly_paired", "pct_mapq0", "secondary_reads",
                    "supplementary_reads", "qcfail_reads", "marked_dup_reads", "pct_marked_dup",
                    "est_dup_reads", "pct_est_dup", "mean_coverage_nodup"]

GRID_BASE = (list(range(1, 21))
            + list(range(25, 101, 5))
            + list(range(120, 301, 20))
            + [350, 400, 500, 600, 800, 1000])


def make_grid(depths):
    """Sorted depths for the report sliders: base and user thresholds"""
    return sorted(set(GRID_BASE) | {int(d) for d in depths if 1 <= d <= HIST_CAP})


(S_SCANNED, S_SECONDARY, S_SUPPL, S_QCFAIL, S_PRIMARY, S_MAPQ0, S_PROPER,
    S_MARKED, S_EST, S_DUPLEN) = range(10)
N_STATS = 10

HIST_CAP = 1000
PROFILE_POINTS = 200
MAX_CONTIGS_IN_REPORT = 100

INTERPRETABLE = "INTERPRETABLE"
NOT_INTERPRETABLE = "NON_INTERPRETABLE"


def output_columns(depths, targeted=False, interpretability=False,
                uniformity=False, read_stats=False, nodup=False):
    """Column names of the final table"""
    if targeted:
        base = ["sample", "gene", "chrom", "start", "end", "length", "n_regions",
                "reads", "percentage_of_mapped", "mean_coverage"]
    else:
        base = ["sample", "total_read", "mapped_read", "percentage_mapped",
                "mean_coverage"]
    cols = base + [f"{d}X" for d in depths]
    if nodup and not targeted:
        cols += [f"{d}X_nodup" for d in depths]
    if uniformity:
        cols += SHAPE_COLUMNS
    if read_stats and not targeted:
        cols += READ_STATS_COLUMNS
    if interpretability:
        cols.append("interpretability")
    return cols


def rule_column(depth, nodup=False):
    """Name of the column the interpretability rule reads"""
    return f"{depth}X_nodup" if nodup else f"{depth}X"


def interpretability_label(row, depth, min_pct, nodup=False):
    """
    INTERPRETABLE if at least `min_pct` percent of the region (the whole genome,
    or one target) is covered at >= `depth` X, NON_INTERPRETABLE otherwise
    """
    return INTERPRETABLE if row[rule_column(depth, nodup)] >= min_pct else NOT_INTERPRETABLE


# ----------------------------------------------------------------------------
# BED reading
# ----------------------------------------------------------------------------
def _merge_intervals(intervals):
    """
    Merges intervals that overlap or touch each other
    """
    intervals.sort()
    merged = [list(intervals[0])]
    for start, end in intervals[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [tuple(m) for m in merged]


def _make_target(name, intervals):
    """
    Builds a target: {name, chrom, start, end, n_regions, intervals}.
    `intervals` = {contig: [(start, end), ...]}
    """
    chroms = list(intervals)
    if len(chroms) == 1:
        ivs = intervals[chroms[0]]
        chrom, start, end = chroms[0], ivs[0][0], ivs[-1][1]
    else:
        chrom, start, end = ",".join(chroms), "NA", "NA"
    return {"name": name, "chrom": chrom, "start": start, "end": end,
            "n_regions": sum(len(v) for v in intervals.values()),
            "intervals": intervals}


def load_bed(path, merge_same_name=True):
    """
    Reads a BED file
    Columns used: 1 chrom, 2 start (0-based), 3 end (exclusive), 4 gene name; without a 4th column the name is "chrom:start-end"

    merge_same_name=True: lines with the SAME NAME form ONE target (their
                        intervals are merged, so no double counting);
    merge_same_name=False: every BED line is an independent target, even if
                            the name is repeated (overlaps are then analysed
                            separately)
    """
    entries = []
    with open(path) as fh:
        for n, line in enumerate(fh, start=1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t") if "\t" in line else line.split()
            if len(fields) < 3:
                raise ValueError(f"{path} line {n}: at least 3 columns required (chrom, start, end)")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                raise ValueError(f"{path} line {n}: start/end are not integers")
            if start < 0 or end <= start:
                raise ValueError(f"{path} line {n}: invalid interval ({start}-{end})")
            name = fields[3].strip() if len(fields) > 3 and fields[3].strip() \
                else f"{fields[0]}:{start}-{end}"
            entries.append((fields[0], start, end, name))
    if not entries:
        raise ValueError(f"{path}: no region found")

    if not merge_same_name:
        return [_make_target(name, {chrom: [(start, end)]})
                for chrom, start, end, name in entries]

    groups = {}
    for chrom, start, end, name in entries:
        groups.setdefault(name, defaultdict(list))[chrom].append((start, end))
    return [_make_target(name, {c: _merge_intervals(iv) for c, iv in per_chrom.items()})
            for name, per_chrom in groups.items()]


def missing_contigs(targets, references):
    """BED contigs absent from the BAM header"""
    wanted = {chrom for t in targets for chrom in t["intervals"]}
    return sorted(wanted - set(references))


# ----------------------------------------------------------------------------
# Depth
# ----------------------------------------------------------------------------
def _chunk_coverage(bam, chrom, start, end, min_mapq, exclude_flags, count_floor,
                    read_stats=False, nodup=False):
    """
    Depth at every position
    `count_floor` prevents counting the same read in several windows

    Estimated duplicate = a primary read whose (start, end, strand, mate start,
    template length) is identical to that of an earlier read
    """

    starts, ends = array("q"), array("q")
    starts_nd, ends_nd = array("q"), array("q")
    n_reads = 0
    st = [0] * N_STATS if read_stats else None
    cur_start, seen = -1, set()

    for read in bam.fetch(chrom, start, end):
        flag = read.flag
        first = read.reference_start
        dup = False

        if st is not None and not flag & FLAG_UNMAPPED:
            counted = first >= count_floor
            primary = not flag & FLAG_NOT_PRIMARY
            marked = bool(flag & FLAG_DUPLICATE)
            est = False
            if primary:
                if first != cur_start:
                    cur_start = first
                    seen.clear()
                key = (read.reference_end, flag & FLAG_REVERSE,
                    read.next_reference_start, read.template_length)
                if key in seen:
                    est = True
                else:
                    seen.add(key)
            dup = marked or est
            if counted:
                st[S_SCANNED] += 1
                if flag & FLAG_QCFAIL:
                    st[S_QCFAIL] += 1
                if primary:
                    st[S_PRIMARY] += 1
                    if read.mapping_quality == 0:
                        st[S_MAPQ0] += 1
                    if (flag & 3) == 3:
                        st[S_PROPER] += 1
                    if marked:
                        st[S_MARKED] += 1
                    if est:
                        st[S_EST] += 1
                elif flag & FLAG_SECONDARY:
                    st[S_SECONDARY] += 1
                else:
                    st[S_SUPPL] += 1

        if flag & exclude_flags or read.mapping_quality < min_mapq:
            continue
        r_end = read.reference_end
        if r_end is None:
            continue
        r_start = first

        if r_start < start:
            r_start = start
        if r_end > end:
            r_end = end
        if r_end > r_start:
            starts.append(r_start - start)
            ends.append(r_end - start)
            if first >= count_floor:
                n_reads += 1
            if dup:
                st[S_DUPLEN] += r_end - r_start
            elif nodup:
                starts_nd.append(r_start - start)
                ends_nd.append(r_end - start)

    if not starts:
        return None, 0, st, None

    size = end - start
    return _depth_from_intervals(starts, ends, size), n_reads, st, \
        (_depth_from_intervals(starts_nd, ends_nd, size) if starts_nd else None)


def _depth_from_intervals(starts, ends, size):
    """Depth per position from read starts/ends"""
    diff = (np.bincount(np.frombuffer(starts, dtype=np.int64), minlength=size + 1)
            - np.bincount(np.frombuffer(ends, dtype=np.int64), minlength=size + 1))
    return np.cumsum(diff[:size])


class _Acc:
    """Running sums of the results of windows"""

    def __init__(self, depths, hist_cap=None, read_stats=False, want_sq=False, nodup=False, profile=False):
        self.depths = list(depths)
        self.length = self.sum_depth = self.reads = self.sum_sq = 0
        self.n_ge = {d: 0 for d in self.depths}
        self.stats = [0] * N_STATS if read_stats else None
        self.hist = np.zeros(hist_cap + 1, dtype=np.int64) if hist_cap else None
        self.n_ge_nd = {d: 0 for d in self.depths} if nodup else None
        self.hist_nd = np.zeros(hist_cap + 1, dtype=np.int64) if (nodup and hist_cap) else None
        self.pieces = [] if profile else None

    def add_empty(self, size):
        """`size` positions without any read (0X)."""
        self.length += size
        if self.hist is not None:
            self.hist[0] += size
        if self.hist_nd is not None:
            self.hist_nd[0] += size
        if self.pieces is not None:
            self.pieces.append(np.zeros(size, dtype=np.int32))

    def add(self, p):
        """Adds the partial result `p` of one window (see _scan_window)."""
        self.length += p["length"]
        self.sum_depth += p["sum_depth"]
        self.sum_sq += p["sum_sq"]
        self.reads += p["reads"]
        for d in self.depths:
            self.n_ge[d] += p["n_ge"][d]
        if self.stats is not None:
            for k, v in enumerate(p["stats"]):
                self.stats[k] += v
        if self.hist is not None:
            self.hist += p["hist"]
        if self.n_ge_nd is not None:
            for d in self.depths:
                self.n_ge_nd[d] += p["n_ge_nd"][d]
            if self.hist_nd is not None:
                self.hist_nd += p["hist_nd"]
        if self.pieces is not None:
            self.pieces.append(p["piece"])


def _scan_window(bam, chrom, w_start, w_end, floor, depths, min_mapq, exclude_flags,
                hist_cap=None, profile=False, read_stats=False, want_sq=False, nodup=False):
    """
    Everything computed for one window
    """
    size = w_end - w_start
    cov, n, st, cov_nd = _chunk_coverage(bam, chrom, w_start, w_end, min_mapq, exclude_flags,
                                        floor, read_stats, nodup)
    acc = _Acc(depths, hist_cap, read_stats, want_sq, nodup, profile)
    acc.reads, acc.length = n, size
    if st is not None:
        acc.stats = st
    if cov is None:
        if acc.hist is not None:
            acc.hist[0] += size
        if acc.hist_nd is not None:
            acc.hist_nd[0] += size
        if acc.pieces is not None:
            acc.pieces.append(np.zeros(size, dtype=np.int32))
    else:
        acc.sum_depth = int(cov.sum())
        if want_sq:
            acc.sum_sq = int(np.dot(cov, cov))
        for d in depths:
            acc.n_ge[d] = int(np.count_nonzero(cov >= d))
        if acc.hist is not None:
            acc.hist += np.bincount(np.minimum(cov, hist_cap), minlength=hist_cap + 1)
        if nodup:
            if cov_nd is None:
                if acc.hist_nd is not None:
                    acc.hist_nd[0] += size
            else:
                for d in depths:
                    acc.n_ge_nd[d] = int(np.count_nonzero(cov_nd >= d))
                if acc.hist_nd is not None:
                    acc.hist_nd += np.bincount(np.minimum(cov_nd, hist_cap), minlength=hist_cap + 1)
        if acc.pieces is not None:
            acc.pieces.append(cov.astype(np.int32))
    return {"length": size, "sum_depth": acc.sum_depth, "sum_sq": acc.sum_sq, "reads": n,
            "n_ge": acc.n_ge, "stats": acc.stats, "hist": acc.hist, "n_ge_nd": acc.n_ge_nd,
            "hist_nd": acc.hist_nd, "piece": acc.pieces[0] if acc.pieces else None}


def _windows(intervals, mapped_per_contig, size=CHUNK_SIZE):
    """
    Plan of the windows over {contig: [(start, end), ...]}
    """
    wins, empty = [], []
    for chrom, ivs in intervals.items():
        if mapped_per_contig.get(chrom, 0) == 0:
            empty.append((chrom, sum(e - s for s, e in ivs)))
            continue
        prev_end = 0
        for r_start, r_end in ivs:
            floor = prev_end
            for w_start in range(r_start, r_end, size):
                w_end = min(w_start + size, r_end)
                wins.append((chrom, w_start, w_end, floor))
                floor = w_end
            prev_end = r_end
    return wins, empty


def _scan_intervals(bam, intervals, mapped_per_contig, depths, min_mapq, exclude_flags,
                    hist_cap=None, profile=False, read_stats=False, want_sq=False,
                    nodup=False):
    """
    Accumulates the depth over a set of intervals {contig: [(start, end), ...]}
    (sorted and disjoint per contig), window after window
    """
    acc = _Acc(depths, hist_cap, read_stats, want_sq, nodup, profile)
    for chrom, ivs in intervals.items():
        wins, empty = _windows({chrom: ivs}, mapped_per_contig)
        for _, span in empty:
            acc.add_empty(span)
        for c, w_start, w_end, floor in wins:
            acc.add(_scan_window(bam, c, w_start, w_end, floor, depths, min_mapq, exclude_flags,
                                hist_cap, profile, read_stats, want_sq, nodup))
    return {"length": acc.length, "sum_depth": acc.sum_depth, "reads": acc.reads, "n_ge": acc.n_ge,
            "hist": acc.hist, "sum_sq": acc.sum_sq, "stats": acc.stats,
            "n_ge_nd": acc.n_ge_nd, "hist_nd": acc.hist_nd,
            "profile": np.concatenate(acc.pieces) if acc.pieces else None}


def depth_shape(hist, total, sum_depth, sum_sq, cap=None):
    """
    Shape of the depth distribution
    """
    na = {c: None for c in SHAPE_COLUMNS}
    if total <= 0 or sum_depth <= 0:
        return na
    mean = sum_depth / total
    cv = math.sqrt(max(sum_sq / total - mean * mean, 0.0)) / mean
    cum = np.cumsum(hist)
    median = int(np.searchsorted(cum, 0.5 * total, side="left"))
    p20 = int(np.searchsorted(cum, 0.2 * total, side="right"))
    thr = math.ceil(0.2 * mean)
    if thr == 0:
        above = 100.0
    elif cap is not None and thr > cap:
        above = None
    elif thr >= len(cum):
        above = 0.0
    else:
        above = (total - int(cum[thr - 1])) / total * 100
    return {"median_coverage": float(median),
            "cv_coverage": round(cv, 2),
            "pct_above_0.2_mean": None if above is None else round(above, 2),
            "fold_80": round(mean / p20, 2) if p20 > 0 else None}


def _read_stats_result(st, sum_depth, genome_len):
    """Read-level columns (READ_STATS_COLUMNS) from _chunk_coverage"""
    prim = st[S_PRIMARY]
    pct = lambda n: round(n / prim * 100, 2) if prim else 0.0
    return {"primary_mapped_reads": prim,
            "pct_properly_paired": pct(st[S_PROPER]), "pct_mapq0": pct(st[S_MAPQ0]),
            "secondary_reads": st[S_SECONDARY], "supplementary_reads": st[S_SUPPL],
            "qcfail_reads": st[S_QCFAIL],
            "marked_dup_reads": st[S_MARKED], "pct_marked_dup": pct(st[S_MARKED]),
            "est_dup_reads": st[S_EST], "pct_est_dup": pct(st[S_EST]),
            "mean_coverage_nodup": round((sum_depth - st[S_DUPLEN]) / genome_len, 2)}


def _grid_hundredths(n_ge, total):
    """Percentage of positions covered, rounded to 2 decimals as in the table"""
    return int(round(round(n_ge / total * 100, 2) * 100))


def _grid_from_hist(hist, total, xs):
    """Grid values from a depth histogram"""
    cum = np.cumsum(hist)
    return [_grid_hundredths(total - int(cum[x - 1]) if x - 1 < len(cum) else 0, total) for x in xs]


def _bin_profile(depth, points):
    """
    Depth profile of a target, reduced to at most `points` values
    """
    n = len(depth)
    if n <= points:
        values = depth.astype(np.float64)
    else:
        edges = np.linspace(0, n, points + 1).astype(np.int64)
        values = np.add.reduceat(depth.astype(np.float64), edges[:-1]) / np.diff(edges)
    return [int(v) if float(v).is_integer() else round(float(v), 1) for v in values]


def _open_indexed(bam_path, threads):
    """Opens a BAM file and checks that it is indexed"""
    bam = pysam.AlignmentFile(bam_path, "rb", threads=threads)
    if not bam.has_index():
        bam.close()
        raise RuntimeError("index missing (.bai/.csi): run `samtools index`")
    return bam


_HANDLES = {}
MIN_WINDOW = 250_000


def _window_size(genome_len, workers):
    """Window size: CHUNK_SIZE, or smaller when several workers must share a small genome"""
    if workers <= 1:
        return CHUNK_SIZE
    return max(MIN_WINDOW, min(CHUNK_SIZE, genome_len // (workers * 8)))


def _window_job(task):
    """Pool worker: one window of one BAM"""
    path, chrom, w_start, w_end, floor, params = task
    bam = _HANDLES.get(path)
    if bam is None:
        bam = _HANDLES[path] = pysam.AlignmentFile(path, "rb")
    return chrom, _scan_window(bam, chrom, w_start, w_end, floor, *params)


# ----------------------------------------------------------------------------
# Whole-genome analysis: 1 result per BAM
# ----------------------------------------------------------------------------
def analyze_bam(bam_path, depths=(1, 5, 10), threads=1, min_mapq=0,
                exclude_flags=0, contigs=None, detail=False,
                uniformity=False, read_stats=False, grid_xs=None, nodup=False, workers=1):
    """
    Analyses a sorted and indexed BAM over the whole genome
    """
    exclude_flags |= FLAG_UNMAPPED

    with _open_indexed(bam_path, threads) as bam:
        # 1. Mapping statistics (from the index file)
        mapped = bam.mapped
        total = mapped + bam.unmapped
        pct_mapped = mapped / total * 100 if total else 0.0
        mapped_per_contig = {s.contig: s.mapped for s in bam.get_index_statistics()}

        # 2. Selected contigs
        wanted = set(contigs) if contigs else None
        if wanted:
            unknown = wanted - set(bam.references)
            if unknown:
                raise ValueError(f"contig(s) not found in the header: {sorted(unknown)}")
        intervals = {c: [(0, n)] for c, n in zip(bam.references, bam.lengths)
                    if n > 0 and (wanted is None or c in wanted)}

    # 3. Depth, window by window
    need_hist = detail or uniformity or bool(grid_xs)
    hist_cap = HIST_CAP if need_hist else None
    with_stats = read_stats or nodup
    make = lambda: _Acc(depths, hist_cap, with_stats, uniformity, nodup)
    total_acc = make()
    per = {c: make() for c in intervals} if detail else None
    wins, empty = _windows(intervals, mapped_per_contig, _window_size(sum(e - s for ivs in intervals.values() for s, e in ivs), workers))
    for chrom, span in empty:
        total_acc.add_empty(span)
        if detail:
            per[chrom].add_empty(span)

    def take(chrom, p):
        total_acc.add(p)
        if detail:
            per[chrom].add(p)

    params = (list(depths), min_mapq, exclude_flags, hist_cap, False, with_stats, uniformity, nodup)
    if workers > 1 and len(wins) > 1:
        tasks = [(bam_path, c, ws, we, fl, params) for c, ws, we, fl in wins]
        with multiprocessing.Pool(min(workers, len(wins))) as pool:
            for chrom, p in pool.imap_unordered(_window_job, tasks, chunksize=1):
                take(chrom, p)
    else:
        with _open_indexed(bam_path, threads) as bam:
            for c, ws, we, fl in wins:
                take(c, _scan_window(bam, c, ws, we, fl, *params))

    target_len, sum_depth, sum_sq = total_acc.length, total_acc.sum_depth, total_acc.sum_sq
    n_ge, hist, hist_nd, n_ge_nd = total_acc.n_ge, total_acc.hist, total_acc.hist_nd, total_acc.n_ge_nd
    stats = total_acc.stats if with_stats else [0] * N_STATS
    per_contig = []
    if detail:
        for chrom, a in per.items():
            per_contig.append({
                "name": chrom, "length": a.length, "mapped": mapped_per_contig.get(chrom, 0),
                "mean": round(a.sum_depth / a.length, 2),
                "pct": [round(a.n_ge[d] / a.length * 100, 2) for d in depths]})

    if target_len == 0:
        raise ValueError("genome length is zero (contigs missing or filtered out)")

    result = {
        "total_read": total,
        "mapped_read": mapped,
        "percentage_mapped": round(pct_mapped, 2),
        "mean_coverage": round(sum_depth / target_len, 2),
    }
    for d in depths:
        result[f"{d}X"] = round(n_ge[d] / target_len * 100, 2)
    if nodup:
        for d in depths:
            result[f"{d}X_nodup"] = round(n_ge_nd[d] / target_len * 100, 2)
    if uniformity:
        result.update(depth_shape(hist, target_len, sum_depth, sum_sq, cap=HIST_CAP))
    if read_stats:
        result.update(_read_stats_result(stats, sum_depth, target_len))

    detail_out = None
    if detail:
        largest = sorted(range(len(per_contig)), key=lambda i: -per_contig[i]["length"])
        keep = sorted(largest[:MAX_CONTIGS_IN_REPORT])
        nonzero = np.flatnonzero(hist)
        detail_out = {
            "hist_cap": HIST_CAP,
            "hist": hist[: (int(nonzero[-1]) + 1 if len(nonzero) else 1)].tolist(),
            "n_contigs": len(per_contig),
            "contigs": [per_contig[i] for i in keep],
        }
        if nodup:
            nz = np.flatnonzero(hist_nd)
            detail_out["hist_nd"] = hist_nd[: (int(nz[-1]) + 1 if len(nz) else 1)].tolist()
        if grid_xs:
            detail_out["grid"] = _grid_from_hist(hist_nd if nodup else hist, target_len, grid_xs)
    return result, detail_out


# ----------------------------------------------------------------------------
# Amplicon / panel analysis: 1 result per target (gene) and per BAM
# ----------------------------------------------------------------------------
def analyze_targets(bam_path, targets, depths=(1, 5, 10), threads=1,
                    min_mapq=0, exclude_flags=0, detail=False, uniformity=False,
                    grid_xs=None):
    """
    Analyses a sorted and indexed BAM over each target (amplicon)
    """
    exclude_flags |= FLAG_UNMAPPED
    rows, profiles, grids = [], [], []

    with _open_indexed(bam_path, threads) as bam:
        mapped = bam.mapped
        mapped_per_contig = {s.contig: s.mapped for s in bam.get_index_statistics()}

        unknown = missing_contigs(targets, bam.references)
        if unknown:
            raise ValueError(
                f"BED contig(s) not found in the BAM: {unknown[:5]}"
                f"{'...' if len(unknown) > 5 else ''} "
                "(different 'chr' prefix between BED and BAM?)")

        for t in targets:
            intervals = {}
            for chrom, ivs in t["intervals"].items():
                length = bam.get_reference_length(chrom)
                clipped = [(s, min(e, length)) for s, e in ivs if s < length]
                if clipped:
                    intervals[chrom] = clipped

            r = _scan_intervals(bam, intervals, mapped_per_contig, depths, min_mapq,
                                exclude_flags, profile=detail or uniformity or bool(grid_xs),
                                want_sq=uniformity)
            tlen = r["length"]
            if tlen == 0:
                raise ValueError(f"target '{t['name']}' lies entirely outside the BAM contigs")

            row = {
                "gene": t["name"], "chrom": t["chrom"], "start": t["start"],
                "end": t["end"], "length": tlen, "n_regions": t["n_regions"],
                "reads": r["reads"],
                "percentage_of_mapped": round(r["reads"] / mapped * 100, 2) if mapped else 0.0,
                "mean_coverage": round(r["sum_depth"] / tlen, 2),
            }
            for d in depths:
                row[f"{d}X"] = round(r["n_ge"][d] / tlen * 100, 2)
            if uniformity:
                prof = r["profile"].astype(np.int64)
                row.update(depth_shape(np.bincount(prof), tlen, r["sum_depth"], r["sum_sq"]))
            rows.append(row)
            if grid_xs:
                srt = np.sort(r["profile"])
                grids.append([_grid_hundredths(tlen - int(np.searchsorted(srt, x, side="left")), tlen)
                            for x in grid_xs])

            if detail:
                region_sizes = [e - s for ivs in intervals.values() for s, e in ivs]
                breaks = [round(float(x) / tlen, 4) for x in np.cumsum(region_sizes)[:-1]]
                profiles.append({"gene": t["name"], "length": tlen, "breaks": breaks,
                                "depth": _bin_profile(r["profile"], PROFILE_POINTS)})
    extra = None
    if detail:
        extra = {"profiles": profiles}
        if grid_xs:
            extra["grid"] = grids
    return rows, extra
