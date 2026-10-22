"""
fastbamqc.cli
=============
Command-line interface: finds the BAM files, runs the analysis in parallel
and writes the result tables.

Examples:
    fastbamqc -i ./bams -f "*_sorted.bam" -o ./results -F tsv -c 16 -t 2
    fastbamqc -i ./bams -f "*_sorted.bam" -o ./results --bed amplicons.bed   # amplicon-seq
"""
from __future__ import annotations

import argparse
import csv
import glob
import logging
import os
import shlex
import sys
from multiprocessing import Pool

import pysam

from . import __version__
from .core import (INTERPRETABLE, analyze_bam, analyze_targets, interpretability_label,
                   load_bed, make_grid, missing_contigs, output_columns)
from .report import write_report

log = logging.getLogger("fastbamqc")


# ----------------------------------------------------------------------------
# Utilities
# ----------------------------------------------------------------------------
def find_bams(input_dir, pattern):
    """Sorted list of the files of the directory matching the pattern (e.g. *_sorted.bam)."""
    return sorted(glob.glob(os.path.join(input_dir, pattern)))


def sample_name(path, pattern):
    """
    Sample name = file name without the fixed part of the pattern.
    "S1_sorted.bam" with "*_sorted.bam" -> "S1".
    (The suffix is removed only at the END of the name, not at every occurrence.)
    """
    name = os.path.basename(path)
    # The suffix only makes sense if the pattern contains a star; without a star
    # (e.g. "S1.bam") it would equal the whole file name and the sample would be empty.
    suffix = pattern.split("*")[-1] if "*" in pattern else ""
    if suffix and name.endswith(suffix) and len(name) > len(suffix):
        return name[: -len(suffix)]
    return os.path.splitext(name)[0]         # fallback: file name without extension


def _depth_list(text):
    """Parse '1,5,10,30' -> [1, 5, 10, 30] (argparse type)."""
    try:
        depths = sorted({int(x) for x in text.split(",") if x.strip()})
    except ValueError:
        raise argparse.ArgumentTypeError("expected integers separated by commas")
    if not depths or depths[0] < 1:
        raise argparse.ArgumentTypeError("thresholds must be integers >= 1")
    return depths


def _worker(task):
    """
    Runs in a child process: analyses ONE BAM and returns (sample, rows, extra, error):
    rows = LIST of rows (whole genome: 1 row; amplicon/panel: 1 row per target),
    extra = detailed data for the report (None unless --report).
    Exceptions are never raised here: the error is sent back to the parent
    process, which decides what to do (one failing BAM must not block the others).
    """
    bam_path, sample, kwargs, targets, interp = task
    try:
        if targets is None:
            result, extra = analyze_bam(bam_path, **kwargs)
            rows = [result]
        else:
            rows, extra = analyze_targets(bam_path, targets, **kwargs)
        if interp is not None:               # (depth, percentage): add the verdict column
            for r in rows:
                r["interpretability"] = interpretability_label(r, *interp)
        return sample, [{"sample": sample, **r} for r in rows], extra, None
    except Exception as exc:                 # noqa: BLE001 (deliberately broad)
        return sample, None, None, f"{type(exc).__name__}: {exc}"


def _run(tasks, jobs):
    """Iterates over the results as they come in."""
    if jobs == 1:                            # no Pool: simpler to debug
        yield from map(_worker, tasks)
    else:
        with Pool(processes=jobs) as pool:
            yield from pool.imap_unordered(_worker, tasks)


def _fmt(value):
    """Decimal numbers -> exactly 2 decimals (98.2 -> '98.20'); None -> 'NA'; else unchanged."""
    if value is None:
        return "NA"
    return f"{value:.2f}" if isinstance(value, float) else value


def write_table(rows, columns, path, fmt):
    """Writes the table as tsv / csv / xlsx (decimals with 2 digits)."""
    if fmt in ("tsv", "csv"):
        with open(path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=columns,
                                    delimiter="\t" if fmt == "tsv" else ",")
            writer.writeheader()
            writer.writerows({k: _fmt(v) for k, v in row.items()} for row in rows)
    elif fmt == "xlsx":
        from openpyxl import Workbook        # late import: optional dependency
        wb = Workbook()
        ws = wb.active
        ws.append(columns)
        for row in rows:
            ws.append([("NA" if row[c] is None else row[c]) for c in columns])
            for cell in ws[ws.max_row]:      # numeric cells, displayed with 2 decimals
                if isinstance(cell.value, float):
                    cell.number_format = "0.00"
        wb.save(path)


# ----------------------------------------------------------------------------
# Help display
# ----------------------------------------------------------------------------
DEFAULT_INTERP_DEPTH = 10      # interpretability rule used unless --interp-depth / --interp-pct / --no-interp
DEFAULT_INTERP_PCT = 75.0

BANNER = r"""
 _____ _    ____ _____ ____    _    __  __  ___   ____
|  ___/ \  / ___|_   _| __ )  / \  |  \/  |/ _ \ / ___|
| |_ / _ \ \___ \ | | |  _ \ / _ \ | |\/| | | | | |
|  _/ ___ \ ___) || | | |_) / ___ \| |  | | |_| | |___
|_|/_/   \_\____/ |_| |____/_/   \_\_|  |_|\__\_\\____|

"""
TAGLINE = "Fast mapping & depth QC for sorted, indexed BAM files"


def _use_color():
    """Colours only on a real terminal, and never if NO_COLOR is set."""
    return sys.stdout.isatty() and "NO_COLOR" not in os.environ


def _banner():
    """Big software name + version + tagline, shown at the very top of the help."""
    art = "\n".join("  " + line for line in BANNER.strip("\n").splitlines())
    sub = f"  v{__version__}  |  {TAGLINE}"
    if _use_color():
        art = f"\033[1;36m{art}\033[0m"      # bold cyan
        sub = f"\033[2m{sub}\033[0m"         # dim
    return f"\n{art}\n\n{sub}\n\n"


class _HelpFormatter(argparse.RawDescriptionHelpFormatter):
    """
    Airy help: upper-case section titles, a blank line after every option,
    line breaks of descriptions/epilog preserved.
    """
    def __init__(self, prog):
        super().__init__(prog, max_help_position=30)

    def start_section(self, heading):
        super().start_section(heading.upper())

    def add_usage(self, usage, actions, groups, prefix=None):
        super().add_usage(usage, actions, groups, prefix or "USAGE:\n  ")

    def _format_action(self, action):
        return super()._format_action(action) + "\n"

    def _format_action_invocation(self, action):
        # "-i, --input-dir DIRECTORY" instead of "-i DIRECTORY, --input-dir DIRECTORY"
        if not action.option_strings or action.nargs == 0:
            return super()._format_action_invocation(action)
        metavar = self._format_args(action, self._get_default_metavar_for_optional(action))
        return ", ".join(action.option_strings) + " " + metavar


class _Parser(argparse.ArgumentParser):
    """ArgumentParser that prints the banner above the help."""
    def format_help(self):
        return _banner() + super().format_help()


EPILOG = """\
EXAMPLES:

  Whole genome, 16 cores in total, 2 threads per BAM (so 8 BAM files at a time):
    fastbamqc -i ./bams -f "*_sorted.bam" -o ./results -c 16 -t 2

  Whole genome, selected chromosomes only, CSV output:
    fastbamqc -i ./bams -o ./results -F csv --contigs chr1,chr2,chr3,chrX

  Amplicon-seq / panel, depth computed on the BED targets only:
    fastbamqc -i ./bams -f "*_sorted.bam" -o ./results --bed amplicons.bed -d 1,10,100,500

  Amplicon-seq, BED lines sharing the same gene name analysed separately:
    fastbamqc -i ./bams -o ./results --bed amplicons.bed --same-gene separate

  Interpretability is on by default: a row is INTERPRETABLE if at least 75% of it is
  covered at >= 10X. Change the rule, or switch it off:
    fastbamqc -i ./bams -o ./results --bed amplicons.bed --interp-depth 20 --interp-pct 90
    fastbamqc -i ./bams -o ./results --no-interp

  Same analysis, plus an interactive HTML report (<out-dir>/<out-name>.html):
    fastbamqc -i ./bams -o ./results --bed amplicons.bed --report


OUTPUT FILES:

  Whole genome: a single table, one row per sample
    sample, total_read, mapped_read, percentage_mapped, mean_coverage, <d>X

  With --bed: one row per gene/amplicon and per sample
    <out-dir>/<out-name>.<fmt>                        all rows of all samples
    <out-dir>/per_sample/<sample>.amplicons.<fmt>     one file per sample
    sample, gene, chrom, start, end, length, n_regions, reads,
    percentage_of_mapped, mean_coverage, <d>X

  Unless --no-interp is given, a last column "interpretability" is added to every
  row (INTERPRETABLE or NON_INTERPRETABLE).

  With --report, a single self-contained HTML file <out-dir>/<out-name>.html:
  overview of all samples + a page per sample (works offline, can be e-mailed),
  with sliders to change the interpretability rule live. In whole-genome mode the
  report also has a duplicates section, which turns on --read-stats.


REQUIREMENTS:

  BAM files must be coordinate-sorted and indexed (samtools index).
  Always quote the -f pattern (e.g. "*.bam") so the shell does not expand it.
"""


# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------
def _normalize_argv(argv):
    """
    Accepts --out_dir as --out-dir: '_' is replaced by '-' in the NAME of long
    options only (never in their values). This avoids declaring every option
    twice, which would clutter the help.
    """
    fixed = []
    for arg in argv:
        if arg.startswith("--"):
            name, sep, value = arg.partition("=")
            arg = name.replace("_", "-") + sep + value
        fixed.append(arg)
    return fixed


def parse_args(argv=None):
    p = _Parser(
        prog="fastbamqc",
        usage="fastbamqc -i DIRECTORY [-f PATTERN] [-b FILE.bed] [options]",
        description="Fast mapping and depth metrics for sorted and indexed BAM files.",
        epilog=EPILOG,
        formatter_class=_HelpFormatter,
        add_help=False,                      # -h is added below, in "other options"
    )

    # --- Required argument ---
    g = p.add_argument_group("required")
    g.add_argument("-i", "--input-dir", dest="input_dir", required=True, metavar="DIRECTORY",
                   help="Directory containing the BAM files.")

    # --- Input files ---
    g = p.add_argument_group("input files")
    g.add_argument("-f", "--pattern", default="*.bam", metavar="PATTERN",
                   help='Pattern of the BAM files to analyse, e.g. "*_sorted.bam" '
                        '(default: "*.bam"). Sample names are the file names without '
                        "the fixed part of the pattern.")

    # --- Analysis mode ---
    g = p.add_argument_group(
        "analysis mode",
        "Default   : WHOLE GENOME. Depth is computed over every contig of the BAM\n"
        "            header; one result row per sample.\n"
        "With --bed: AMPLICON / PANEL. Depth is computed over the BED targets only;\n"
        "            one result row per gene (BED column 4) and per sample.\n")
    g.add_argument("-b", "--bed", metavar="FILE.bed",
                   help="Enable amplicon/panel mode. BED columns: 1 chrom, 2 start (0-based), "
                        "3 end, 4 gene name. Without a 4th column, targets are named "
                        "chrom:start-end.")
    g.add_argument("--same-gene", choices=["merge", "separate"], metavar="{merge,separate}",
                   help="With --bed: how to handle several BED lines sharing the same gene "
                        "name. 'merge' (default) groups them into ONE target whose regions "
                        "are merged; 'separate' treats every line as an independent target.")
    g.add_argument("--contigs", metavar="CONTIGS",
                   type=lambda s: [c for c in s.split(",") if c],
                   help="Whole-genome mode only: restrict the analysis to these "
                        "comma-separated contigs, e.g. chr1,chr2,chrX (default: all contigs "
                        "of the BAM header). Cannot be combined with --bed.")

    # --- Output ---
    g = p.add_argument_group("output")
    g.add_argument("-o", "--out-dir", dest="out_dir", default=".", metavar="DIRECTORY",
                   help="Output directory, created if needed (default: current directory).")
    g.add_argument("-F", "--out-format", dest="out_format", default="tsv",
                   choices=["tsv", "csv", "xlsx"], metavar="{tsv,csv,xlsx}",
                   help="Format of the output tables (default: tsv; xlsx requires openpyxl).")
    g.add_argument("-n", "--out-name", dest="out_name", default="fastbamqc_summary",
                   metavar="NAME",
                   help="Name of the main output file, without extension "
                        "(default: fastbamqc_summary).")

    # --- Depth and filters ---
    g = p.add_argument_group("depth and filters")
    g.add_argument("-d", "--depths", type=_depth_list, default=[1, 5, 10], metavar="THRESHOLDS",
                   help="Comma-separated depth thresholds. One <d>X column is reported per "
                        "threshold: the percentage of positions covered by at least d reads "
                        "(default: 1,5,10).")
    g.add_argument("-q", "--min-mapq", type=int, default=0, metavar="MAPQ",
                   help="Minimum mapping quality (MAPQ) for a read to be counted in the "
                        "depth (default: 0).")
    g.add_argument("--exclude-flags", type=lambda s: int(s, 0), default=0, metavar="FLAGS",
                   help="SAM flags to exclude from the depth, decimal or 0x.. format "
                        "(default: 0; unmapped reads are always excluded). E.g. 1024 = "
                        "duplicates, 3584 = duplicates + secondary + supplementary. "
                        "Avoid excluding duplicates on amplicon data, where they are expected.")

    # --- Extra metrics ---
    g = p.add_argument_group(
        "extra metrics",
        "Optional columns, off by default (the default output and speed are unchanged).\n")
    g.add_argument("--read-stats", action="store_true",
                   help="Whole genome only. Adds read-level QC columns: %% properly paired, "
                        "%% MAPQ 0, secondary / supplementary / QC-fail reads, duplicates "
                        "flagged in the BAM (marked_dup_*) AND duplicates estimated from the "
                        "alignments (est_dup_*: same start, end, strand and mate position as "
                        "an earlier read), and mean_coverage_nodup (mean depth without "
                        "duplicates). Needs one extra check per read: slower (see the README). "
                        "Turned on automatically by --report in whole-genome mode.")
    g.add_argument("--uniformity", action="store_true",
                   help="Adds depth-evenness columns (per sample, or per target in amplicon "
                        "mode): median_coverage, cv_coverage (std/mean), pct_above_0.2_mean "
                        "(%% of positions with depth >= 0.2 x mean) and fold_80 (mean divided "
                        "by the depth reached by 80%% of the positions; 1 = perfectly even).")

    # --- Interpretability ---
    g = p.add_argument_group(
        "interpretability",
        "On by default. Flags every result row (one sample in whole-genome mode, one gene in\n"
        "amplicon mode) as INTERPRETABLE or NON_INTERPRETABLE in an extra column.\n"
        "A row is INTERPRETABLE when at least --interp-pct percent of its region is covered\n"
        "by at least --interp-depth reads (defaults: 75% at 10X).\n")
    g.add_argument("--interp-depth", type=int, metavar="DEPTH", default=None,
                   help="Depth (in X) a position must reach to count as covered "
                        "(default: 10). Added to the <d>X columns if missing from --depths.")
    g.add_argument("--interp-pct", type=float, metavar="PERCENT", default=None,
                   help="Minimum percentage of the region (or of the genome) that must reach "
                        "--interp-depth (default: 75). The decision uses the 2-decimal value "
                        "shown in the table.")
    g.add_argument("--no-interp", action="store_true",
                   help="Do not compute the interpretability column (and no rule sliders in "
                        "the report).")

    # --- Report ---
    g = p.add_argument_group("report")
    g.add_argument("--report", action="store_true",
                   help="Also write a self-contained HTML report with tables and charts: "
                        "an overview of all samples and a detailed section per sample "
                        "(depth distribution and per-contig depth in whole-genome mode, "
                        "depth profile of every target in amplicon mode), sliders to change the "
                        "interpretability rule live, and (whole genome) a duplicates section. "
                        "Written to <out-dir>/<out-name>.html.")
    g.add_argument("--report-title", metavar="TEXT", default="fastbamqc report",
                   help='Title shown at the top of the report (default: "fastbamqc report").')

    # --- Performance ---
    g = p.add_argument_group("performance")
    g.add_argument("-c", "--max-cores", type=int, default=4, metavar="N",
                   help="Total number of CPU cores to use (default: 4).")
    g.add_argument("-t", "--threads-per-bam", type=int, default=1, metavar="N",
                   help="Decompression threads per BAM file (default: 1). "
                        "BAM files processed in parallel = max-cores // threads-per-bam.")

    # --- Other ---
    g = p.add_argument_group("other options")
    g.add_argument("--log-file", metavar="FILE",
                   help="Also write the progress messages to this file.")
    g.add_argument("-V", "--version", action="version", version=f"%(prog)s {__version__}",
                   help="Show the version and exit.")
    g.add_argument("-h", "--help", action="help", help="Show this help message and exit.")

    if argv is None:
        argv = sys.argv[1:]
    return p.parse_args(_normalize_argv(argv))


# ----------------------------------------------------------------------------
# Main program
# ----------------------------------------------------------------------------
def main(argv=None):
    argv_list = sys.argv[1:] if argv is None else list(argv)
    args = parse_args(argv_list)

    # Logging: stderr (+ optional file)
    handlers = [logging.StreamHandler(sys.stderr)]
    if args.log_file:
        handlers.append(logging.FileHandler(args.log_file))
    logging.basicConfig(level=logging.INFO, handlers=handlers,
                        format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S")

    # --- Preliminary checks (fail early, before the long computations) ---
    if not os.path.isdir(args.input_dir):
        log.error("Directory not found: %s", args.input_dir)
        return 2
    if args.max_cores < 1 or args.threads_per_bam < 1:
        log.error("--max-cores and --threads-per-bam must be >= 1")
        return 2
    if args.bed and args.contigs:
        log.error("--bed and --contigs are mutually exclusive (the BED already defines the targets)")
        return 2
    if args.read_stats and args.bed:
        log.error("--read-stats is only available in whole-genome mode (without --bed): "
                  "duplicates are expected in amplicon data")
        return 2
    if args.same_gene and not args.bed:
        log.error("--same-gene only makes sense together with --bed")
        return 2
    if args.no_interp and (args.interp_depth is not None or args.interp_pct is not None):
        log.error("--no-interp cannot be combined with --interp-depth / --interp-pct")
        return 2
    interp, depths = None, args.depths
    if not args.no_interp:
        interp_depth = DEFAULT_INTERP_DEPTH if args.interp_depth is None else args.interp_depth
        interp_pct = DEFAULT_INTERP_PCT if args.interp_pct is None else args.interp_pct
        if interp_depth < 1:
            log.error("--interp-depth must be an integer >= 1")
            return 2
        if not 0 < interp_pct <= 100:
            log.error("--interp-pct must be a percentage in ]0, 100]")
            return 2
        interp = (interp_depth, interp_pct)
        if interp_depth not in args.depths:        # the verdict needs this <d>X column
            depths = sorted(set(args.depths) | {interp_depth})
            log.info("Depth threshold %dX added to the depth columns "
                     "(needed for the interpretability rule)", interp_depth)
    if args.out_format == "xlsx":
        try:
            import openpyxl  # noqa: F401
        except ImportError:
            log.error("xlsx output requires openpyxl (pip install 'fastbamqc[xlsx]')")
            return 2

    bams = find_bams(args.input_dir, args.pattern)
    if not bams:
        log.error("No file matching '%s' in %s", args.pattern, args.input_dir)
        return 2

    # --- Amplicon/panel mode: read and validate the BED (once) ---
    targets = None
    if args.bed:
        merge = (args.same_gene or "merge") == "merge"
        try:
            targets = load_bed(args.bed, merge_same_name=merge)
            # Check the contig names against the header of the first BAM: fail right
            # away instead of getting N identical failures (e.g. 'chr1' vs '1').
            with pysam.AlignmentFile(bams[0], "rb") as first:
                unknown = missing_contigs(targets, first.references)
            if unknown:
                raise ValueError(f"BED contig(s) not found in the BAM: {unknown[:5]} "
                                 "(different 'chr' prefix between BED and BAM?)")
        except (OSError, ValueError) as exc:
            log.error("Invalid BED: %s", exc)
            return 2
        n_bp = sum(e - st for t in targets for ivs in t["intervals"].values() for st, e in ivs)
        log.info("Amplicon/panel mode: %d target(s), %d bp (%s)", len(targets), n_bp,
                 "same-name BED lines merged" if merge else "one target per BED line")
    else:
        log.info("Whole-genome mode")
    if interp:
        log.info("Interpretability rule: INTERPRETABLE if at least %g%% of the region "
                 "is covered at >= %dX", interp[1], interp[0])

    # --- Parallelism: N BAM files at a time, each with T threads ---
    jobs = max(1, min(args.max_cores // args.threads_per_bam, len(bams)))
    log.info("%d BAM file(s) found | %d in parallel x %d thread(s)",
             len(bams), jobs, args.threads_per_bam)

    # Analysis parameters sent to every worker
    kwargs = dict(depths=depths, threads=args.threads_per_bam,
                  min_mapq=args.min_mapq, exclude_flags=args.exclude_flags,
                  detail=args.report,       # extra data only needed by the report
                  uniformity=args.uniformity)
    read_stats = args.read_stats or (args.report and targets is None)
    if args.report and targets is None and not args.read_stats:
        log.info("--report in whole-genome mode: read statistics (duplicates...) enabled "
                 "for the duplicates section")
    if targets is None:
        kwargs["read_stats"] = read_stats
    grid_xs = make_grid(list(depths) + ([interp[0]] if interp else [])) if (args.report and interp) else None
    if grid_xs:                              # the report sliders need the exact % at every grid depth
        kwargs["grid_xs"] = grid_xs
    if targets is None:
        kwargs["contigs"] = args.contigs     # --contigs: whole-genome mode only
    tasks = [(b, sample_name(b, args.pattern), kwargs, targets, interp) for b in bams]

    # --- Execution; the parent is the only one writing progress messages ---
    results, details, failed = {}, {}, []    # results: {sample: [rows]}
    for done, (sample, rows, extra, error) in enumerate(_run(tasks, jobs), start=1):
        if error:
            failed.append(sample)
            log.error("[%d/%d] FAILED %s: %s", done, len(bams), sample, error)
        else:
            results[sample] = rows
            details[sample] = extra
            detail = (f"{len(rows)} target(s)" if targets is not None
                      else f"mean coverage {rows[0]['mean_coverage']:.2f}X")
            if interp:                       # verdict summary
                if targets is not None:
                    n_ok = sum(r["interpretability"] == INTERPRETABLE for r in rows)
                    detail += f", {n_ok}/{len(rows)} interpretable"
                else:
                    detail += f", {rows[0]['interpretability']}"
            log.info("[%d/%d] OK %s (%s)", done, len(bams), sample, detail)

    # --- Writing (samples sorted: deterministic order) ---
    if not results:
        log.error("No results to write.")
        return 1
    columns = output_columns(depths, targeted=targets is not None,
                             interpretability=interp is not None,
                             uniformity=args.uniformity, read_stats=read_stats)
    fmt = args.out_format
    samples = sorted(results)
    os.makedirs(args.out_dir, exist_ok=True)

    # Concatenated table: all rows of all samples
    all_rows = [row for sample in samples for row in results[sample]]
    out_path = os.path.join(args.out_dir, f"{args.out_name}.{fmt}")
    write_table(all_rows, columns, out_path, fmt)
    log.info("Main table (%d row(s)): %s", len(all_rows), out_path)

    # Amplicon mode: also one file per sample (1 row per gene)
    if targets is not None:
        per_dir = os.path.join(args.out_dir, "per_sample")
        os.makedirs(per_dir, exist_ok=True)
        for sample in samples:
            write_table(results[sample], columns,
                        os.path.join(per_dir, f"{sample}.amplicons.{fmt}"), fmt)
        log.info("%d per-sample file(s) in %s", len(samples), per_dir)

    # Optional HTML report. Sliders data: exact % covered at every grid depth, one list per
    # table row (same order as all_rows), taken out of the per-sample details.
    grid = None
    if grid_xs:
        grid_rows = []
        for sample in samples:
            g = (details.get(sample) or {}).pop("grid", None)
            grid_rows.extend([g] if targets is None else (g or []))
        if len(grid_rows) == len(all_rows):
            grid = {"x": grid_xs, "rows": grid_rows}
    if args.report:
        report_path = os.path.join(args.out_dir, f"{args.out_name}.html")
        try:
            write_report(
                report_path, title=args.report_title,
                mode="amplicon" if targets is not None else "genome",
                columns=columns, rows=all_rows, samples=samples, details=details,
                depths=depths,
                interp={"depth": interp[0], "pct": interp[1]} if interp else None,
                grid=grid,
                params={k: v for k, v in {
                    "pattern": args.pattern, "min_mapq": args.min_mapq,
                    "exclude_flags": args.exclude_flags,
                    "read_stats": True if read_stats else None,
                    "uniformity": True if args.uniformity else None,
                    "bed": os.path.basename(args.bed) if args.bed else None,
                    "same_gene": (args.same_gene or "merge") if args.bed else None,
                    "contigs": ",".join(args.contigs) if args.contigs else None}.items()
                    if v is not None},
                command="fastbamqc " + " ".join(shlex.quote(a) for a in argv_list),
                failed=failed)
            log.info("HTML report: %s", report_path)
        except Exception as exc:             # noqa: BLE001: never lose the tables for a report error
            log.error("Could not write the HTML report: %s: %s", type(exc).__name__, exc)
            failed.append("<report>")

    if failed:                               # non-zero exit code if at least one failure
        log.warning("%d sample(s) failed: %s", len(failed), ", ".join(sorted(failed)))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
