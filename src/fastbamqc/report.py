"""
fastbamqc.report
================
Writes the self-contained HTML report (--report).

The page is a single file: the HTML/CSS/JS template (templates/report.html) with
the results embedded as JSON. It has no external dependency, so it opens offline
and can be e-mailed. All charts are drawn by the page itself (plain SVG).

Embedded data
-------------
  report-data   general data: the rows of the main table, depth thresholds,
                interpretability rule, run parameters, slider grid...
  detail-<i>    extra data of sample i (depth histogram and per-contig metrics
                in whole-genome mode; depth profiles in amplicon mode). It is
                decoded by the page only when needed.

Compact payload (keeps the file small)
--------------------------------------
Every block is base64( uint32 raw size + deflate( uint32 JSON size + JSON + integers ) ).
The integers (depth profiles, slider grid) are stored outside the JSON as
zigzag-varint differences between neighbours; the page decodes everything with a
small inflate written in the template, so the file stays self-contained.
"""
from __future__ import annotations

import base64
import html
import json
import re
import struct
import zlib
from datetime import datetime
from importlib import resources

from . import __version__


def _varints(values):
    """Integers -> bytes: difference with the previous value, zigzag, then LEB128 varint."""
    out, prev = bytearray(), 0
    for v in values:
        d = int(v) - prev
        prev = int(v)
        z = (d << 1) ^ (d >> 63)
        while z >= 128:
            out.append((z & 127) | 128)
            z >>= 7
        out.append(z)
    return bytes(out)


def _pack(obj, int_runs=()):
    """JSON `obj` + integer runs (each one delta-coded on its own) -> deflated base64 text."""
    def default(o):                           # numpy scalars -> plain Python numbers
        return o.item() if hasattr(o, "item") else str(o)
    js = json.dumps(obj, separators=(",", ":"), default=default).encode("utf-8")
    raw = struct.pack("<I", len(js)) + js + b"".join(_varints(run) for run in int_runs)
    co = zlib.compressobj(9, zlib.DEFLATED, -15)
    return base64.b64encode(struct.pack("<I", len(raw)) + co.compress(raw) + co.flush()).decode("ascii")


def _pack_detail(detail):
    """Per-sample details: depth profiles (0.1 resolution) go to the integer part."""
    if not detail or not detail.get("profiles"):
        return _pack(detail)
    detail = dict(detail)
    runs, slim = [], []
    for p in detail["profiles"]:
        runs.append([int(round(v * 10)) for v in p["depth"]])
        slim.append({k: v for k, v in p.items() if k != "depth"} | {"n": len(p["depth"])})
    detail["profiles"] = slim
    return _pack(detail, runs)


def write_report(path, *, title, mode, columns, rows, samples, details, depths,
                 interp, params, command, failed=(), grid=None):
    """
    path     output HTML file
    mode     "genome" or "amplicon"
    columns  column names of the main table; rows = list of dicts using these keys
    samples  sample names (sorted); details = {sample: extra data or None}
    interp   {"depth": int, "pct": float} or None
    grid     {"x": [depths], "rows": [[% x 100 at each depth] per row of `rows`]} or None;
             with `interp`, it lets the page recompute the verdicts when the sliders move
    params   {name: value} shown in the page footer
    command  the command line, for reproducibility
    """
    grid_rows = None
    if grid and grid.get("rows"):
        grid_rows = grid["rows"]
        grid = {"x": list(grid["x"]), "n": len(grid_rows)}
    data = {
        "tool": "fastbamqc", "version": __version__, "title": title, "mode": mode,
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "command": command, "columns": columns, "depths": list(depths),
        "interp": interp, "grid": grid, "params": params, "samples": list(samples),
        "failed": sorted(failed),
        "rows": [[row.get(c) for c in columns] for row in rows],   # compact: arrays, not dicts
    }
    detail_tags = "\n".join(
        f'<script id="detail-{i}" type="text/plain">{_pack_detail(details.get(s))}</script>'
        for i, s in enumerate(samples))

    template = resources.files("fastbamqc").joinpath("templates/report.html").read_text(encoding="utf-8")
    values = {"TITLE": html.escape(title), "DATA": _pack(data, grid_rows or ()), "DETAILS": detail_tags}
    # one single pass: inserted text is never re-scanned for placeholders
    page = re.sub(r"__(TITLE|DATA|DETAILS)__", lambda m: values[m.group(1)], template)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(page)
