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
                interpretability rule, run parameters...
  detail-<i>    extra data of sample i (depth histogram and per-contig metrics
                in whole-genome mode; depth profiles in amplicon mode). It is
                parsed by the page only when needed.
"""
from __future__ import annotations

import html
import json
import re
from datetime import datetime
from importlib import resources

from . import __version__


def _json_for_script(obj):
    """JSON safe to embed inside <script>: '<' is escaped so '</script>' cannot appear."""
    def default(o):                           # numpy scalars -> plain Python numbers
        return o.item() if hasattr(o, "item") else str(o)
    return json.dumps(obj, separators=(",", ":"), default=default).replace("<", "\\u003c")


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
    data = {
        "tool": "fastbamqc", "version": __version__, "title": title, "mode": mode,
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "command": command, "columns": columns, "depths": list(depths),
        "interp": interp, "grid": grid, "params": params, "samples": list(samples),
        "failed": sorted(failed),
        "rows": [[row.get(c) for c in columns] for row in rows],   # compact: arrays, not dicts
    }
    detail_tags = "\n".join(
        f'<script id="detail-{i}" type="application/json">{_json_for_script(details.get(s))}</script>'
        for i, s in enumerate(samples))

    template = resources.files("fastbamqc").joinpath("templates/report.html").read_text(encoding="utf-8")
    values = {"TITLE": html.escape(title), "DATA": _json_for_script(data), "DETAILS": detail_tags}
    # one single pass: inserted text is never re-scanned for placeholders
    page = re.sub(r"__(TITLE|DATA|DETAILS)__", lambda m: values[m.group(1)], template)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(page)
