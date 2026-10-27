# Changelog

All notable changes to fastbamqc are documented in this file.
The format follows Keep a Changelog and the project uses Semantic Versioning.

## [1.0.1] - 2026-10-27

### Fixed 
- Syntax error in `pyproject.toml` (unclosed quote) that prevent building the package in 1.0.0. 
## [1.0.0] - 2026-10-27

First public release.

### Features
- Whole-genome mode: one row per sample (total_read, mapped_read, percentage_mapped,
  mean_coverage, <d>X = % of positions covered at ≥ d reads).
- Amplicon / panel mode (--bed): one row per target and per sample, with --same-gene {merge,separate}.
- Exact per-position depth computed over independent windows (start/end difference array).
- Interpretability verdict on by default (75 % of the region at ≥ 10X); --interp-depth,
  --interp-pct, --no-interp, and --interp-nodup (rule evaluated without duplicates, whole genome only).
- Self-contained HTML report (--report) with live sliders for the interpretability rule and,
  for whole genomes, a duplicates section.
- Optional extra metrics: --uniformity and --read-stats (read-level QC, flagged and estimated duplicates).
- Multi-core processing: -c (total cores) and -t (cores per BAM); a single BAM is split into windows
  shared between processes, with results identical whatever the number of cores.
- Output formats: TSV and CSV.

### Notes
- Versions before 1.0.0 were internal development versions and are not documented here.