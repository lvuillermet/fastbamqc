# Changelog

All notable changes to fastbamqc are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[Semantic Versioning](https://semver.org/).

## [1.0.3]
### Changed
- HTML report file size: the data embedded in the page are now stored as integer differences and deflate-compressed
  (decoded by a small decoder inside the page, which stays self-contained and works offline). Nothing is lost: the
  charts and the verdicts are identical. On the synthetic test sets, amplicon mode with 96 samples x 15 targets went
  from 10.6 MB to 1.4 MB, 1,920 samples x 30 targets from 63 MB to 8.2 MB (opening in about 1.3 s), and whole-genome
  mode with 300 samples from 0.62 MB to 0.23 MB.
- HTML report, targets x samples heatmap: when there are too many samples to fit the page and the grid has to be
  scrolled sideways, the target names now stay fixed on the left.
- HTML report, "Mean coverage per target" chart: new "Log scale" checkbox, useful when one or a few targets are
  much deeper than the others.


## [1.0.2a]
### Changed
- Incorporated the html template after test of the version


## [1.0.2]
### Fixed
- An unreadable BAM file (empty, truncated, not a BAM) found first in the input directory no longer produces the
  misleading message "Invalid BED: file does not contain alignment data": the BED is checked against the first
  readable BAM, and the unreadable files are reported as `FAILED` without stopping the other samples.

### Changed
- HTML report on large data sets: tables are paginated (100 rows per page by default, sort and filter act on all
  the rows, the CSV download contains all the filtered rows) and the targets x samples heatmap is drawn on a single
  canvas. A report with 57,600 rows (1,920 samples x 30 targets) now opens in about 1.6 s on the test machine,
  with about 10,000 page elements instead of an estimated one million.


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