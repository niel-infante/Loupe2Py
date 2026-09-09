# cloupe_extract - Claude Code Guide

## Project Overview

`cloupe_extract` is the shared `.cloupe`-parsing core underlying both `loupe2py` (this repo's Python/AnnData package, one directory up) and [`Loupe2R`](https://github.com/niel-infante/Loupe2R) (the sibling R/Seurat package, via `reticulate`). It is framework-agnostic: it parses a `.cloupe` file and writes generic, SpaceRanger-convention-like intermediate files (`matrix.mtx`, `barcodes.tsv`, `features.tsv`, `tissue_positions.csv`, `scalefactors_json.json`, an optional tissue image, and optional `cell_boundaries`/clustering outputs) — it has no knowledge of `AnnData` or `Seurat` themselves.

`.cloupe` is a proprietary, undocumented binary format. Parsing is built on a vendored, pinned copy of [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe) (Wellcome Sanger Institute) at `src/cloupe_extract/_vendor/cloupe.py` — a low-level, reverse-engineered binary-container reader. `src/cloupe_extract/extract.py` adds everything format-specific on top of that: Visium HD bin-position recovery, cell-segmentation coordinate handling, tissue-image tile-pyramid reconstruction, Space Ranger clustering decode, and object assembly.

This package was split out of `loupe2py` into its own installable package so that `Loupe2R` could depend on it directly by name, instead of installing a package literally called `loupe2py` to get Seurat support — see `NEWS.md` in both `loupe2py` and `Loupe2R` for the split itself.

## Key Files

| File | Purpose |
|---|---|
| `src/cloupe_extract/_vendor/cloupe.py` | Vendored `cellgeni/cloupe` binary parser — do not modify directly; if a fix is needed, either patch locally with a clear comment explaining the divergence from upstream, or fix upstream and re-vendor at a new pinned commit. |
| `src/cloupe_extract/extract.py` | All extraction logic novel to this project: `extract_cloupe()` (main entry point), `get_spatial_projection()` / `get_cellseg_projection()` (bin/cell-seg position recovery), `parse_array_position()` (barcode-encoded grid position), `stitch_tiles()` (tissue image reconstruction), `read_clusterings()`, `check_format_version()`, `exclude_synthetic_totals()`. |
| `tests/test_extract_helpers.py` | Unit tests for the pure-logic helpers in `extract.py` — lightweight fakes, no real `.cloupe` file needed. |

Real-data integration tests (`test_integration_visium_hd.py`, `test_integration_cellseg.py`) live in `loupe2py`'s own `tests/` directory one level up, not here — they exercise the full `loupe2py` → `cloupe_extract` path end to end, gated by env vars since the fixtures live outside the repo.

## Format-version guard

`check_format_version()` compares a file's internal `container`/`run`/`matrix`/`projection` format-version fields against `_TESTED_VERSIONS` in `extract.py`. An unrecognized version triggers a warning, not a failure — the parser may well still be correct, it just hasn't been checked. Extend `_TESTED_VERSIONS` whenever a new combination is actually validated against real paired SpaceRanger output, not just assumed compatible.

## When a new SpaceRanger version is released

Re-validate against a real sample from the new version before assuming compatibility, and specifically check for:

- **New or renamed columns** in `tissue_positions.parquet`/`.csv`, `scalefactors_json.json`, and `metrics_summary.csv`. SpaceRanger's own output conventions have changed across versions before (e.g. `.csv` → `.parquet` for tissue positions in the HD era), and this package's own written output (`tissue_positions.csv`, `scalefactors_json.json`) mirrors whatever convention it was built against — not a live spec that tracks SpaceRanger automatically.
- **New or renamed rows in the `.cloupe` `Matrices` section.** Confirmed (2026-09) that every `.cloupe` file appends synthetic `type_sum_<FeatureType>`/`genome_sum_<Reference>` aggregate rows after the real genes, which were silently written out as if they were genes until `exclude_synthetic_totals()` was added. A future SpaceRanger/Loupe Browser version could plausibly add a *different* kind of synthetic row under a similar convention; check for this specifically rather than assuming the current prefix match (`type_sum_` / `genome_sum_`) covers everything a newer file might add.
- **New optional `.cloupe` sections** (beyond `Matrices`/`Projections`/`CellSegs`/`Clusterings`/`Analyses`/`Runs`/`Metrics`) that a newer release might start populating.

Add a new version to `_TESTED_VERSIONS` only after directly confirming exact-match concordance against that version's real paired SpaceRanger output (barcode overlap, per-barcode count values, array position, pixel coordinates) — not just that extraction runs without error.

## Future work

**Expose full cell-segmentation polygon boundaries, not just centroids.** `get_cellseg_projection()` reads `CellSegs.GeoJSON` to compute each cell's true area centroid — exact, but still collapses every cell down to one point. The full polygon boundary is already parsed in memory at that point and then discarded. Both downstream packages have a native, purpose-built way to consume real cell shapes instead of a single point — Seurat's `FOV`/`Segmentation` object (see `Loupe2R`'s own `CLAUDE.md`) and, on the Python side, `spatialdata`'s `ShapesModel`. Scoped but not started:

- `extract_cloupe()` would need a new opt-in output (e.g. `include_boundaries=True`, matching the `include_image=True` pattern) writing a long-format `cell_boundaries.csv`: `barcode, x, y`, one row per polygon vertex — a direct flatten of the same `GeoJSON` rings `get_cellseg_projection()` already reads. Real-data scale: ~148k cells × ~15–20 vertices each on the Human Kidney FFPE dataset, so a few million rows — comparable to `tissue_positions.csv` for HD binned data, not free but tractable.
- Only meaningful for cell-segmentation-mode files (binned HD has no `CellSegs`).
- The rest of this work is downstream-specific (Seurat object construction on the `Loupe2R` side, an `AnnData`/`spatialdata` decision on the `loupe2py` side) — see each package's own `CLAUDE.md`, not this one.
