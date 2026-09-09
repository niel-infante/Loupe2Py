# Loupe2Py - Claude Code Guide

## Project Overview

Loupe2Py imports 10x Genomics `.cloupe` files (including Visium HD, binned and cell-segmentation modes) into squidpy-ready `AnnData` objects. It's the Python sibling of `Loupe2R` (Seurat), which calls this package's extraction core via `reticulate`.

`.cloupe` is a proprietary, undocumented binary format. Parsing is built on a vendored, pinned copy of [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe) (Wellcome Sanger Institute) in `src/loupe2py/_vendor/cloupe.py` — a low-level, reverse-engineered binary-container reader. `src/loupe2py/extract.py` adds everything format-specific on top of that: Visium HD bin-position recovery, cell-segmentation coordinate handling, tissue-image tile-pyramid reconstruction, Space Ranger clustering decode, and object assembly.

## Key Files

| File | Purpose |
|---|---|
| `src/loupe2py/_vendor/cloupe.py` | Vendored `cellgeni/cloupe` binary parser — do not modify directly; if a fix is needed, either patch locally with a clear comment explaining the divergence from upstream, or fix upstream and re-vendor at a new pinned commit. |
| `src/loupe2py/extract.py` | All extraction logic novel to this project: `extract_cloupe()` (main entry point), `get_spatial_projection()` / `get_cellseg_projection()` (bin/cell-seg position recovery), `parse_array_position()` (barcode-encoded grid position), `stitch_tiles()` (tissue image reconstruction), `read_clusterings()`, `check_format_version()`, `exclude_synthetic_totals()`. |
| `src/loupe2py/anndata_io.py` | Builds the `AnnData` object from `extract_cloupe()`'s output files. |
| `tests/test_extract_helpers.py` | Unit tests for the pure-logic helpers in `extract.py` — lightweight fakes, no real `.cloupe` file needed. |
| `tests/test_integration_visium_hd.py` / `test_integration_cellseg.py` | Opt-in regression tests against real paired `.cloupe`/SpaceRanger `outs/` data, gated by env vars (`LOUPE2PY_TEST_DIR`) since the fixtures live outside the repo. |

## Format-version guard

`check_format_version()` compares a file's internal `container`/`run`/`matrix`/`projection` format-version fields against `_TESTED_VERSIONS` in `extract.py`. An unrecognized version triggers a warning, not a failure — the parser may well still be correct, it just hasn't been checked. Extend `_TESTED_VERSIONS` whenever a new combination is actually validated against real paired SpaceRanger output, not just assumed compatible.

## When a new SpaceRanger version is released

Re-validate against a real sample from the new version before assuming compatibility, and specifically check for:

- **New or renamed columns** in `tissue_positions.parquet`/`.csv`, `scalefactors_json.json`, and `metrics_summary.csv`. SpaceRanger's own output conventions have changed across versions before (e.g. `.csv` → `.parquet` for tissue positions in the HD era), and this package's own written output (`tissue_positions.csv`, `scalefactors_json.json`) mirrors whatever convention it was built against — not a live spec that tracks SpaceRanger automatically.
- **New or renamed rows in the `.cloupe` `Matrices` section.** Confirmed (2026-09) that every `.cloupe` file appends synthetic `type_sum_<FeatureType>`/`genome_sum_<Reference>` aggregate rows after the real genes, which were silently written out as if they were genes until `exclude_synthetic_totals()` was added — see `NEWS.md` 0.2.1. A future SpaceRanger/Loupe Browser version could plausibly add a *different* kind of synthetic row under a similar convention; check for this specifically rather than assuming the current prefix match (`type_sum_` / `genome_sum_`) covers everything a newer file might add.
- **New optional `.cloupe` sections** (beyond `Matrices`/`Projections`/`CellSegs`/`Clusterings`/`Analyses`/`Runs`/`Metrics`) that a newer release might start populating.

Add a new version to `_TESTED_VERSIONS` only after directly confirming exact-match concordance against that version's real paired SpaceRanger output (barcode overlap, per-barcode count values, array position, pixel coordinates) — not just that extraction runs without error.

## Future work

**Expose full cell-segmentation polygon boundaries, not just centroids.** `get_cellseg_projection()` (fixed in 0.2.2) reads `CellSegs.GeoJSON` to compute each cell's true area centroid — exact, but still collapses every cell down to one point. The full polygon boundary is already parsed in memory at that point and then discarded. Seurat has a native, purpose-built object for real cell shapes (`FOV`/`Segmentation`, confirmed via `CreateSegmentation()`'s actual source — see `Loupe2R`'s own `CLAUDE.md` for the exact API), and `Loupe2R::cloupe_to_seurat()` could build one instead of the plain centroid-only spatial assay it does today. Scoped but not started:

- `extract_cloupe()` would need a new opt-in output (e.g. `include_boundaries=True`, matching the `include_image=True` pattern) writing a long-format `cell_boundaries.csv`: `barcode, x, y`, one row per polygon vertex — a direct flatten of the same `GeoJSON` rings `get_cellseg_projection()` already reads. Real-data scale: ~148k cells × ~15–20 vertices each on the Human Kidney FFPE dataset, so a few million rows — comparable to `tissue_positions.csv` for HD binned data, not free but tractable.
- Only meaningful for cell-segmentation-mode files (binned HD has no `CellSegs`).
- The harder half of this work is entirely on the `Loupe2R` side (Seurat object construction), not here — see that repo's `CLAUDE.md`.

**Python/AnnData equivalent is a separate, larger decision, not a smaller version of the above.** Plain `AnnData` has no native polygon-boundary slot. The scverse-idiomatic answer is `spatialdata` (`ShapesModel`, GeoPandas-backed) — the same package this project's own README already points Visium HD users toward for anything beyond single-resolution `AnnData`. Producing real polygon output from `Loupe2Py` properly likely means an additional `spatialdata`-object output path, not bolting polygons onto `AnnData` in some ad hoc, non-idiomatic way. Treat this as its own scoping exercise when it comes up, not an assumed extension of the Seurat work above.
