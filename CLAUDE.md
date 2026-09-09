# Loupe2Py - Claude Code Guide

## Project Overview

Loupe2Py imports 10x Genomics `.cloupe` files (including Visium HD, binned and cell-segmentation modes) into squidpy-ready `AnnData` objects. All `.cloupe`-format parsing lives in [`cloupe_extract`](cloupe_extract/), a separate Python package that lives inside this same repo but is independently pip-installable — `Loupe2R` (the Seurat sibling package) depends on `cloupe_extract` directly, via `reticulate`, without going through `loupe2py` at all. `loupe2py` itself has no format-specific parsing code of its own; `src/loupe2py/anndata_io.py` builds an `AnnData` object from the generic files `cloupe_extract.extract_cloupe()` writes.

## Key Files

| File | Purpose |
|---|---|
| `src/loupe2py/anndata_io.py` | Builds the `AnnData` object from `cloupe_extract.extract_cloupe()`'s output files — this package's own, novel logic. |
| `src/loupe2py/__init__.py` | Re-exports `cloupe_to_anndata` (own) and `extract_cloupe` (from `cloupe_extract`, for convenience/backward compatibility) as the public API. |
| `tests/test_anndata_io_helpers.py` | Unit tests for `anndata_io.py`'s pure-logic helpers. |
| `tests/test_integration_visium_hd.py` / `test_integration_cellseg.py` | Opt-in regression tests against real paired `.cloupe`/SpaceRanger `outs/` data, gated by env vars (`LOUPE2PY_TEST_DIR`) since the fixtures live outside the repo. Exercise the full `loupe2py` → `cloupe_extract` path end to end. |
| `cloupe_extract/` | The extraction core itself — see [its own `CLAUDE.md`](cloupe_extract/CLAUDE.md). |

## When a new SpaceRanger version is released

Because all `.cloupe`-format parsing lives in `cloupe_extract`, check [its `CLAUDE.md`](cloupe_extract/CLAUDE.md) first — that's where the format-version guard, the `Matrices`-section synthetic-row handling, and the SpaceRanger-output-column checks belong.

Separately, re-run `test_integration_visium_hd.py`/`test_integration_cellseg.py` against a real sample from the new version: they read official SpaceRanger output directly for cross-checking, so a new SpaceRanger version adding, renaming, or restructuring those columns could silently make the comparisons less thorough than they look, even if the tests still pass without error.

## Future work

**Python/AnnData equivalent of real cell-segmentation polygon boundaries, not just centroids.** `cloupe_extract`'s `get_cellseg_projection()` computes each cell's true area centroid from `CellSegs.GeoJSON` — exact, but still one point per cell, discarding the actual shape. Plain `AnnData` has no native polygon-boundary slot; the scverse-idiomatic answer is `spatialdata` (`ShapesModel`, GeoPandas-backed) — the same package this project's own README already points Visium HD users toward for anything beyond single-resolution `AnnData`. Producing real polygon output from `Loupe2Py` properly likely means an additional `spatialdata`-object output path, not bolting polygons onto `AnnData` in some ad hoc, non-idiomatic way. This is a separate, larger decision from the Seurat-side equivalent (see `Loupe2R`'s `CLAUDE.md` and `cloupe_extract`'s `CLAUDE.md`) — not a smaller version of it. Not scoped in detail yet.
