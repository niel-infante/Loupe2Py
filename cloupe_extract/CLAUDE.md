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

## What's vendored vs. novel, per extracted item

`cellgeni/cloupe` (vendored in `_vendor/cloupe.py`) is a low-level binary-container reader — it decodes the header, byte-offset index, and raw section contents generically. Everything about *which* sections matter for Visium HD, how to interpret them, and how to turn them into usable output is this project's own work. Per item actually extracted:

| Extracted item | `.cloupe` section(s) | `cellgeni/cloupe` provides | Novel (this project) |
|---|---|---|---|
| Count matrix (barcodes × features, sparse) | `Matrices` | **Yes** — full parse: `Barcodes`, `FeatureIds`, `FeatureNames`, `UMICounts`, CSR/CSC arrays | Feature×barcode → obs×var transpose; MTX/AnnData serialization; `exclude_synthetic_totals()` |
| Barcode / feature ID and name lists | `Matrices` | **Yes** — same section | Re-serialization as `barcodes.tsv.gz`/`features.tsv.gz` (10x convention) |
| UMAP / other non-spatial embeddings | `Projections` | **Yes**, generically | Filtering `Spatial`/`Fiducials` out; writing the rest as `obsm['X_<name>']`/Seurat `DimReducObject`s |
| User-created Loupe Browser cell tracks | `CellTracks` | **Yes** — full parse | Re-serialized as `celltracks.csv`/`obs` columns |
| Visium HD spatial pixel coordinates (bin-level) | `Projections` (`Spatial`) | **Partial** — same generic decode as UMAP | Recognizing `Spatial` specifically; 3rd array as `bin_size_px`; `CellSegs` fallback |
| Visium HD array row/column (grid position) | *(barcode string)* | No | `parse_array_position()` parses it directly from the barcode string; fixes a disagreement (up to ~1,400 bins) from the earlier pixel-rounding approach |
| Cell-segmentation spatial coordinates | `CellSegs` | **No** — section untouched by their code | Exact polygon-boundary area centroid via `CellSegs.GeoJSON` (`loupe2py` 0.2.2, pre-split); also fixed a diagonal-scrambling `Centers` layout bug (`loupe2py` 0.1.1, pre-split) |
| Tissue image (full-resolution, stitched) | `SpatialImageTiles` | **No** — absent from shipped code | `stitch_tiles()`: tile-pyramid zoom selection, positioning/cropping, RGB normalization |
| Space Ranger graph/k-means clusterings | `Clusterings` | **No** — never read by their code | `read_clusterings()`: int16 group-key decode against a groups metadata list |
| Scale factors (spot diameter, hires/lowres, bin size) | *(derived)* | No | Derived from `bin_size_px`/microns-per-pixel/image dimensions |
| Format-version compatibility check | header, per-section `FormatVersion` | **Partial** — fields are present, but nothing checks them | `check_format_version()` against `_TESTED_VERSIONS`, enforced by default (raises `UnvalidatedFormatVersionError`, see below) |
| Seurat/AnnData object assembly | *(assembly)* | **No** — their `to_anndata()` is bare, no spatial/image convention | `cloupe_to_anndata()` (`loupe2py`), `cloupe_to_seurat()` (`Loupe2R`) — downstream of this package, listed here for completeness |
| Multi-resolution combination | *(assembly)* | N/A | `combine_cloupe_bins()` (`Loupe2R` only) — downstream of this package, listed here for completeness |

This was originally a table in the JBT paper draft (`loupe2r.md`); it was cut for length and moved here, since this is where it stays accurate as the code changes.

## Format-version guard

`check_format_version()` compares a file's internal `container`/`run`/`matrix`/`projection` format-version fields against `_TESTED_VERSIONS` in `extract.py` and returns the detected versions plus any warnings — it never raises itself, it's a pure inspection function. `extract_cloupe()` is what acts on the result: by default (`version_check=True`) it raises `UnvalidatedFormatVersionError` before extracting anything if any warning was produced, since the parser may well still be correct but that hasn't been checked, and a data-recovery tool guessing silently is the wrong default. `version_check=False` downgrades this back to a printed warning and lets extraction proceed — the caller is then explicitly responsible for verifying the result. `format_info.json` (detected versions + warnings) is written either way, including on the raising path, so provenance survives a refusal. Extend `_TESTED_VERSIONS` whenever a new combination is actually validated against real paired SpaceRanger output, not just assumed compatible.

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
