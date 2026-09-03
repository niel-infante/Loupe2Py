# Loupe2Py

Import 10x Genomics `.cloupe` files (including Visium HD) into [squidpy](https://squidpy.readthedocs.io/)-ready [AnnData](https://anndata.readthedocs.io/) objects.

R/Seurat user? See the companion package, [**Loupe2R**](https://github.com/niel-infante/Loupe2R). Looking for the full function list (including internal, unexported functions)? See [`REFERENCE.md`](REFERENCE.md).

squidpy has no object type of its own — it operates directly on `AnnData`, following the same `.obsm['spatial']` / `.uns['spatial'][library_id]` convention [scanpy](https://scanpy.readthedocs.io/) established. `Loupe2Py` extracts the count matrix, spatial coordinates, tissue image, UMAP embedding, and Space Ranger cluster labels from a `.cloupe` file and assembles an `AnnData` object in that convention, ready for squidpy's spatial analysis functions directly.

This is the Python sibling of [Loupe2R](https://github.com/niel-infante/Loupe2R), which does the same job for Seurat. Both are built on the same underlying `.cloupe`-parsing core.

**Read [Limitations & Risks](#limitations--risks) before using this on anything you plan to publish.** This package parses an undocumented, proprietary file format using an unofficial, reverse-engineered parser (vendored — see [Credits](#credits)). It has been validated carefully on real data (see [Validation](#validation) below), but it is not a substitute for 10x Genomics' own SpaceRanger output when that's available.

## Installation

```bash
pip install git+https://github.com/niel-infante/Loupe2Py.git
```

No separate install step for the `.cloupe` parser — unlike Loupe2R, `Loupe2Py` vendors a pinned copy of it directly (see [Credits](#credits)), so there's no external path to configure.

## Quick start

```python
from loupe2py import cloupe_to_anndata

adata = cloupe_to_anndata("path/to/sample.cloupe")

import squidpy as sq
sq.gr.spatial_neighbors(adata)
sq.pl.spatial_scatter(adata, color="total_counts")
```

## What this does and does not do

- Extracts counts, spatial coordinates, the tissue image, UMAP/other embeddings (`obsm['X_<name>']`), Space Ranger clusterings (`obs['sr_<name>']`), and any user-created Loupe Browser cell tracks.
- Captures **one** bin resolution per `.cloupe` file — see the multi-resolution note below.
- Reflects whatever default, unfiltered pipeline output Loupe Browser was showing — not your own QC, filtering, or normalization choices.

## Limitations & Risks

**The `.cloupe` format is proprietary and undocumented.** It's a custom binary container (JSON header + byte-offset index, not SQLite) with no public spec. The vendored parser this package is built on, [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe), is an unofficial, reverse-engineered tool with no version-compatibility guarantees across Loupe Browser/CellRanger/SpaceRanger releases. `Loupe2Py` checks the file's internal format-version fields against a known-tested set on every extraction and warns (via Python's `warnings`) — while still proceeding — if a file reports an unrecognized version. The detected versions are always available afterward via `adata.uns['cloupe_format_info']`.

**Loupe-derived annotations reflect default, unfiltered pipeline output, not your own analysis.** In [satijalab/seurat#9269](https://github.com/satijalab/seurat/issues/9269), a Seurat maintainer cautions that data exported from Loupe Browser reflects only the default parameters SpaceRanger/CellRanger ran with — no manual QC, filtering, or normalization decisions are captured. This applies regardless of which downstream framework you're targeting. Treat `sr_*` clusterings as exploratory defaults, not a substitute for your own analysis.

**A `.cloupe` file cannot represent Visium HD's full multi-resolution structure — and squidpy itself has no native answer for that either.** squidpy's own maintainers redirect Visium HD users to a different scverse package, [`spatialdata`](https://spatialdata.scverse.org/), for genuine multi-resolution work; plain `AnnData` (what this package and `squidpy.read.visium()` both produce) only ever represents one bin resolution. `cloupe_to_anndata()` matches that same scope — one `.cloupe` file, one resolution. If you need multiple resolutions combined, that's a `spatialdata`-shaped problem this package doesn't attempt to solve.

**Visium HD's cell-segmentation output uses a coarser spatial fallback than its binned output.** SpaceRanger 4.0+ can process a Visium HD run in two ways — binned (square bins at a chosen resolution) or cell-segmentation mode (transcripts assigned to individual cells via image-based nucleus/cell segmentation) — each producing its own separate `.cloupe` file; cell segmentation is a Visium HD analysis mode, not a separate non-HD product. Binned-mode barcodes encode their own exact grid position (see [Validation](#validation)); cell-segmentation-mode barcodes don't use that convention, so their coordinates fall back to averaging the file's per-cell-segment centroids for each barcode, which is inherently less precise than a single measured point but not corrupted — a real bug in this averaging step (Centers read with the wrong array layout, producing a scrambled diagonal artifact for 30-40% of cells) was found and fixed in 0.1.1; see `NEWS.md`. `array_row`/`array_col` are unconditionally `(0, 0)` for cell-segmentation-mode `.cloupe` files — there's no fixed grid for individual cells to round to — which doesn't affect plotting (that uses the real pixel coordinates) but is worth knowing if you read those columns directly.

**Prefer official SpaceRanger output when you have it — but note `squidpy.read.visium()` itself can't read modern Visium HD directories.** Confirmed directly against a real sample: `squidpy.read.visium()` (as of squidpy 1.8.2) only looks for `spatial/tissue_positions.csv` (or the older `tissue_positions_list.csv`), but current SpaceRanger HD output writes `spatial/tissue_positions.parquet` instead — `squidpy.read.visium()` raises `FileNotFoundError` outright on a real, current HD `outs/` directory. This is consistent with squidpy's own maintainers pointing HD users at `spatialdata` instead (see above), and is why this package's own validation (below) reads the official `.h5`/`.parquet` files directly rather than through `squidpy.read.visium()`.

## Credits

The hard part — reverse-engineering the proprietary `.cloupe` binary format at all (the header layout, byte-offset index block, matrix/projection encoding, and tiled image storage) — is not this package's work. It's [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe), written by **Martin Prete** and **Nithin Mathew Joseph** of the Wellcome Sanger Institute's Cellular Genetics Informatics (cellgeni) team, AGPL-3.0 licensed. A pinned copy is vendored directly into this package at `src/loupe2py/_vendor/cloupe.py` (see that file's header for the exact commit) rather than kept as an external dependency, to eliminate a class of path-configuration bugs and pin an exact, tested parser version. Full license text: [`THIRD_PARTY_LICENSES/cellgeni-cloupe-AGPL-3.0.txt`](THIRD_PARTY_LICENSES/cellgeni-cloupe-AGPL-3.0.txt).

**Note:** vendoring this AGPL-3.0 code makes this package a combined work under AGPL-3.0 terms, so `Loupe2Py` is itself licensed **AGPL-3.0-or-later** — see [License](#license) below.

Also indebted to the [scverse](https://scverse.org/) ecosystem ([AnnData](https://anndata.readthedocs.io/), [scanpy](https://scanpy.readthedocs.io/), [squidpy](https://squidpy.readthedocs.io/)) for the conventions and reference outputs this package's Visium HD support was validated against, and to [10x Genomics' SpaceRanger](https://www.10xgenomics.com/support/software/space-ranger) for the official output used as ground truth (see [Validation](#validation)).

## Validation

`cloupe_to_anndata()` was cross-validated against the official SpaceRanger output for the same Visium HD sample used to validate Loupe2R (8µm bins, 19,072 features × 472,752 barcodes). Results:

| Check | Result |
|---|---|
| Barcode overlap (cloupe vs. official filtered matrix) | 100% |
| Count correlation on shared barcodes | 1.000000 |
| `array_row`/`array_col` exact match vs. official `tissue_positions.parquet` | 100% (0 max difference) |
| `pxl_col`/`pxl_row` correlation | 1.000000 |

Identical results to Loupe2R's own validation on the same sample, as expected — same extraction logic underneath. Also smoke-tested against real squidpy functions (`squidpy.gr.spatial_neighbors()`, `squidpy.pl.spatial_scatter()`), not just hand-written assertions.

The `array_row`/`array_col` result relies on the same fix validated in Loupe2R: Visium HD barcodes encode their exact grid position directly (e.g. `s_008um_00269_00526-1` → row 269, column 526), and `loupe2py.extract.parse_array_position()` parses that directly rather than approximating it from pixel coordinates (which was found to disagree with the official grid by hundreds of bins due to a pixel-origin mismatch between the stitched tissue image and SpaceRanger's own coordinate frame). Non-HD barcodes fall back to the pixel-based approximation.

Reproduce this validation yourself (see `tests/test_integration_visium_hd.py`) by setting `LOUPE2PY_TEST_DIR` to a SpaceRanger Visium HD `outs/` directory containing `cloupe_008um.cloupe` and `binned_outputs/square_008um/`, and running `pytest` with the `test` extras installed.

## License

AGPL-3.0-or-later (see `LICENSE`). This package vendors AGPL-3.0 code ([Credits](#credits) above) directly into its own source, which under AGPL-3.0's terms makes the combined work AGPL-3.0 as well. `THIRD_PARTY_LICENSES/` contains the full license text for the vendored code specifically.
