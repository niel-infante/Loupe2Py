# Loupe2Py

Import 10x Genomics `.cloupe` files (including Visium HD) into [squidpy](https://squidpy.readthedocs.io/)-ready [AnnData](https://anndata.readthedocs.io/) objects.

R/Seurat user? See the companion package, [**Loupe2R**](https://github.com/niel-infante/Loupe2R). Looking for the full function list (including internal, unexported functions)? See [`REFERENCE.md`](REFERENCE.md).

squidpy has no object type of its own — it operates directly on `AnnData`, following the same `.obsm['spatial']` / `.uns['spatial'][library_id]` convention [scanpy](https://scanpy.readthedocs.io/) established. `Loupe2Py` extracts the count matrix, spatial coordinates, tissue image, UMAP embedding, and Space Ranger cluster labels from a `.cloupe` file and assembles an `AnnData` object in that convention, ready for squidpy's spatial analysis functions directly.

This is the Python sibling of [Loupe2R](https://github.com/niel-infante/Loupe2R), which does the same job for Seurat. Both are built on the same underlying `.cloupe`-parsing core.

**Read [Limitations & Risks](#limitations--risks) before using this on anything you plan to publish.** This package parses an undocumented, proprietary file format via [`cloupe_extract`](https://github.com/niel-infante/Loupe2Py/tree/main/cloupe_extract), which wraps an unofficial, reverse-engineered parser (vendored — see [Credits](#credits)). It has been validated carefully on real data (see [Validation](#validation) below), but it is not a substitute for 10x Genomics' own SpaceRanger output when that's available.

## Installation

```bash
pip install git+https://github.com/niel-infante/Loupe2Py.git
```

No separate install step for the `.cloupe` parser — `pip install` above pulls in [`cloupe_extract`](cloupe_extract/) (which vendors the parser itself, see [Credits](#credits)) automatically as a dependency, unlike Loupe2R, where you install `cloupe_extract` yourself since R packages can't declare a Python dependency.

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

**The `.cloupe` format is proprietary and undocumented.** It's a custom binary container (JSON header + byte-offset index, not SQLite) with no public spec. `Loupe2Py` depends on [`cloupe_extract`](https://github.com/niel-infante/Loupe2Py/tree/main/cloupe_extract) for the actual parsing, which vendors [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe), an unofficial, reverse-engineered tool with no version-compatibility guarantees across Loupe Browser/CellRanger/SpaceRanger releases. Every extraction checks the file's internal format-version fields against a known-tested set. By default, a file reporting an unrecognized version raises `cloupe_extract.UnvalidatedFormatVersionError` before anything is extracted — the parser may well still be correct, but that hasn't been checked, so `cloupe_to_anndata()`/`extract_cloupe()` refuse to guess silently. Pass `version_check=False` to proceed anyway (a warning is printed and the file is extracted); you're then responsible for independently verifying the result. When extraction does proceed, the detected versions are available via `adata.uns['cloupe_format_info']`; if it instead raises, the specific unvalidated version(s) are listed directly in the exception message (`format_info.json` is also written to `outdir` first, but won't survive past the exception under the default auto-cleaned temp directory — pass `outdir`/`keep_files=True` if you want it to persist).

**Loupe-derived annotations reflect default, unfiltered pipeline output, not your own analysis.** In [satijalab/seurat#9269](https://github.com/satijalab/seurat/issues/9269), a Seurat maintainer cautions that data exported from Loupe Browser reflects only the default parameters SpaceRanger/CellRanger ran with — no manual QC, filtering, or normalization decisions are captured. This applies regardless of which downstream framework you're targeting. Treat `sr_*` clusterings as exploratory defaults, not a substitute for your own analysis.

**A `.cloupe` file cannot represent Visium HD's full multi-resolution structure — and squidpy itself has no native answer for that either.** squidpy's own maintainers redirect Visium HD users to a different scverse package, [`spatialdata`](https://spatialdata.scverse.org/), for genuine multi-resolution work; plain `AnnData` (what this package and `squidpy.read.visium()` both produce) only ever represents one bin resolution. `cloupe_to_anndata()` matches that same scope — one `.cloupe` file, one resolution. If you need multiple resolutions combined, that's a `spatialdata`-shaped problem this package doesn't attempt to solve.

**Every `.cloupe` file's count matrix carried 2 bogus "genes" with inflated totals until 0.2.1.** The `.cloupe` format itself appends synthetic `type_sum_<FeatureType>`/`genome_sum_<Reference>` rows to the `Matrices` section — each holding a barcode's *entire* total for that category, not a real gene's count. `extract_cloupe()`/`cloupe_to_anndata()` wrote these out as if they were ordinary genes, which inflated `total_counts`/`percent.mt`-style metrics (3x, for the common single-feature-type/single-genome case) and added 2 meaningless entries to `.var`. Fixed in 0.2.1 — see `NEWS.md`. Re-extract if you used an earlier version; the real per-gene counts themselves were always correct, only these 2 extra rows were the problem.

**Visium HD's cell-segmentation output uses the file's exact cell boundary data for position.** SpaceRanger 4.0+ can process a Visium HD run in two ways — binned (square bins at a chosen resolution) or cell-segmentation mode (transcripts assigned to individual cells via image-based nucleus/cell segmentation) — each producing its own separate `.cloupe` file; cell segmentation is a Visium HD analysis mode, not a separate non-HD product. Binned-mode barcodes encode their own exact grid position (see [Validation](#validation)); cell-segmentation-mode barcodes don't use that convention. Since 0.2.2, position is computed as the true area centroid of each cell's embedded polygon boundary (`CellSegs.GeoJSON`) — confirmed coordinate-for-coordinate identical to SpaceRanger's own segmented output — which matches official ground truth exactly. Earlier versions averaged a coarser rectangular approximation of the cell mask instead (kept today only as a fallback for the rare barcode without polygon data): 0.1.1 fixed a real bug in that averaging (`Centers` read with the wrong array layout, scrambling ~30-40% of cells onto a diagonal), and 0.2.2 replaced the averaging approach itself, which still had a small residual (~2 µm) even once correctly read. `array_row`/`array_col` are unconditionally `(0, 0)` for cell-segmentation-mode `.cloupe` files — there's no fixed grid for individual cells to round to — which doesn't affect plotting (that uses the real pixel coordinates) but is worth knowing if you read those columns directly.

**Prefer official SpaceRanger output when you have it — but note `squidpy.read.visium()` itself can't read modern Visium HD directories.** Confirmed directly against a real sample: `squidpy.read.visium()` (as of squidpy 1.8.2) only looks for `spatial/tissue_positions.csv` (or the older `tissue_positions_list.csv`), but current SpaceRanger HD output writes `spatial/tissue_positions.parquet` instead — `squidpy.read.visium()` raises `FileNotFoundError` outright on a real, current HD `outs/` directory. This is consistent with squidpy's own maintainers pointing HD users at `spatialdata` instead (see above), and is why this package's own validation (below) reads the official `.h5`/`.parquet` files directly rather than through `squidpy.read.visium()`.

### What a `.cloupe` file can't give you

Even for common cases, `.cloupe` doesn't carry everything in a full SpaceRanger output tree. Confirmed directly against real files:

| Not in `.cloupe` | What it's for |
|---|---|
| Raw (unfiltered) count matrix | Ambient-RNA correction, background-signal QC, custom re-filtering — `.cloupe`'s barcode set matches SpaceRanger's *filtered* matrix only |
| Molecule-level data (`molecule_info.h5`) | Per-UMI/per-read records for saturation analysis |
| Aligned reads (`possorted_genome_bam.bam`) | Variant calling, allele-specific expression, splice-junction/isoform analysis |
| Probe-level data (`probe_set.csv`, `raw_probe_bc_matrix.h5`) | Anything where gene-level aggregation is inappropriate (probe-based/FFPE chemistry) |
| `feature_slice.h5` | HD-specific rasterized per-gene expression image data |
| `barcode_mappings.parquet` | The exact HD bin-to-bin crosswalk (which 2 µm bins compose which 8 µm bin) — `combine_cloupe_bins()` (Loupe2R) merges resolutions as independent assays with no containment relationship between their bins, a practical but less precise substitute |
| Additional image products (`cytassist_image.tiff`, `aligned_fiducials.jpg`, `detected_tissue_image.jpg`) | Redoing tissue detection or checking fiducial-alignment quality |
| Full QC/metrics reporting (`web_summary.html`, `metrics_summary.csv`) | The `.cloupe` sections that would carry this (`Metrics`, `Analyses`) exist in the format but are confirmed empty in every real file checked — not partially present, genuinely absent |
| Differential expression / PCA tables | Confirmed absent, not just unextracted — the `Analyses` section contains only a null placeholder entry in every file checked |

## Credits

The hard part — reverse-engineering the proprietary `.cloupe` binary format at all (the header layout, byte-offset index block, matrix/projection encoding, and tiled image storage) — is not this package's work. It's [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe), written by **Martin Prete** and **Nithin Mathew Joseph** of the Wellcome Sanger Institute's Cellular Genetics Informatics (cellgeni) team, AGPL-3.0 licensed. A pinned copy is vendored directly into [`cloupe_extract`](cloupe_extract/), the extraction core `Loupe2Py` depends on, at `cloupe_extract/src/cloupe_extract/_vendor/cloupe.py` (see that file's header for the exact commit) rather than kept as an external dependency, to eliminate a class of path-configuration bugs and pin an exact, tested parser version. Full license text: [`THIRD_PARTY_LICENSES/cellgeni-cloupe-AGPL-3.0.txt`](THIRD_PARTY_LICENSES/cellgeni-cloupe-AGPL-3.0.txt).

**Note:** vendoring this AGPL-3.0 code makes `cloupe_extract` a combined work under AGPL-3.0 terms, and `Loupe2Py` depends on it inseparably in practice (there's no meaningful way to use `Loupe2Py` without it), so `Loupe2Py` is itself licensed **AGPL-3.0-or-later** too — see [License](#license) below.

Also indebted to the [scverse](https://scverse.org/) ecosystem ([AnnData](https://anndata.readthedocs.io/), [scanpy](https://scanpy.readthedocs.io/), [squidpy](https://squidpy.readthedocs.io/)) for the conventions and reference outputs this package's Visium HD support was validated against, and to [10x Genomics' SpaceRanger](https://www.10xgenomics.com/support/software/space-ranger) for the official output used as ground truth (see [Validation](#validation)).

## Validation

`cloupe_to_anndata()` was cross-validated against the official SpaceRanger output for the same Visium HD sample used to validate Loupe2R (8µm bins, 19,072 features × 472,752 barcodes). Results:

| Check | Result |
|---|---|
| Barcode overlap (cloupe vs. official filtered matrix) | 100% |
| Count exact match vs. official (not just correlated) | 100% |
| `array_row`/`array_col` exact match vs. official `tissue_positions.parquet` | 100% (0 max difference) |
| `pxl_col`/`pxl_row` correlation | 1.000000 |

Identical results to Loupe2R's own validation on the same sample, as expected — same extraction logic underneath. Also smoke-tested against real squidpy functions (`squidpy.gr.spatial_neighbors()`, `squidpy.pl.spatial_scatter()`), not just hand-written assertions.

The `array_row`/`array_col` result relies on the same fix validated in Loupe2R: Visium HD barcodes encode their exact grid position directly (e.g. `s_008um_00269_00526-1` → row 269, column 526), and `cloupe_extract.extract.parse_array_position()` parses that directly rather than approximating it from pixel coordinates (which was found to disagree with the official grid by hundreds of bins due to a pixel-origin mismatch between the stitched tissue image and SpaceRanger's own coordinate frame). Non-HD barcodes fall back to the pixel-based approximation.

Reproduce this validation yourself (see `tests/test_integration_visium_hd.py`) by setting `LOUPE2PY_TEST_DIR` to a SpaceRanger Visium HD `outs/` directory containing `cloupe_008um.cloupe` and `binned_outputs/square_008um/`, and running `pytest` with the `test` extras installed.

## License

AGPL-3.0-or-later (see `LICENSE`). This package depends on `cloupe_extract`, which vendors AGPL-3.0 code ([Credits](#credits) above) directly into its own source — under AGPL-3.0's terms that makes `cloupe_extract` a combined work, and `Loupe2Py` inherits the same license as a practical matter, since it's not usable without `cloupe_extract`. `THIRD_PARTY_LICENSES/` contains the full license text for the vendored code specifically.
