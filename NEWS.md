# loupe2py 0.2.2

## Bug fix

- **Fixed imprecise cell-segmentation positions by reading the exact boundary data the file already contained.** `get_cellseg_projection()` computed each cell's position by averaging the centers of `CellSegs`' `Centers`/`Sizes` rects — a coarse rectangular decomposition of the cell mask. On a real paired sample (Visium HD Human Kidney, SpaceRanger 4.0.1), that rect decomposition covered only ~69% of the 2 µm bins SpaceRanger itself assigns to each cell, biasing the averaged position by a small but real amount (mean ~7 px / ~2 µm off official ground truth).
- `CellSegs` also carries a `GeoJSON` field — the exact polygon boundary for every cell, confirmed coordinate-for-coordinate identical to SpaceRanger's own `segmented_outputs/cell_segmentations.geojson` — that was never read. `get_cellseg_projection()` now computes each cell's true area centroid from this polygon data via the shoelace formula, and falls back to the rect-average only for the rare barcode (or whole file) without `GeoJSON`. Re-validated against the same real kidney sample: position now matches official ground truth **exactly** (0.0 px mean/median/max distance across all 148,056 cells), not just closely.
- Confirmed all 6 real cell-segmented `.cloupe` files checked across this project (5 private, 1 public) have `GeoJSON` present, so this is expected to be the common path, not an edge case.
- If you've used `cloupe_to_seurat()`/`cloupe_to_anndata()` on cell-segmentation-mode `.cloupe` files, re-extract for exact positions; the previous result was close (small px-scale bias) but not exact.

# loupe2py 0.2.1

## Bug fix

- **Fixed inflated count totals from synthetic per-file rows leaking into the count matrix.** Every `.cloupe` file's `Matrices` section appends built-in aggregate rows after the real genes — one `type_sum_<FeatureType>` row per feature type and one `genome_sum_<Reference>` row per reference genome — each holding that barcode's *entire* total for the category, not a real gene's count. `extract_cloupe()` previously wrote these out as if they were ordinary genes, in `features.tsv.gz`/`matrix.mtx.gz` and therefore in every `AnnData`/`Seurat` object built from them. For a typical single-feature-type, single-genome sample (the common case), this meant every cell's real gene sum, a `type_sum` row equal to that same sum, and a `genome_sum` row equal to it again — totaling exactly **3x** the correct value for standard per-cell metrics (`nCount_Spatial`/`total_counts`), and silently deflating `percent.mt` by the same factor.
- New `exclude_synthetic_totals()` drops these rows (matched by ID prefix, not position, so it generalizes to files with multiple feature types or reference genomes) before the count matrix, feature list, or any downstream metric is built.
- **If you've run `cloupe_to_seurat()` or `cloupe_to_anndata()` on any `.cloupe` file with a `loupe2py` version before this one, re-extract** — per-cell total counts and `percent.mt` were wrong, and two bogus "genes" (`type_sum_...`, `genome_sum_...`) were present in the feature list. Confirmed via paired official SpaceRanger output on two independent real Visium HD samples: after this fix, every real (gene, barcode) count matches exactly (not just correlates) — this was purely an extra-rows issue, not corruption of the real per-gene values themselves.
- This affects any `.cloupe` file, not just probe-based/FFPE chemistry specifically — it's a property of the `.cloupe` format's own `Matrices` section, not something specific to the samples it was first found on.
- Also extends `_TESTED_VERSIONS` to cover container format `8.0.0` and matrix format `6.2.0`, validated on the same sample that surfaced this bug.

# loupe2py 0.2.0

## License change

- `loupe2py` is now formally licensed **AGPL-3.0-or-later**, resolving the "not yet finalized" note that had been carried since initial development. This package vendors AGPL-3.0-licensed code from [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe) directly into its own source (`src/loupe2py/_vendor/cloupe.py`), which makes the combined work AGPL-3.0 under that license's terms. See the README's [License](https://github.com/niel-infante/Loupe2Py#license) section.

# loupe2py 0.1.1

## Bug fix

- **Fixed spatial coordinate corruption for cell-segmentation `.cloupe` files.** `get_cellseg_projection()` read the `CellSegs` `Centers` array as interleaved `(x, y)` pairs; it's actually stored as struct-of-arrays (all x values, then all y values). This silently scrambled pixel coordinates for 30-40% of cells in every real cell-segmented `.cloupe` file checked, producing a visible diagonal artifact (`pxl_row_in_fullres ≈ pxl_col_in_fullres`) in spatial plots. Also removed a second, previously-unreported bug in the same function: `Centers` values are already in full-resolution image pixel coordinates and were being incorrectly divided by `MicronsPerPixel` on top, mis-scaling every cell-segmented sample (not just the diagonal-affected fraction).
- If you've run `cloupe_to_seurat()` (via Loupe2R) or `cloupe_to_anndata()` on any `_cell.cloupe` file with a `loupe2py` version before this one, re-extract — the spatial coordinates were wrong.
- This does not affect Visium HD (binned) `.cloupe` files, which use a different code path.
- Clarified (no behavior change): `array_row`/`array_col` are unconditionally `(0, 0)` for cell-segmentation data, since there's no fixed grid for cell segments to round to — not an approximation, just not a meaningful value for this data type. Real pixel coordinates (what drives plot placement) are unaffected.

# loupe2py 0.1.0

Initial release: `cloupe_to_anndata()`.
