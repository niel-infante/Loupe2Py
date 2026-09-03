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
