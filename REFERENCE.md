# loupe2py function reference

Every function in the package, public and internal. If you're calling something from outside the package, use only the **Public API** section — everything under **Internal** can change or disappear without notice and isn't part of any compatibility guarantee, even the entries that aren't underscore-prefixed (Python doesn't enforce privacy the way R's export list does — "not in `__all__`" is the actual signal here, not the name).

## Public API

In `loupe2py.__all__`, importable as `from loupe2py import ...`, stable.

### `cloupe_to_anndata()`

```python
cloupe_to_anndata(
    cloupe_path,
    library_id=None,
    include_image=True,
    outdir=None,
    keep_files=False,
)
```

Imports one `.cloupe` file — Visium HD or standard Visium with cell segmentation — into a squidpy-ready `AnnData` object, following the same `obsm['spatial']` / `uns['spatial'][library_id]` convention `scanpy`/`squidpy` use. Calls `extract_cloupe()` internally, then builds the count matrix (as `layers["counts"]` plus the working `.X`), `obs` (library identity, QC metrics, `array_row`/`array_col`, Space Ranger clusterings as `sr_`-prefixed categoricals, user cell tracks), `var`, `obsm['spatial']`, `obsm['X_<name>']` for other embeddings, and `uns['spatial'][library_id]` (images, scale factors). Also stashes `bin_size_um` and `cloupe_format_info` in `uns`.

- `cloupe_path` — path to the `.cloupe` file.
- `library_id` — identifier for this sample; becomes the `uns['spatial']` key and `obs['library_id']`. Defaults to the `.cloupe` filename with its extension stripped.
- `include_image` — whether to reconstruct and embed the tissue image.
- `outdir`, `keep_files` — control where intermediate extracted files go; `None`/`False` uses an auto-deleted temp directory.

Returns an `anndata.AnnData` object.

### `extract_cloupe()`

```python
extract_cloupe(cloupe_path, outdir, include_image=True)
```

Lower-level: parses a `.cloupe` file and writes the extracted data as generic, SpaceRanger-convention-like files into `outdir` (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`, `tissue_positions.csv`, `scalefactors_json.json`, `tissue_hires_image.png`, `projections.csv`, `clusterings.csv`, `celltracks.csv`, `format_info.json`) — no `AnnData` object is built. `cloupe_to_anndata()` calls this internally; call it directly if you want the raw files instead (e.g., to point another tool at them, or to build a different object type entirely). This is also what `Loupe2R::cloupe_to_seurat()` calls via `reticulate`.

Returns `outdir`.

### CLI

```bash
python -m loupe2py.extract <cloupe_path> <outdir>
```

Command-line wrapper around `extract_cloupe()`, for use outside Python (a shell script, a pipeline step).

## Internal

Not in `__all__`, no compatibility guarantee.

### `loupe2py.extract` — importable but not part of the public interface

| Function | What it does |
|---|---|
| `check_format_version(cloupe_obj)` | Compares the file's internal format-version fields (container, run, matrix, projection) against a known-tested table; returns detected versions and any warnings. Called automatically by `extract_cloupe()`. |
| `stitch_tiles(cloupe_obj, target_level=None)` | Reconstructs the full-resolution tissue image from the file's internal tile pyramid; auto-selects the highest zoom level and normalizes to RGB. |
| `get_spatial_projection(cloupe_obj)` | Returns the Visium HD `"Spatial"` projection's pixel coordinates and bin size; falls back to `get_cellseg_projection()` when the file has no HD projection. |
| `get_cellseg_projection(cloupe_obj)` | Returns per-barcode pixel coordinates for cell-segmented (non-HD) `.cloupe` files, by averaging per-cell-segment centroids from the `CellSegs` data block. |
| `parse_array_position(barcode, pxl_row, pxl_col, bin_px)` | Returns `(array_row, array_col)` for one barcode: parsed exactly from a Visium HD barcode string when it matches that format, else rounded from pixel coordinates, else `(0, 0)`. |
| `read_clusterings(cloupe_obj)` | Returns Space Ranger's graph/k-means cluster label assignments per barcode. |

### `loupe2py.anndata_io` — underscore-prefixed, private by Python convention

Direct ports of the equivalent `Loupe2R` R helpers (`R/utils.R`) — same logic, same edge cases.

| Function | What it does |
|---|---|
| `_detect_mt_pattern(var_names)` | Returns `"MT-"` if any gene name starts with it (human convention), else `"mt-"` (mouse convention) — even if nothing matches, in which case `pct_counts_mt` ends up 0. |
| `_derive_sample_name(cloupe_path)` | Strips the directory and final extension from a `.cloupe` path to derive the default `library_id`. |
| `_align_positions(df, cell_names)` | Reindexes a `DataFrame` with a `barcode` column to a given cell-name order. |
