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

Imports one `.cloupe` file — Visium HD, in either binned or cell-segmentation mode — into a squidpy-ready `AnnData` object, following the same `obsm['spatial']` / `uns['spatial'][library_id]` convention `scanpy`/`squidpy` use. Calls `extract_cloupe()` internally, then builds the count matrix (as `layers["counts"]` plus the working `.X`), `obs` (library identity, QC metrics, `array_row`/`array_col`, Space Ranger clusterings as `sr_`-prefixed categoricals, user cell tracks), `var`, `obsm['spatial']`, `obsm['X_<name>']` for other embeddings, and `uns['spatial'][library_id]` (images, scale factors). Also stashes `bin_size_um` and `cloupe_format_info` in `uns`.

- `cloupe_path` — path to the `.cloupe` file.
- `library_id` — identifier for this sample; becomes the `uns['spatial']` key and `obs['library_id']`. Defaults to the `.cloupe` filename with its extension stripped.
- `include_image` — whether to reconstruct and embed the tissue image.
- `outdir`, `keep_files` — control where intermediate extracted files go; `None`/`False` uses an auto-deleted temp directory.

Returns an `anndata.AnnData` object.

### `extract_cloupe()`

```python
extract_cloupe(cloupe_path, outdir, include_image=True)
```

Lower-level: parses a `.cloupe` file and writes the extracted data as generic, SpaceRanger-convention-like files into `outdir` (`matrix.mtx.gz`, `barcodes.tsv.gz`, `features.tsv.gz`, `tissue_positions.csv`, `scalefactors_json.json`, `tissue_hires_image.png`, `projections.csv`, `clusterings.csv`, `celltracks.csv`, `format_info.json`) — no `AnnData` object is built. `cloupe_to_anndata()` calls this internally; call it directly if you want the raw files instead (e.g., to point another tool at them, or to build a different object type entirely). Implemented in [`cloupe_extract`](cloupe_extract/), the standalone extraction package this package depends on, and re-exported here (`from loupe2py import extract_cloupe`) for convenience — `from cloupe_extract import extract_cloupe` gets you the identical function. This is also what `Loupe2R::cloupe_to_seurat()` calls directly via `reticulate`, without going through `loupe2py` at all.

Returns `outdir`.

### CLI

```bash
python -m cloupe_extract.extract <cloupe_path> <outdir>
```

Command-line wrapper around `extract_cloupe()`, for use outside Python (a shell script, a pipeline step). Provided by `cloupe_extract`, not `loupe2py` itself.

## Internal

Not in `__all__`, no compatibility guarantee.

### `loupe2py.anndata_io` — underscore-prefixed, private by Python convention

Direct ports of the equivalent `Loupe2R` R helpers (`R/utils.R`) — same logic, same edge cases.

| Function | What it does |
|---|---|
| `_detect_mt_pattern(var_names)` | Returns `"MT-"` if any gene name starts with it (human convention), else `"mt-"` (mouse convention) — even if nothing matches, in which case `pct_counts_mt` ends up 0. |
| `_derive_sample_name(cloupe_path)` | Strips the directory and final extension from a `.cloupe` path to derive the default `library_id`. |
| `_align_positions(df, cell_names)` | Reindexes a `DataFrame` with a `barcode` column to a given cell-name order. |

### `cloupe_extract.extract` — implementation detail of the dependency, not this package

The functions `extract_cloupe()` actually calls (`check_format_version()`, `stitch_tiles()`, `get_spatial_projection()`, `get_cellseg_projection()`, `parse_array_position()`, `read_clusterings()`, ...) live in `cloupe_extract`, not `loupe2py` — see [`cloupe_extract`'s own source](cloupe_extract/src/cloupe_extract/extract.py) for the full list and behavior. Listed here only as a pointer, since `loupe2py`'s own `anndata_io.py` calls into them indirectly through `extract_cloupe()`.
