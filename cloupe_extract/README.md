# cloupe_extract

Shared `.cloupe` extraction core underlying both [`loupe2py`](https://github.com/niel-infante/Loupe2Py) (Python, builds `AnnData` objects) and [`Loupe2R`](https://github.com/niel-infante/Loupe2R) (R, builds `Seurat` objects via `reticulate`). Not meant to be installed directly by end users — install `loupe2py` or `Loupe2R` instead, which depend on this package automatically.

Parses a `.cloupe` file — 10x Genomics' proprietary, undocumented format for Loupe Browser — and writes generic, SpaceRanger-convention-like intermediate files (count matrix, barcodes, features, spatial positions, tissue image, clusterings, cell tracks). Framework-agnostic: no knowledge of `AnnData` or `Seurat` itself.

Built on a vendored, pinned copy of [`cellgeni/cloupe`](https://github.com/cellgeni/cloupe) (Wellcome Sanger Institute, Cellular Genetics Informatics) in `src/cloupe_extract/_vendor/cloupe.py` — see [Loupe2Py's README](https://github.com/niel-infante/Loupe2Py#credits) for full credits, licensing rationale, and validation details, which cover this package too.

By default, `extract_cloupe()` raises `UnvalidatedFormatVersionError` and extracts nothing if a file reports an internal `.cloupe` format version outside the known-tested set — the parser may well still handle it correctly, but that hasn't been checked, and this package would rather refuse than guess silently. Pass `version_check=False` to proceed anyway (a warning is printed instead); you're then responsible for independently verifying the result.

## License

AGPL-3.0-or-later (see `LICENSE`) — this package vendors AGPL-3.0 code directly (`cellgeni/cloupe`), which makes the combined work AGPL-3.0 under that license's terms. `THIRD_PARTY_LICENSES/` contains the full license text for the vendored code specifically.
