# cloupe_extract 0.2.0

## Breaking changes

- **An unrecognized `.cloupe` format version now aborts `extract_cloupe()` by default, instead of just warning.** Previously, a file reporting a format version outside `_TESTED_VERSIONS` printed a warning and extraction proceeded anyway — silently guessing on an undocumented, reverse-engineered format is the wrong default for a data-recovery tool. `extract_cloupe()` now raises `UnvalidatedFormatVersionError` before extracting anything, unless called with the new `version_check=False`, which restores the previous warn-and-proceed behavior and puts the responsibility for verifying the result on the caller. `check_format_version()` itself is unchanged — it only ever inspects and returns; `extract_cloupe()` is what acts on the result.
- `format_info.json` (detected versions + warnings) is still written even when `version_check=True` raises, so the specific detected versions remain on record — though note it won't survive past the exception if you're using the default auto-cleaned temp directory (pass `outdir=`/`keep_files=True` in `loupe2py`/`Loupe2R`, or your own `outdir` here, if you want it to persist).
- The CLI (`python -m cloupe_extract.extract`) gained a `--no-version-check` flag, the equivalent of `version_check=False`.
- New: `UnvalidatedFormatVersionError`, exported from `cloupe_extract` (`from cloupe_extract import UnvalidatedFormatVersionError`) and re-exported from `loupe2py` for convenience.

# cloupe_extract 0.1.0

Initial release: split out of `loupe2py` into its own independently pip-installable package, so `Loupe2R` (the Seurat sibling package) could depend on something more aptly named than a package literally called `loupe2py`. No behavior change from `loupe2py` 0.2.2 — same code, moved. See `loupe2py`'s and `Loupe2R`'s own `NEWS.md` for the split itself.
