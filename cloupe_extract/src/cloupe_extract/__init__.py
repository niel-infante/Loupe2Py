"""Shared .cloupe extraction core underlying both loupe2py (Python, AnnData)
and Loupe2R (R, Seurat, via reticulate). Framework-agnostic: parses a
.cloupe file and writes generic, SpaceRanger-convention-like intermediate
files; has no knowledge of AnnData or Seurat itself.
"""

from .extract import extract_cloupe, UnvalidatedFormatVersionError

__all__ = ["extract_cloupe", "UnvalidatedFormatVersionError"]

__version__ = "0.2.0"
