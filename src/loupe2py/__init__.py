"""Import 10x Genomics .cloupe files (including Visium HD) into squidpy-ready AnnData objects."""

from .anndata_io import cloupe_to_anndata
from .extract import extract_cloupe

__all__ = ["cloupe_to_anndata", "extract_cloupe"]

__version__ = "0.1.1"
