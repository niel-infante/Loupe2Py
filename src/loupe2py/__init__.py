"""Import 10x Genomics .cloupe files (including Visium HD) into squidpy-ready AnnData objects."""

from cloupe_extract import extract_cloupe

from .anndata_io import cloupe_to_anndata

__all__ = ["cloupe_to_anndata", "extract_cloupe"]

__version__ = "0.3.0"
