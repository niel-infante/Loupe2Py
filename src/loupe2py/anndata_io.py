"""
Assemble a squidpy-ready AnnData object from an extracted .cloupe file.

squidpy has no object type of its own -- it operates directly on AnnData
following the convention scanpy's read_visium() established:
    adata.obsm['spatial']                         spot pixel coordinates
    adata.uns['spatial'][library_id]['images']     'hires'/'lowres' arrays
    adata.uns['spatial'][library_id]['scalefactors']
    adata.uns['spatial'][library_id]['metadata']

This module's job is exactly what R/cloupe_to_seurat.R does for Seurat:
read back the generic files cloupe_extract.extract_cloupe() already wrote,
and assemble them into that convention.
"""

import json
import os
import re
import shutil
import tempfile
import warnings

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
from PIL import Image

from cloupe_extract import extract as _extract


def _detect_mt_pattern(var_names):
    """Return 'MT-' if any human-style mitochondrial genes are present,
    otherwise 'mt-' (mouse-style) -- matching the fallback assumed even when
    neither pattern actually matches (pct_counts_mt is then simply 0)."""
    return "MT-" if any(str(v).startswith("MT-") for v in var_names) else "mt-"


def _derive_sample_name(cloupe_path):
    """Basename with the final extension stripped (keeps internal dots,
    e.g. 'sample.v2.cloupe' -> 'sample.v2')."""
    base = os.path.basename(cloupe_path)
    return re.sub(r"\.[^.]+$", "", base)


def _align_positions(df, cell_names):
    """Reindex df (must have a 'barcode' column) to cell_names, matching the
    R side's .align_positions() -- rows for names absent from df become NaN."""
    return df.set_index("barcode").reindex(cell_names)


def cloupe_to_anndata(
    cloupe_path,
    library_id=None,
    include_image=True,
    outdir=None,
    keep_files=False,
):
    """Import a 10x Genomics .cloupe file into a squidpy-ready AnnData object.

    Calls cloupe_extract.extract.extract_cloupe() to parse the .cloupe binary,
    then assembles an AnnData object with the count matrix, spatial
    coordinates, tissue image, UMAP embedding, and Space Ranger cluster
    labels -- following the same obsm['spatial']/uns['spatial'] convention
    scanpy's read_visium() (and squidpy.read.visium()) use.

    Parameters
    ----------
    cloupe_path : str
        Path to the .cloupe file.
    library_id : str, optional
        Identifier for this sample. Becomes the uns['spatial'] dict key and
        obs['library_id']. Defaults to the .cloupe filename (extension
        stripped).
    include_image : bool
        Whether to extract and embed the tissue image. Requires Pillow.
    outdir : str, optional
        Directory for intermediate extracted files. None (default) uses an
        auto-created temp dir, deleted on return unless keep_files=True.
    keep_files : bool
        If True and outdir is None, keep the intermediate files (useful for
        debugging or cross-validation).

    Returns
    -------
    anndata.AnnData
        With obs columns library_id, in_tissue, array_row, array_col,
        total_counts, n_genes_by_counts, pct_counts_mt, plus any Space
        Ranger clusterings (prefixed sr_) and user cell tracks; obsm
        'spatial' and any embeddings (X_umap, etc.); uns['spatial'],
        uns['cloupe_format_info'], uns['bin_size_um'].
    """
    cloupe_path = os.path.abspath(os.path.expanduser(cloupe_path))
    if not os.path.exists(cloupe_path):
        raise FileNotFoundError(cloupe_path)

    cleanup = outdir is None and not keep_files
    if outdir is None:
        outdir = tempfile.mkdtemp(prefix="cloupe_")

    try:
        print(f"Extracting from {os.path.basename(cloupe_path)} "
              "(may take several minutes)...")
        _extract.extract_cloupe(cloupe_path, outdir, include_image=include_image)

        print("Building AnnData object...")

        # ---- Count matrix: mtx (features x barcodes) -> obs x var ----
        mat = scipy.io.mmread(os.path.join(outdir, "matrix.mtx.gz"))
        X = mat.T.tocsr().astype(np.int32)

        barcodes = pd.read_csv(
            os.path.join(outdir, "barcodes.tsv.gz"), header=None
        )[0].astype(str).values
        feats = pd.read_csv(
            os.path.join(outdir, "features.tsv.gz"), header=None, sep="\t",
            names=["feature_id", "feature_name", "feature_type"],
        )

        var = pd.DataFrame(
            {
                "gene_ids": feats["feature_id"].values,
                "feature_types": feats["feature_type"].values,
            },
            index=pd.Index(feats["feature_name"].astype(str), name=None),
        )
        obs = pd.DataFrame(index=pd.Index(barcodes, name=None))

        adata = ad.AnnData(X=X, obs=obs, var=var)
        adata.var_names_make_unique()
        adata.layers["counts"] = adata.X.copy()

        # ---- Sample identity ----
        lib = library_id or _derive_sample_name(cloupe_path)
        adata.obs["library_id"] = lib

        # ---- percent-mito / basic QC (no scanpy dependency) ----
        mt_pattern = _detect_mt_pattern(adata.var_names)
        mt_mask = adata.var_names.str.startswith(mt_pattern)
        total_counts = np.asarray(adata.X.sum(axis=1)).ravel()
        mt_counts = np.asarray(adata.X[:, mt_mask].sum(axis=1)).ravel()
        adata.obs["total_counts"] = total_counts
        adata.obs["n_genes_by_counts"] = adata.X.getnnz(axis=1)
        adata.obs["pct_counts_mt"] = np.divide(
            mt_counts, total_counts,
            out=np.zeros_like(mt_counts, dtype=float),
            where=total_counts > 0,
        ) * 100
        # Column names deliberately match sc.pp.calculate_qc_metrics()'s own
        # output names, so results stay consistent if the user re-runs it.

        # ---- Spatial coordinates ----
        pos = pd.read_csv(os.path.join(outdir, "tissue_positions.csv"))
        pos = _align_positions(pos, adata.obs_names)
        adata.obs["in_tissue"] = pos["in_tissue"].astype(int).values
        adata.obs["array_row"] = pos["array_row"].values
        adata.obs["array_col"] = pos["array_col"].values
        # Select by column NAME (tissue_positions.csv's header is trustworthy,
        # unlike historical headerless SpaceRanger CSVs) rather than
        # replicating squidpy's position-based relabeling workaround for
        # that -- same final (x, y) convention, more transparent.
        adata.obsm["spatial"] = pos[
            ["pxl_col_in_fullres", "pxl_row_in_fullres"]
        ].to_numpy()

        # ---- uns['spatial'][library_id]: images / scalefactors / metadata ----
        with open(os.path.join(outdir, "scalefactors_json.json")) as f:
            sf = json.load(f)

        images = {}
        img_path = os.path.join(outdir, "tissue_hires_image.png")
        if include_image and os.path.exists(img_path):
            print("Loading tissue image...")
            # float32 in [0,1], NOT Pillow's default uint8 [0,255] --
            # matplotlib.image.imread() (what scanpy/squidpy actually use to
            # load these PNGs) returns float32 [0,1]; storing uint8 here
            # would render spatial plots ~255x too bright.
            arr = np.asarray(Image.open(img_path).convert("RGB"), dtype=np.float32) / 255.0
            # Both scalefactors are 1.0 by design (see extract.py) since only
            # one raster is ever stored; register it under both keys so
            # nothing defaulting to img_key="lowres" hits a KeyError.
            images["hires"] = arr
            images["lowres"] = arr
        elif include_image:
            warnings.warn("No tissue image found. Install Pillow: pip install Pillow")

        adata.uns["spatial"] = {lib: {
            "images": images,
            "scalefactors": sf,
            "metadata": {},
        }}

        # ---- format_info / bin_size_um -> uns (Misc() analog) ----
        with open(os.path.join(outdir, "format_info.json")) as f:
            fmt_info = json.load(f)
        if fmt_info["warnings"]:
            warnings.warn(
                "Unvalidated .cloupe format version(s) detected:\n  "
                + "\n  ".join(fmt_info["warnings"])
                + "\nExtraction proceeded, but treat results with extra scrutiny."
            )
        adata.uns["cloupe_format_info"] = fmt_info["versions"]
        if sf.get("bin_size_um") is not None:
            adata.uns["bin_size_um"] = sf["bin_size_um"]

        # ---- Non-spatial projections (UMAP, tSNE, ...) -> obsm['X_<name>'] ----
        proj_path = os.path.join(outdir, "projections.csv")
        if os.path.exists(proj_path):
            proj = pd.read_csv(proj_path)
            proj = _align_positions(proj, adata.obs_names)
            names = sorted({re.sub(r"_\d+$", "", c) for c in proj.columns})
            for name in names:
                cols = [c for c in proj.columns
                        if re.fullmatch(rf"{re.escape(name)}_\d+", c)]
                adata.obsm[f"X_{name.lower()}"] = proj[cols].to_numpy()
            print(f"Added reductions: {', '.join(n.lower() for n in names)}")

        # ---- Space Ranger clusterings -> obs['sr_<name>'] (categorical) ----
        cl_path = os.path.join(outdir, "clusterings.csv")
        if os.path.exists(cl_path):
            cl = pd.read_csv(cl_path)
            cl = _align_positions(cl, adata.obs_names)
            for col in cl.columns:
                adata.obs[f"sr_{col}"] = pd.Categorical(cl[col])
            print(f"Added {len(cl.columns)} Space Ranger clustering(s) (prefix 'sr_')")

        # ---- User cell tracks -> obs ----
        ct_path = os.path.join(outdir, "celltracks.csv")
        if os.path.exists(ct_path):
            ct = pd.read_csv(ct_path)
            ct = _align_positions(ct, adata.obs_names)
            for col in ct.columns:
                adata.obs[col] = ct[col].values
            print(f"Added {len(ct.columns)} user cell track(s)")

        print(
            f"Done!  {adata.n_obs} spots x {adata.n_vars} genes  |  "
            f"library_id: {lib}  |  "
            f"image: {'yes' if images else 'none'}"
        )
        return adata
    finally:
        if cleanup:
            shutil.rmtree(outdir, ignore_errors=True)
