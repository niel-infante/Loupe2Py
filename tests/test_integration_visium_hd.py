# Opt-in: cross-validates cloupe_to_anndata() against the official
# SpaceRanger output for the same sample. Not run by default (no real
# .cloupe file is committed to this repo) -- set LOUPE2PY_TEST_DIR to a
# SpaceRanger Visium HD `outs/` directory containing cloupe_008um.cloupe
# and binned_outputs/square_008um/ to run it, e.g.:
#   LOUPE2PY_TEST_DIR=/path/to/outs pytest tests/test_integration_visium_hd.py

import hashlib
import json
import os

import numpy as np
import pandas as pd
import pytest

test_dir = os.environ.get("LOUPE2PY_TEST_DIR", "")
pytestmark = pytest.mark.skipif(
    not test_dir, reason="Set LOUPE2PY_TEST_DIR to run the local Visium HD integration test"
)

FIXTURE_PATH = os.path.join(os.path.dirname(__file__), "fixtures", "cloupe_008um_expected.json")


@pytest.fixture(scope="module")
def cloupe_file():
    path = os.path.join(test_dir, "cloupe_008um.cloupe")
    if not os.path.exists(path):
        pytest.skip(f"cloupe_008um.cloupe not found in {test_dir}")
    return path


@pytest.fixture(scope="module")
def official_dir():
    path = os.path.join(test_dir, "binned_outputs", "square_008um")
    if not os.path.isdir(path):
        pytest.skip(f"binned_outputs/square_008um not found in {test_dir}")
    return path


@pytest.fixture(scope="module")
def extracted(tmp_path_factory, cloupe_file):
    from loupe2py import cloupe_to_anndata

    outdir = tmp_path_factory.mktemp("cloupe_extract")
    adata = cloupe_to_anndata(cloupe_file, outdir=str(outdir), keep_files=True)
    return adata, str(outdir)


def test_extraction_matches_regression_fixture(extracted):
    adata, outdir = extracted
    with open(FIXTURE_PATH) as f:
        fixture = json.load(f)

    assert adata.n_vars == fixture["n_features"]
    assert adata.n_obs == fixture["n_barcodes"]
    assert int(adata.layers["counts"].sum()) == fixture["total_umi"]
    assert adata.uns["bin_size_um"] == fixture["bin_size_um"]

    lib = adata.obs["library_id"].iloc[0]
    assert adata.uns["spatial"][lib]["scalefactors"]["spot_diameter_fullres"] == pytest.approx(
        fixture["spot_diameter_fullres"], abs=1e-5
    )

    bc_digest = hashlib.sha256(
        "\n".join(sorted(adata.obs_names)).encode("utf-8")
    ).hexdigest()
    assert bc_digest == fixture["barcode_set_sha256"]


def test_no_unvalidated_format_version_for_this_file(extracted):
    adata, _ = extracted
    with open(FIXTURE_PATH) as f:
        fixture = json.load(f)

    fmt = adata.uns["cloupe_format_info"]
    assert fmt["container"] == fixture["cloupe_format_info"]["container"]
    assert fmt["matrix"] == fixture["cloupe_format_info"]["matrix"]


def test_agrees_with_official_spaceranger_output(extracted, official_dir):
    import scanpy as sc

    adata, outdir = extracted

    # Deliberately NOT using squidpy.read.visium() here: it only looks for
    # spatial/tissue_positions.csv (or the older tissue_positions_list.csv),
    # but modern SpaceRanger Visium HD output writes tissue_positions.parquet
    # instead -- confirmed on this real sample, squidpy.read.visium() raises
    # FileNotFoundError on it outright. That's a genuine squidpy limitation
    # (consistent with squidpy's own maintainers pointing HD users at
    # spatialdata instead), not something to route around silently -- read
    # the official filtered/raw .h5 matrices and tissue_positions.parquet
    # directly instead, which is still genuine official SpaceRanger ground
    # truth, just without going through squidpy's own (HD-incompatible)
    # directory reader.
    official = sc.read_10x_h5(os.path.join(official_dir, "filtered_feature_bc_matrix.h5"))
    raw = sc.read_10x_h5(os.path.join(official_dir, "raw_feature_bc_matrix.h5"))

    cloupe_bc = set(adata.obs_names)
    overlap_filtered = len(cloupe_bc & set(official.obs_names)) / len(cloupe_bc)
    overlap_raw = len(cloupe_bc & set(raw.obs_names)) / len(cloupe_bc)
    print(f"Barcode overlap vs filtered official: {overlap_filtered:.4f}, "
          f"vs raw official: {overlap_raw:.4f}")
    assert max(overlap_filtered, overlap_raw) > 0.95

    shared = sorted(cloupe_bc & set(official.obs_names))
    assert len(shared) > 0

    nc_cloupe = np.asarray(adata[shared].layers["counts"].sum(axis=1)).ravel()
    nc_official = np.asarray(official[shared].X.sum(axis=1)).ravel()
    corr_counts = np.corrcoef(nc_cloupe, nc_official)[0, 1]
    print(f"nCount correlation: {corr_counts:.6f}")
    assert corr_counts > 0.99

    pos_official = pd.read_parquet(
        os.path.join(official_dir, "spatial", "tissue_positions.parquet")
    ).set_index("barcode").loc[shared]

    row_diff = np.abs(
        adata[shared].obs["array_row"].to_numpy() - pos_official["array_row"].to_numpy()
    )
    col_diff = np.abs(
        adata[shared].obs["array_col"].to_numpy() - pos_official["array_col"].to_numpy()
    )
    print(f"array_row exact match: {(row_diff == 0).mean():.4f} (max diff {row_diff.max()}), "
          f"array_col exact match: {(col_diff == 0).mean():.4f} (max diff {col_diff.max()})")
    assert (row_diff == 0).mean() > 0.99
    assert (col_diff == 0).mean() > 0.99
    assert row_diff.max() <= 2
    assert col_diff.max() <= 2

    sp = pd.DataFrame(adata[shared].obsm["spatial"], index=shared, columns=["pxl_col", "pxl_row"])
    cor_col = np.corrcoef(sp["pxl_col"], pos_official.loc[shared, "pxl_col_in_fullres"])[0, 1]
    cor_row = np.corrcoef(sp["pxl_row"], pos_official.loc[shared, "pxl_row_in_fullres"])[0, 1]
    cor_row_vs_col = np.corrcoef(sp["pxl_row"], pos_official.loc[shared, "pxl_col_in_fullres"])[0, 1]
    print(f"pxl_col correlation: {cor_col:.6f}, pxl_row correlation: {cor_row:.6f}")
    assert cor_col > 0.999
    assert cor_row > 0.999
    # axis-swap sanity check: row shouldn't correlate with col
    assert abs(cor_row_vs_col) < 0.5
