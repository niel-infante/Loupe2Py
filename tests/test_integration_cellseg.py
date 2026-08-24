# Opt-in: regression guard for the cell-segmentation diagonal-collapse bug
# (Centers read as interleaved (x,y) pairs instead of struct-of-arrays,
# producing pxl_row ~= pxl_col for 30-40% of cells in every real
# cell-segmented .cloupe file checked). There's no official ground-truth
# output to cross-validate cell-segmentation positions against the way HD
# binned data has (see test_integration_visium_hd.py), so this test asserts
# the same diagnostic used to find and confirm the fix: real 2D tissue
# positions should not show a near-perfect x~=y correlation.
#
# Not run by default (no real .cloupe file is committed to this repo) --
# set LOUPE2PY_CELLSEG_TEST_FILE to a real cell-segmentation .cloupe file to
# run it, e.g.:
#   LOUPE2PY_CELLSEG_TEST_FILE="/path/to/sample_cell.cloupe" pytest tests/test_integration_cellseg.py

import os

import numpy as np
import pandas as pd
import pytest

test_file = os.environ.get("LOUPE2PY_CELLSEG_TEST_FILE", "")
pytestmark = pytest.mark.skipif(
    not test_file, reason="Set LOUPE2PY_CELLSEG_TEST_FILE to run the cell-segmentation regression test"
)


def test_cellseg_positions_are_not_scrambled_onto_a_diagonal(tmp_path):
    from loupe2py import extract_cloupe

    if not os.path.exists(test_file):
        pytest.skip(f"{test_file} not found")

    extract_cloupe(test_file, str(tmp_path), include_image=False)
    pos = pd.read_csv(tmp_path / "tissue_positions.csv")

    row = pos["pxl_row_in_fullres"].to_numpy()
    col = pos["pxl_col_in_fullres"].to_numpy()

    corr = np.corrcoef(row, col)[0, 1]
    coeffs = np.polyfit(col, row, 1)
    resid = np.abs(row - np.polyval(coeffs, col))
    diag_frac = (resid < 5).mean()

    print(f"n={len(pos)}  corr(pxl_row, pxl_col)={corr:.4f}  near-diagonal fraction={diag_frac:.4f}")

    # Real 2D tissue positions should not correlate this strongly; a bug
    # that scrambles Centers onto (roughly) the diagonal produces |corr| well
    # above 0.9 and a near-diagonal fraction in the 30-40% range (see the
    # bug report this test guards against). Healthy real data observed so
    # far tops out well under these thresholds.
    assert abs(corr) < 0.5
    assert diag_frac < 0.05

    # array_row/array_col are unconditionally 0 for cell-segmentation data
    # by design (no fixed grid exists) -- confirm that's still true rather
    # than silently changing, since downstream code/docs assume it.
    assert (pos["array_row"] == 0).all()
    assert (pos["array_col"] == 0).all()
