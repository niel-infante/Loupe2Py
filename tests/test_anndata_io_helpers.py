import pandas as pd

from loupe2py.anndata_io import _align_positions, _derive_sample_name, _detect_mt_pattern


def test_detect_mt_pattern_human():
    assert _detect_mt_pattern(["MT-ND1", "GAPDH", "MT-CO1"]) == "MT-"


def test_detect_mt_pattern_mouse():
    assert _detect_mt_pattern(["mt-Nd1", "Gapdh", "mt-Co1"]) == "mt-"


def test_detect_mt_pattern_fallback():
    # No MT- genes -> falls through to "mt-", even if that also matches
    # nothing (pct_counts_mt ends up 0 either way). Matches the R side.
    assert _detect_mt_pattern(["GAPDH", "ACTB"]) == "mt-"


def test_detect_mt_pattern_prefers_human():
    assert _detect_mt_pattern(["MT-ND1", "mt-Nd1"]) == "MT-"


def test_derive_sample_name_strips_dir_and_final_extension():
    assert _derive_sample_name("/path/to/sample1.cloupe") == "sample1"


def test_derive_sample_name_only_strips_final_extension():
    assert _derive_sample_name("sample.v2.cloupe") == "sample.v2"


def test_derive_sample_name_no_extension():
    assert _derive_sample_name("sample_no_ext") == "sample_no_ext"


def test_align_positions_reindexes_to_requested_order():
    pos = pd.DataFrame({"barcode": ["bc2", "bc1", "bc3"], "value": [20, 10, 30]})
    aligned = _align_positions(pos, ["bc1", "bc2", "bc3"])
    assert list(aligned.index) == ["bc1", "bc2", "bc3"]
    assert list(aligned["value"]) == [10, 20, 30]


def test_align_positions_nan_for_missing_names():
    pos = pd.DataFrame({"barcode": ["bc1", "bc2"], "value": [10, 20]})
    aligned = _align_positions(pos, ["bc1", "bc_missing"])
    assert aligned["value"].iloc[0] == 10
    assert pd.isna(aligned["value"].iloc[1])
