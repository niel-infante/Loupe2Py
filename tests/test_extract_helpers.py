from loupe2py.extract import check_format_version, parse_array_position


# ---------------------------------------------------------------------------
# parse_array_position
# ---------------------------------------------------------------------------

def test_parse_array_position_parses_hd_barcode_exactly():
    ar, ac = parse_array_position("s_008um_00269_00526-1", pxl_row=9999, pxl_col=9999, bin_px=1)
    assert (ar, ac) == (269, 526)


def test_parse_array_position_ignores_pixel_args_when_barcode_matches():
    # The whole point of barcode parsing is that it's exact regardless of
    # (potentially pixel-origin-shifted) pixel coordinates.
    ar1, ac1 = parse_array_position("s_016um_00001_00002-1", pxl_row=0, pxl_col=0, bin_px=1)
    ar2, ac2 = parse_array_position("s_016um_00001_00002-1", pxl_row=99999, pxl_col=99999, bin_px=1)
    assert (ar1, ac1) == (ar2, ac2) == (1, 2)


def test_parse_array_position_falls_back_to_pixel_rounding_for_non_hd_barcode():
    ar, ac = parse_array_position("cellid_000000001-1", pxl_row=17.6, pxl_col=8.4, bin_px=2.0)
    assert (ar, ac) == (round(17.6 / 2.0), round(8.4 / 2.0))


def test_parse_array_position_returns_zero_zero_without_bin_px_or_hd_barcode():
    ar, ac = parse_array_position("cellid_000000001-1", pxl_row=17.6, pxl_col=8.4, bin_px=None)
    assert (ar, ac) == (0, 0)


# ---------------------------------------------------------------------------
# check_format_version
# ---------------------------------------------------------------------------

class _FakeCloupe:
    def __init__(self, header, index_block):
        self.header = header
        self.index_block = index_block


def test_check_format_version_no_warnings_for_known_good_versions():
    cl = _FakeCloupe(
        header={"version": "9.0.0"},
        index_block={
            "Runs": [{"FormatVersion": "3.0.0"}],
            "Matrices": [{"FormatVersion": "6.3.0"}],
            "Projections": [{"FormatVersion": "4.1.0"}],
        },
    )
    result = check_format_version(cl)
    assert result["warnings"] == []
    assert result["versions"]["container"] == "9.0.0"


def test_check_format_version_warns_on_unrecognized_versions():
    cl = _FakeCloupe(
        header={"version": "10.0.0"},
        index_block={
            "Runs": [{"FormatVersion": "3.0.0"}],
            "Matrices": [{"FormatVersion": "99.9.9"}],
            "Projections": [{"FormatVersion": "4.1.0"}, {"FormatVersion": "5.0.0"}],
        },
    )
    result = check_format_version(cl)
    assert any("container" in w for w in result["warnings"])
    assert any("matrix" in w for w in result["warnings"])
    assert any("projection" in w for w in result["warnings"])
    # the known-good run/one projection version shouldn't be flagged
    assert not any("run" in w for w in result["warnings"])


def test_check_format_version_handles_missing_sections():
    cl = _FakeCloupe(header={"version": "9.0.0"}, index_block={})
    result = check_format_version(cl)
    assert result["versions"]["run"] is None
    assert result["versions"]["matrix"] is None
    assert result["versions"]["projection"] == []
