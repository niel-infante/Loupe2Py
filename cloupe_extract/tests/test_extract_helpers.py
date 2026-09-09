import json
import tempfile

import numpy as np
import pytest
import scipy.sparse

import cloupe_extract.extract as extract_module
from cloupe_extract.extract import (
    UnvalidatedFormatVersionError,
    _cellseg_positions_from_geojson,
    _polygon_centroid,
    check_format_version,
    exclude_synthetic_totals,
    extract_cloupe,
    get_cellseg_projection,
    parse_array_position,
)


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


# ---------------------------------------------------------------------------
# extract_cloupe() -- version_check gate (raises by default, opt out with
# version_check=False; the caller then owns the risk, not the library)
# ---------------------------------------------------------------------------

class _FakeCloupeUnrecognizedVersion(_FakeCloupe):
    """A stand-in for _vendor.cloupe.Cloupe reporting an unvalidated version.

    matrices=[] lets a version_check=False run prove it got *past* the
    version gate (it fails later, on `cl.matrices[0]`, with an ordinary
    IndexError -- not UnvalidatedFormatVersionError) without needing a full
    fake matrix/barcode pipeline just to test this one control-flow branch.
    """

    def __init__(self):
        super().__init__(
            header={"version": "10.0.0"},
            index_block={
                "Runs": [{"FormatVersion": "3.0.0"}],
                "Matrices": [{"FormatVersion": "6.3.0"}],
                "Projections": [{"FormatVersion": "4.1.0"}],
            },
        )
        self.matrices = []


def test_extract_cloupe_raises_by_default_on_unvalidated_version(monkeypatch):
    monkeypatch.setattr(extract_module, "Cloupe", lambda path, load_csr=True: _FakeCloupeUnrecognizedVersion())
    with tempfile.TemporaryDirectory() as outdir:
        with pytest.raises(UnvalidatedFormatVersionError, match="container"):
            extract_cloupe("fake.cloupe", outdir)
        # Provenance survives the failure: the caller can still see exactly
        # what was detected, even though extraction refused to proceed.
        with open(f"{outdir}/format_info.json") as f:
            fmt_info = json.load(f)
        assert any("container" in w for w in fmt_info["warnings"])


def test_extract_cloupe_version_check_false_proceeds_past_the_gate(monkeypatch):
    monkeypatch.setattr(extract_module, "Cloupe", lambda path, load_csr=True: _FakeCloupeUnrecognizedVersion())
    with tempfile.TemporaryDirectory() as outdir:
        with pytest.raises(IndexError):
            extract_cloupe("fake.cloupe", outdir, version_check=False)


# ---------------------------------------------------------------------------
# get_cellseg_projection -- regression test for the diagonal-collapse bug
# (Centers is struct-of-arrays: all x values then all y values, NOT
# interleaved (x,y) pairs -- see the comment above the .reshape(2, n_cells).T
# call in get_cellseg_projection() for the full story).
# ---------------------------------------------------------------------------

class _FakeCloupeCellSegs:
    """Minimal fake supporting only what get_cellseg_projection() touches."""

    def __init__(self, index_block, n_barcodes, blocks, barcodes=None):
        self.index_block = index_block
        barcodes = barcodes or [f"cellid_{i:09d}-1" for i in range(n_barcodes)]
        self.matrices = [{"BarcodeCount": n_barcodes, "Barcodes": barcodes}]
        self._blocks = blocks  # {(start, end): bytes}

    def read_block(self, start, end, as_json=False, verbose=False):
        return self._blocks[(start, end)]


def test_get_cellseg_projection_unscrambles_struct_of_arrays_centers():
    # 4 segments: segments 0,1 -> barcode 0; segments 2,3 -> barcode 1.
    # x = [10, 20, 100, 200], y = [1000, 2000, 5, 15], stored as
    # struct-of-arrays: [x0,x1,x2,x3, y0,y1,y2,y3]. If this were (wrongly)
    # read as interleaved (x,y) pairs instead, segment 0 would come out as
    # (10, 20) rather than the true (10, 1000) -- a completely different,
    # scrambled result. This is exactly the failure mode the real bug
    # produced against real files.
    x = np.array([10.0, 20.0, 100.0, 200.0])
    y = np.array([1000.0, 2000.0, 5.0, 15.0])
    centers_bytes = np.concatenate([x, y]).astype(np.float64).tobytes()
    bi_bytes = np.array([0, 0, 1, 1], dtype=np.int32).tobytes()

    centers_range = (0, len(centers_bytes))
    bi_range = (len(centers_bytes), len(centers_bytes) + len(bi_bytes))

    cl = _FakeCloupeCellSegs(
        index_block={
            "CellSegs": [{
                "MicronsPerPixel": 1.5,  # must NOT affect the result (no /mpp anymore)
                "RectCount": 4,
                "Centers": {"Start": centers_range[0], "End": centers_range[1]},
                "BarcodeIndices": {"Start": bi_range[0], "End": bi_range[1]},
            }]
        },
        n_barcodes=2,
        blocks={centers_range: centers_bytes, bi_range: bi_bytes},
    )

    result = get_cellseg_projection(cl)

    assert result["pxl_col"] == [15.0, 150.0]  # mean(10,20), mean(100,200)
    assert result["pxl_row"] == [1500.0, 10.0]  # mean(1000,2000), mean(5,15)
    assert result["bin_size_px"] is None
    assert result["microns_per_pixel"] == 1.5


# ---------------------------------------------------------------------------
# _polygon_centroid / _cellseg_positions_from_geojson / get_cellseg_projection
# GeoJSON path -- CellSegs' embedded GeoJSON is the exact cell boundary data
# SpaceRanger's own pipeline uses (confirmed byte-identical, coordinate for
# coordinate, to a real paired segmented_outputs/cell_segmentations.geojson).
# The rect-average fallback above is a coarser approximation that only covers
# ~69% of the 2 um bins SpaceRanger itself assigns to a cell on a real
# dataset checked; the true polygon centroid has no such gap and is what
# get_cellseg_projection() should prefer whenever it's available.
# ---------------------------------------------------------------------------

def test_polygon_centroid_is_area_weighted_not_a_vertex_average():
    # L-shaped hexagon: a 4x1 rectangle union a 1x2 rectangle sharing the
    # corner at (1,1). True area centroid (by decomposition into those two
    # rectangles, area-weighted): (1.5, 1.0). A naive vertex average of the
    # 6 corners gives a different point, (1.667, 1.333) -- exactly the kind
    # of bias this function exists to avoid for non-convex/irregular cells.
    ring = [[0, 0], [4, 0], [4, 1], [1, 1], [1, 3], [0, 3]]
    cx, cy = _polygon_centroid(ring)
    assert (round(cx, 6), round(cy, 6)) == (1.5, 1.0)

    naive_x = sum(p[0] for p in ring) / len(ring)
    naive_y = sum(p[1] for p in ring) / len(ring)
    assert (cx, cy) != (naive_x, naive_y)


def test_polygon_centroid_of_a_simple_square():
    ring = [[0, 0], [2, 0], [2, 2], [0, 2]]
    cx, cy = _polygon_centroid(ring)
    assert (cx, cy) == (1.0, 1.0)


def test_polygon_centroid_winding_order_does_not_flip_sign():
    # Same square, vertices listed clockwise instead of counter-clockwise --
    # the centroid must come out identical either way (numerator and
    # denominator both carry the signed area and cancel).
    ring = [[0, 0], [0, 2], [2, 2], [2, 0]]
    cx, cy = _polygon_centroid(ring)
    assert (cx, cy) == (1.0, 1.0)


def _fake_geojson_block(features):
    body = json.dumps({"type": "FeatureCollection", "features": features}).encode()
    return {"Start": 0, "End": len(body)}, {(0, len(body)): body}


def test_cellseg_positions_from_geojson_uses_true_centroid_and_matches_barcode():
    # Two cells: a plain square (cell_id 5) and the L-shape from above
    # (cell_id 12), matched to barcodes by the numeric id embedded in each
    # barcode string, not by list position.
    features = [
        {"type": "Feature", "geometry": {"type": "Polygon",
            "coordinates": [[[0, 0], [2, 0], [2, 2], [0, 2]]]},
         "properties": {"cell_id": 5}},
        {"type": "Feature", "geometry": {"type": "Polygon",
            "coordinates": [[[0, 0], [4, 0], [4, 1], [1, 1], [1, 3], [0, 3]]]},
         "properties": {"cell_id": 12}},
    ]
    geojson_field, blocks = _fake_geojson_block(features)

    class _Fake:
        def read_block(self, start, end, as_json=False, verbose=False):
            return blocks[(start, end)]

    barcodes = ["cellid_000000012-1", "cellid_000000005-1"]  # deliberately out of order
    pxl_col, pxl_row, found = _cellseg_positions_from_geojson(
        _Fake(), geojson_field, barcodes, n_barcodes=2
    )

    assert found.tolist() == [True, True]
    assert (round(pxl_col[0], 6), round(pxl_row[0], 6)) == (1.5, 1.0)  # barcode 0 -> cell 12 (L-shape)
    assert (pxl_col[1], pxl_row[1]) == (1.0, 1.0)                       # barcode 1 -> cell 5 (square)


def test_get_cellseg_projection_prefers_geojson_over_rects_when_both_present():
    # cell_id 5's rect-average (below) would give a biased position; the
    # GeoJSON polygon centroid should win instead.
    features = [{"type": "Feature", "geometry": {"type": "Polygon",
                 "coordinates": [[[0, 0], [2, 0], [2, 2], [0, 2]]]},
                 "properties": {"cell_id": 5}}]
    geojson_field, geojson_blocks = _fake_geojson_block(features)

    # A single rect for the same barcode centered somewhere else entirely --
    # if the fallback fired, the result would be (100.0, 100.0), not (1, 1).
    centers_bytes = np.array([100.0, 100.0], dtype=np.float64).tobytes()  # x=[100], y=[100]
    bi_bytes = np.array([0], dtype=np.int32).tobytes()
    centers_range = (1000, 1000 + len(centers_bytes))
    bi_range = (2000, 2000 + len(bi_bytes))

    blocks = dict(geojson_blocks)
    blocks[centers_range] = centers_bytes
    blocks[bi_range] = bi_bytes

    cl = _FakeCloupeCellSegs(
        index_block={
            "CellSegs": [{
                "MicronsPerPixel": 1.0,
                "RectCount": 1,
                "Centers": {"Start": centers_range[0], "End": centers_range[1]},
                "BarcodeIndices": {"Start": bi_range[0], "End": bi_range[1]},
                "GeoJSON": geojson_field,
            }]
        },
        n_barcodes=1,
        blocks=blocks,
        barcodes=["cellid_000000005-1"],
    )

    result = get_cellseg_projection(cl)
    assert (result["pxl_col"][0], result["pxl_row"][0]) == (1.0, 1.0)


def test_get_cellseg_projection_falls_back_to_rects_for_barcodes_missing_from_geojson():
    # cell_id 5 has a GeoJSON polygon; barcode 1 (cell_id 7) does not, and
    # must fall back to its rect-average instead of being left as (0, 0).
    features = [{"type": "Feature", "geometry": {"type": "Polygon",
                 "coordinates": [[[0, 0], [2, 0], [2, 2], [0, 2]]]},
                 "properties": {"cell_id": 5}}]
    geojson_field, geojson_blocks = _fake_geojson_block(features)

    # rect 0 -> barcode 0 (cell 5, ignored since GeoJSON covers it);
    # rect 1 -> barcode 1 (cell 7, no GeoJSON entry, must use this rect).
    x = np.array([9999.0, 50.0])
    y = np.array([9999.0, 60.0])
    centers_bytes = np.concatenate([x, y]).astype(np.float64).tobytes()
    bi_bytes = np.array([0, 1], dtype=np.int32).tobytes()
    centers_range = (1000, 1000 + len(centers_bytes))
    bi_range = (2000, 2000 + len(bi_bytes))

    blocks = dict(geojson_blocks)
    blocks[centers_range] = centers_bytes
    blocks[bi_range] = bi_bytes

    cl = _FakeCloupeCellSegs(
        index_block={
            "CellSegs": [{
                "MicronsPerPixel": 1.0,
                "RectCount": 2,
                "Centers": {"Start": centers_range[0], "End": centers_range[1]},
                "BarcodeIndices": {"Start": bi_range[0], "End": bi_range[1]},
                "GeoJSON": geojson_field,
            }]
        },
        n_barcodes=2,
        blocks=blocks,
        barcodes=["cellid_000000005-1", "cellid_000000007-1"],
    )

    result = get_cellseg_projection(cl)
    assert (result["pxl_col"][0], result["pxl_row"][0]) == (1.0, 1.0)   # from GeoJSON
    assert (result["pxl_col"][1], result["pxl_row"][1]) == (50.0, 60.0)  # from rect fallback


# ---------------------------------------------------------------------------
# exclude_synthetic_totals -- regression test for the "3x inflated counts"
# bug: every .cloupe file's Matrices section appends "type_sum_<FeatureType>"
# and "genome_sum_<Reference>" rows holding each barcode's FULL total, which
# extract_cloupe() previously wrote out as if they were real genes.
# ---------------------------------------------------------------------------

def test_exclude_synthetic_totals_drops_type_and_genome_sum_rows():
    feature_ids = ["ENSG00000001", "ENSG00000002", "type_sum_Gene Expression", "genome_sum_GRCh38"]
    feature_names = ["GeneA", "GeneB", "Gene Expression Sum", "GRCh38 Sum"]
    # barcode 0: GeneA=1, GeneB=2, sums=3,3 ; barcode 1: GeneA=4, GeneB=0, sums=4,4
    csr = scipy.sparse.csr_matrix(np.array([
        [1, 4],
        [2, 0],
        [3, 4],
        [3, 4],
    ]))

    ids, names, filtered_csr, excluded = exclude_synthetic_totals(feature_ids, feature_names, csr)

    assert ids == ["ENSG00000001", "ENSG00000002"]
    assert names == ["GeneA", "GeneB"]
    assert excluded == ["type_sum_Gene Expression", "genome_sum_GRCh38"]
    assert np.array_equal(filtered_csr.toarray(), [[1, 4], [2, 0]])
    # the real per-barcode total no longer includes the synthetic rows
    assert np.asarray(filtered_csr.sum(axis=0)).ravel().tolist() == [3, 4]


def test_exclude_synthetic_totals_is_a_no_op_when_nothing_to_drop():
    feature_ids = ["ENSG00000001", "ENSG00000002"]
    feature_names = ["GeneA", "GeneB"]
    csr = scipy.sparse.csr_matrix(np.array([[1, 4], [2, 0]]))

    ids, names, filtered_csr, excluded = exclude_synthetic_totals(feature_ids, feature_names, csr)

    assert ids == feature_ids
    assert names == feature_names
    assert excluded == []
    assert np.array_equal(filtered_csr.toarray(), csr.toarray())


def test_get_cellseg_projection_drops_unassigned_segments():
    # Segment 1 has bi = -1 (unassigned); should be excluded from the
    # average for barcode 0, not crash np.bincount().
    x = np.array([10.0, 9999.0, 100.0])
    y = np.array([1000.0, 9999.0, 5.0])
    centers_bytes = np.concatenate([x, y]).astype(np.float64).tobytes()
    bi_bytes = np.array([0, -1, 1], dtype=np.int32).tobytes()

    centers_range = (0, len(centers_bytes))
    bi_range = (len(centers_bytes), len(centers_bytes) + len(bi_bytes))

    cl = _FakeCloupeCellSegs(
        index_block={
            "CellSegs": [{
                "MicronsPerPixel": 1.0,
                "RectCount": 3,
                "Centers": {"Start": centers_range[0], "End": centers_range[1]},
                "BarcodeIndices": {"Start": bi_range[0], "End": bi_range[1]},
            }]
        },
        n_barcodes=2,
        blocks={centers_range: centers_bytes, bi_range: bi_bytes},
    )

    result = get_cellseg_projection(cl)
    assert result["pxl_col"] == [10.0, 100.0]
    assert result["pxl_row"] == [1000.0, 5.0]
