"""
Extract data from 10x Genomics .cloupe files into a set of generic,
SpaceRanger-convention-like intermediate files.

Usage (library):
    from cloupe_extract.extract import extract_cloupe
    extract_cloupe("sample.cloupe", "/path/to/outdir")

Usage (CLI):
    python -m cloupe_extract.extract <cloupe_path> <outdir>

Outputs written to outdir:
    matrix.mtx.gz           - sparse count matrix (Market Exchange format)
    barcodes.tsv.gz          - barcode list (one per line)
    features.tsv.gz          - feature table (id, name, type)
    tissue_positions.csv     - spatial coordinates
    scalefactors_json.json   - scale factors for image alignment
    tissue_hires_image.png   - stitched tissue image (if tiles present)
    projections.csv          - UMAP and other embeddings
    clusterings.csv          - Spaceranger graph and k-means cluster labels
    celltracks.csv           - user-created annotations (if any)
    format_info.json         - detected .cloupe format versions + any
                                unvalidated-version warnings

Requires: scipy, numpy, Pillow. The .cloupe binary parser (cellgeni/cloupe)
is vendored in cloupe_extract._vendor.cloupe -- no external install or path
configuration needed.

This module is deliberately framework-agnostic: it has no knowledge of
Seurat, AnnData, or any other downstream object model. loupe2py (Python,
AnnData) and Loupe2R (R, Seurat, via reticulate) both depend on this
package and build on top of it by reading the files it writes.
"""

import gzip
import io
import json
import os
import re
import struct
import sys

import numpy as np
import scipy.io
import scipy.sparse

try:
    from PIL import Image
    _HAS_PILLOW = True
except ImportError:
    _HAS_PILLOW = False

from ._vendor.cloupe import Cloupe


# ---------------------------------------------------------------------------
# Format version guard
# ---------------------------------------------------------------------------

# Known-good format versions, seeded from real Visium HD .cloupe files in
# both binned and cell-segmentation modes (multiple samples/labs; same
# values validated for Loupe2R). Grow these sets as new files are
# validated; an unrecognized version triggers a warning, not a failure,
# since the parser may well still be correct -- it just hasn't been checked.
_TESTED_VERSIONS = {
    "container": {"9.0.0", "8.0.0"},
    "run": {"3.0.0"},
    "matrix": {"6.3.0", "6.2.0"},
    "projection": {"4.1.0"},
}


def check_format_version(cloupe_obj):
    """Compare this file's internal format versions against the known-tested set.

    Returns a dict with 'versions' (what was found) and 'warnings' (a list of
    human-readable strings for anything outside _TESTED_VERSIONS).
    """
    runs = cloupe_obj.index_block.get("Runs", [])
    matrices = cloupe_obj.index_block.get("Matrices", [])
    projections = cloupe_obj.index_block.get("Projections", [])

    versions = {
        "container": cloupe_obj.header.get("version"),
        "run": runs[0].get("FormatVersion") if runs else None,
        "matrix": matrices[0].get("FormatVersion") if matrices else None,
        "projection": sorted({
            p.get("FormatVersion") for p in projections if p.get("FormatVersion")
        }),
    }

    warnings = []
    for key in ("container", "run", "matrix"):
        v = versions.get(key)
        if v is not None and v not in _TESTED_VERSIONS[key]:
            warnings.append(
                f"{key} format version '{v}' has not been validated "
                f"(tested: {sorted(_TESTED_VERSIONS[key])})"
            )
    for v in versions.get("projection", []):
        if v not in _TESTED_VERSIONS["projection"]:
            warnings.append(
                f"projection format version '{v}' has not been validated "
                f"(tested: {sorted(_TESTED_VERSIONS['projection'])})"
            )

    return {"versions": versions, "warnings": warnings}


# ---------------------------------------------------------------------------
# Image reconstruction
# ---------------------------------------------------------------------------

def stitch_tiles(cloupe_obj, target_level=None):
    """Reconstruct full tissue image from the cloupe tile pyramid.

    Returns a PIL.Image, or None if no tiles are present or Pillow is missing.
    """
    if not _HAS_PILLOW:
        print("[cloupe] Pillow not installed - skipping image extraction "
              "(pip install Pillow to enable)", flush=True)
        return None

    tile_sets = cloupe_obj.index_block.get("SpatialImageTiles", [])
    if not tile_sets:
        return None

    ts = tile_sets[0]
    tiles = ts.get("Tiles", {})
    if not tiles:
        return None

    img_w, img_h = ts["Dims"]
    tile_size = ts["TileSize"]
    tile_overlap = ts.get("TileOverlap", 0)
    stride = tile_size - tile_overlap

    # Auto-select highest zoom level (most tiles = highest detail)
    levels = sorted({int(k.split("/")[0]) for k in tiles})
    if target_level is None:
        target_level = max(levels)

    level_tiles = {k: v for k, v in tiles.items()
                   if k.startswith(f"{target_level}/")}

    # Detect image mode from first tile
    first_key = sorted(level_tiles)[0]
    first_bytes = cloupe_obj.read_block(
        level_tiles[first_key]["Start"], level_tiles[first_key]["End"]
    )
    first_img = Image.open(io.BytesIO(first_bytes))
    mode = first_img.mode

    canvas = Image.new(mode, (img_w, img_h))

    for tile_key, tile_info in sorted(level_tiles.items()):
        _, col_row = tile_key.split("/")
        col_s, row_s = col_row.replace(".png", "").split("_")
        col, row = int(col_s), int(row_s)

        x0 = col * stride
        y0 = row * stride

        tile_bytes = cloupe_obj.read_block(tile_info["Start"], tile_info["End"])
        tile_img = Image.open(io.BytesIO(tile_bytes))

        # Crop tile to canvas bounds (last tile may extend beyond image edge)
        x1 = min(x0 + tile_img.width, img_w)
        y1 = min(y0 + tile_img.height, img_h)
        tile_crop = tile_img.crop((0, 0, x1 - x0, y1 - y0))
        canvas.paste(tile_crop, (x0, y0))

    # Tile mode is auto-detected from the source tiles above and could be
    # palette ("P") or grayscale ("L"); force a standard 3-channel RGB image
    # so downstream consumers always get a consistent array shape regardless
    # of the source tile mode.
    return canvas.convert("RGB")


# ---------------------------------------------------------------------------
# Spatial coordinates
# ---------------------------------------------------------------------------

def _polygon_centroid(ring):
    """True area centroid of a simple polygon ring via the shoelace formula.

    ring: sequence of [x, y] pairs (GeoJSON exterior ring; closed or not --
    both work, since a redundant closing vertex contributes a zero-length
    edge). Falls back to a plain vertex average for degenerate (zero-area)
    rings, e.g. a 1-2 point sliver.
    """
    pts = np.asarray(ring, dtype=np.float64)
    x, y = pts[:, 0], pts[:, 1]
    x2, y2 = np.roll(x, -1), np.roll(y, -1)
    cross = x * y2 - x2 * y
    area2 = cross.sum()
    if abs(area2) < 1e-9:
        return float(x.mean()), float(y.mean())
    cx = ((x + x2) * cross).sum() / (3 * area2)
    cy = ((y + y2) * cross).sum() / (3 * area2)
    return float(cx), float(cy)


def _cellseg_positions_from_geojson(cloupe_obj, geojson_field, barcodes, n_barcodes):
    """Per-cell true area centroids from CellSegs' embedded GeoJSON boundaries.

    This is the exact, complete cell-boundary data SpaceRanger's own official
    pipeline uses -- confirmed byte-identical (coordinate-for-coordinate) to
    real paired segmented_outputs/cell_segmentations.geojson on a real Visium
    HD sample. Distinct from the coarser Centers/Sizes rect decomposition
    below: rects are a compact, lossy tiling of each cell's mask (~2.1M rects
    for ~148k cells in that sample) that only covers ~69% of the 2 um bins
    SpaceRanger itself assigns to a cell; a naive rect-center average is
    biased by that missing coverage. The polygon boundary has no such gap.

    Returns (pxl_col, pxl_row, found) as float arrays of length n_barcodes;
    found[i] is False where barcode i has no matching polygon (caller should
    fall back to the rect-based estimate for those).
    """
    raw = cloupe_obj.read_block(geojson_field["Start"], geojson_field["End"])
    gj = json.loads(raw)

    # Map each barcode's own encoded cell id -> its index, parsed from the
    # barcode string itself rather than reconstructed, so this doesn't
    # depend on assuming a fixed zero-pad width or "-1" suffix.
    cellid_to_idx = {}
    for i, bc in enumerate(barcodes):
        m = re.match(r"^cellid_(\d+)-\d+$", bc)
        if m:
            cellid_to_idx[int(m.group(1))] = i

    pxl_col = np.full(n_barcodes, np.nan)
    pxl_row = np.full(n_barcodes, np.nan)
    found = np.zeros(n_barcodes, dtype=bool)

    for feat in gj.get("features", []):
        if feat.get("geometry", {}).get("type") != "Polygon":
            continue
        idx = cellid_to_idx.get(feat.get("properties", {}).get("cell_id"))
        if idx is None:
            continue
        cx, cy = _polygon_centroid(feat["geometry"]["coordinates"][0])
        pxl_col[idx] = cx
        pxl_row[idx] = cy
        found[idx] = True

    return pxl_col, pxl_row, found


def _cellseg_positions_from_rects(cloupe_obj, cs, n_barcodes):
    """Per-cell position as the unweighted average of CellSegs rect centers.

    Coarser fallback for files/barcodes without GeoJSON coverage -- see
    get_cellseg_projection()'s docstring for why GeoJSON is preferred.

    Returns (pxl_col, pxl_row, found) as float arrays of length n_barcodes.
    """
    n_rects = cs["RectCount"]

    raw_c = cloupe_obj.read_block(cs["Centers"]["Start"], cs["Centers"]["End"])
    # Struct-of-arrays layout: [x0, x1, ..., xN-1, y0, y1, ..., yN-1], NOT
    # interleaved (x0, y0, x1, y1, ...). This matches the same convention the
    # Cloupe class's own Projections parser already uses elsewhere in this
    # file format (see how "coordinates" is built in the vendored cloupe.py:
    # each axis's full run of values comes first, then the next axis's).
    # An earlier version of this function used .reshape(n_cells, 2), which
    # silently paired up two nearby same-axis values as a fake (x, y) pair
    # for a large fraction of segments -- producing near-perfect
    # pxl_row=pxl_col diagonal artifacts for ~30-40% of cells in every real
    # cell-segmented .cloupe file checked. Confirmed by an independent
    # signal too: per-segment corr(x, y) is ~0.99 under the wrong reshape
    # (implausible for real 2D tissue positions) and drops to a normal,
    # weak ~-0.15 under this one.
    centers = np.frombuffer(raw_c, dtype=np.float64).reshape(2, n_rects).T
    # The CellSegs "Sizes" field (per-rect width/height) shares this same
    # struct-of-arrays layout; used by the GeoJSON path above, not here.

    raw_bi = cloupe_obj.read_block(
        cs["BarcodeIndices"]["Start"], cs["BarcodeIndices"]["End"]
    )
    bi = np.frombuffer(raw_bi, dtype=np.int32)

    # -1 (or other negative values) marks an unassigned cell segment in 10x's
    # index-array convention; np.bincount() raises on negative input, so drop
    # those segments before centroid-averaging.
    assigned = bi >= 0
    bi, centers = bi[assigned], centers[assigned]

    counts = np.bincount(bi, minlength=n_barcodes)
    safe_cnt = np.maximum(counts, 1)
    pxl_col = np.bincount(bi, weights=centers[:, 0], minlength=n_barcodes) / safe_cnt
    pxl_row = np.bincount(bi, weights=centers[:, 1], minlength=n_barcodes) / safe_cnt
    found = counts > 0

    return pxl_col, pxl_row, found


def get_cellseg_projection(cloupe_obj):
    """Return spatial coordinates from CellSegs (Visium HD, cell-segmentation mode --
    a SpaceRanger 4.0+ analysis mode, distinct from binned mode, that assigns
    transcripts to individual cells via image-based nucleus/cell segmentation;
    each mode produces its own separate .cloupe file).

    Prefers the embedded GeoJSON polygon boundaries (see
    _cellseg_positions_from_geojson()) -- the exact cell shapes, not an
    approximation -- and computes each cell's true area centroid from them.
    Falls back to averaging the coarser Centers/Sizes rect decomposition
    (_cellseg_positions_from_rects()) for any barcode GeoJSON doesn't cover,
    or for the whole file if it has no GeoJSON field at all (older/unusual
    files; not observed in any real file checked so far, but not assumed
    universal either).

    Returns dict with keys: pxl_col, pxl_row, bin_size_px, microns_per_pixel.
    """
    cell_segs = cloupe_obj.index_block.get("CellSegs", [])
    if not cell_segs:
        return None

    cs = cell_segs[0]
    mpp = cs.get("MicronsPerPixel")
    if not mpp:
        return None

    barcodes = cloupe_obj.matrices[0]["Barcodes"]
    n_barcodes = cloupe_obj.matrices[0]["BarcodeCount"]

    pxl_col = np.full(n_barcodes, np.nan)
    pxl_row = np.full(n_barcodes, np.nan)
    found = np.zeros(n_barcodes, dtype=bool)

    geojson_field = cs.get("GeoJSON")
    if geojson_field is not None:
        pxl_col, pxl_row, found = _cellseg_positions_from_geojson(
            cloupe_obj, geojson_field, barcodes, n_barcodes
        )

    if not found.all():
        rect_col, rect_row, rect_found = _cellseg_positions_from_rects(
            cloupe_obj, cs, n_barcodes
        )
        need = ~found
        pxl_col[need] = rect_col[need]
        pxl_row[need] = rect_row[need]
        found = found | rect_found

    # Any barcode with no position from either source (never observed, but
    # not structurally impossible) reports (0, 0) rather than NaN, matching
    # this function's long-standing contract for "no data" cells.
    pxl_col = np.nan_to_num(pxl_col, nan=0.0)
    pxl_row = np.nan_to_num(pxl_row, nan=0.0)

    return {
        "pxl_col": pxl_col.tolist(),
        "pxl_row": pxl_row.tolist(),
        "bin_size_px": None,
        "microns_per_pixel": mpp,
    }


def get_spatial_projection(cloupe_obj):
    """Return spatial coordinate arrays and scale metadata.

    Tries the Visium HD 'Spatial' projection (binned mode) first, then falls
    back to CellSegs coordinates (cell-segmentation mode).

    Returns dict with keys: pxl_col, pxl_row, bin_size_px, microns_per_pixel,
    or None if no spatial information is found.
    All lists are aligned to the barcode order in cloupe_obj.matrices[0].
    """
    mpp = None
    for proj in cloupe_obj.index_block.get("Projections", []):
        if proj["Name"] == "Spatial":
            mpp = proj.get("MicronsPerPixel")
            break

    spatial = next(
        (p for p in cloupe_obj.projections if p["name"] == "Spatial"), None
    )
    if spatial is not None:
        coords = spatial["coordinates"]
        return {
            "pxl_col": coords[0],
            "pxl_row": coords[1],
            "bin_size_px": coords[2][0] if len(coords) > 2 else None,
            "microns_per_pixel": mpp,
        }

    return get_cellseg_projection(cloupe_obj)


# ---------------------------------------------------------------------------
# array_row / array_col
# ---------------------------------------------------------------------------

# Visium HD barcodes encode their own grid position, e.g. "s_008um_00269_00526-1"
# -> row 269, col 526. Parsing this directly is exact (verified against a real
# paired SpaceRanger sample: 100% agreement across 472,752 spots), unlike
# round(pixel / bin_size_px), which was found to disagree with SpaceRanger's
# own array_row/array_col by several hundred bins -- SpaceRanger's grid isn't
# anchored at this package's stitched-image pixel origin. Non-HD barcodes
# (e.g. cell-segmentation "cellid_..." barcodes) don't match this pattern.
# For cell-segmentation data specifically, there IS no pixel-rounding
# fallback in practice: get_cellseg_projection() always returns
# bin_size_px=None (cell segments aren't laid out on any fixed grid), so
# array_row/array_col are unconditionally (0, 0) for every cell in a
# cell-segmented .cloupe file -- not an approximation, just not meaningful
# for this data type. Real spatial positions (pxl_row/pxl_col, obsm['spatial'])
# are unaffected and are what actually drives spot placement in plots.
_HD_BARCODE_RE = re.compile(r"^s_\d+um_(\d+)_(\d+)-\d+$")


def parse_array_position(barcode, pxl_row, pxl_col, bin_px):
    """Return (array_row, array_col) for one barcode.

    Prefers parsing the row/col directly out of a Visium HD barcode name;
    falls back to round(pixel / bin_size_px) for barcodes that aren't in
    that format (e.g. cell-segmentation data), or (0, 0) if neither is
    available.
    """
    m = _HD_BARCODE_RE.match(barcode)
    if m:
        return int(m.group(1)), int(m.group(2))
    if bin_px:
        return int(round(pxl_row / bin_px)), int(round(pxl_col / bin_px))
    return 0, 0


# ---------------------------------------------------------------------------
# Clusterings
# ---------------------------------------------------------------------------

def read_clusterings(cloupe_obj):
    """Return dict: clustering_name -> list of cluster labels (str or 'Unassigned')."""
    result = {}
    for cl in cloupe_obj.index_block.get("Clusterings", []):
        name = cl.get("Name", "Clustering")
        assign_info = cl.get("Assignments")
        if not assign_info:
            continue
        try:
            raw = cloupe_obj.read_block(
                assign_info["Start"], assign_info["End"]
            )
            n = assign_info["ArraySize"]
            keys = struct.unpack(f"{n}h", raw)

            meta = cloupe_obj.read_block(
                cl["Metadata"]["Start"], cl["Metadata"]["End"], as_json=True
            )
            groups = meta.get("groups", [])

            labels = [
                groups[k] if (k >= 0 and k < len(groups)) else "Unassigned"
                for k in keys
            ]
            result[name] = labels
        except Exception as exc:
            print(f"[cloupe] Warning: could not read clustering '{name}': {exc}",
                  flush=True)
    return result


# ---------------------------------------------------------------------------
# Synthetic total rows
# ---------------------------------------------------------------------------

def exclude_synthetic_totals(feature_ids, feature_names, csr):
    """Drop the .cloupe format's built-in per-file total rows.

    Every .cloupe file's Matrices section appends synthetic rows after the
    real features -- one "type_sum_<FeatureType>" row per feature type and
    one "genome_sum_<Reference>" row per reference genome -- each holding
    that barcode's FULL total for the category, not a real gene's count.
    Left in, these silently inflate every downstream total (nCount /
    percent.mt) and show up as bogus "genes" (e.g. "type_sum_Gene
    Expression"). Verified against paired official SpaceRanger output: once
    these rows are excluded and the rest ID-aligned, every real (gene,
    barcode) count matches exactly -- this is a pure extra-rows issue, not a
    value-corruption one.

    Returns (feature_ids, feature_names, csr, excluded_ids) -- all three
    inputs filtered (or returned as-is, with excluded_ids == [], when there
    is nothing to drop).
    """
    real_mask = np.array([
        not (fid.startswith("type_sum_") or fid.startswith("genome_sum_"))
        for fid in feature_ids
    ])
    if real_mask.all():
        return feature_ids, feature_names, csr, []

    excluded = [fid for fid, keep in zip(feature_ids, real_mask) if not keep]
    feature_ids = [fid for fid, keep in zip(feature_ids, real_mask) if keep]
    feature_names = [fn for fn, keep in zip(feature_names, real_mask) if keep]
    return feature_ids, feature_names, csr[real_mask, :], excluded


# ---------------------------------------------------------------------------
# Main extraction function
# ---------------------------------------------------------------------------

def extract_cloupe(cloupe_path, outdir, include_image=True):
    """Extract everything from a .cloupe file into outdir.

    Parameters
    ----------
    cloupe_path : str
        Path to the .cloupe file.
    outdir : str
        Output directory (created if it does not exist).
    include_image : bool
        Whether to stitch and save the tissue image.

    Returns
    -------
    str
        Path to outdir.
    """
    import logging
    logging.getLogger().setLevel(logging.WARNING)

    os.makedirs(outdir, exist_ok=True)

    print(f"[cloupe] Loading {os.path.basename(cloupe_path)} ...", flush=True)
    cl = Cloupe(cloupe_path, load_csr=True)

    # Check internal format versions against the known-tested set before
    # trusting anything else in the file; write the result for callers to
    # surface as a warning (and to stash as object provenance) regardless.
    fmt_info = check_format_version(cl)
    for w in fmt_info["warnings"]:
        print(f"[cloupe] WARNING: {w}", flush=True)
    with open(os.path.join(outdir, "format_info.json"), "w") as f:
        json.dump(fmt_info, f, indent=2)

    matrix = cl.matrices[0]
    barcodes = matrix["Barcodes"]
    feature_ids = matrix["FeatureIds"]
    feature_names = matrix["FeatureNames"]
    n_features = matrix["FeatureCount"]
    n_barcodes = matrix["BarcodeCount"]

    # Structural sanity checks: catch a mis-parsed file (wrong offsets from an
    # unexpected format revision) as an immediate, specific error rather than
    # a downstream object that looks fine until someone notices the biology
    # is wrong.
    assert len(feature_ids) == n_features, \
        f"FeatureCount ({n_features}) != len(FeatureIds) ({len(feature_ids)})"
    assert len(feature_names) == n_features, \
        f"FeatureCount ({n_features}) != len(FeatureNames) ({len(feature_names)})"
    assert len(barcodes) == n_barcodes, \
        f"BarcodeCount ({n_barcodes}) != len(Barcodes) ({len(barcodes)})"
    assert matrix["CSR"].shape == (n_features, n_barcodes), \
        f"CSR shape {matrix['CSR'].shape} != (FeatureCount, BarcodeCount) ({n_features}, {n_barcodes})"

    feature_ids, feature_names, csr_real, excluded = exclude_synthetic_totals(
        feature_ids, feature_names, matrix["CSR"]
    )
    if excluded:
        print(f"[cloupe] Excluding {len(excluded)} synthetic total row(s): "
              f"{excluded}", flush=True)
        n_features = len(feature_ids)

    print(f"[cloupe] {n_features:,} features x {n_barcodes:,} barcodes", flush=True)

    # Count matrix (Features x Barcodes -> int32)
    csr_int = csr_real.astype(np.int32)
    mtx_path = os.path.join(outdir, "matrix.mtx")
    scipy.io.mmwrite(mtx_path, csr_int)
    with open(mtx_path, "rb") as fin, \
         gzip.open(mtx_path + ".gz", "wb") as fout:
        fout.write(fin.read())
    os.remove(mtx_path)

    # Barcodes
    with gzip.open(os.path.join(outdir, "barcodes.tsv.gz"), "wt") as f:
        f.write("\n".join(barcodes) + "\n")

    # Features
    with gzip.open(os.path.join(outdir, "features.tsv.gz"), "wt") as f:
        for fid, fname in zip(feature_ids, feature_names):
            f.write(f"{fid}\t{fname}\tGene Expression\n")

    # Spatial coordinates -> tissue_positions.csv
    spatial = get_spatial_projection(cl)
    if spatial:
        assert len(spatial["pxl_col"]) == n_barcodes and len(spatial["pxl_row"]) == n_barcodes, \
            "Spatial projection coordinate count does not match BarcodeCount"
    mpp = spatial["microns_per_pixel"] if spatial else None
    bin_px = spatial["bin_size_px"] if spatial else None

    with open(os.path.join(outdir, "tissue_positions.csv"), "w") as f:
        f.write("barcode,in_tissue,array_row,array_col,"
                "pxl_row_in_fullres,pxl_col_in_fullres\n")
        if spatial:
            pxl_col = spatial["pxl_col"]
            pxl_row = spatial["pxl_row"]
            for i, bc in enumerate(barcodes):
                pc = pxl_col[i]
                pr = pxl_row[i]
                ar, ac = parse_array_position(bc, pr, pc, bin_px)
                f.write(f"{bc},1,{ar},{ac},{pr:.4f},{pc:.4f}\n")
        else:
            for bc in barcodes:
                f.write(f"{bc},1,0,0,0.0,0.0\n")

    # Scale factors -> scalefactors_json.json
    img_dims = None
    for ts in cl.index_block.get("SpatialImageTiles", []):
        img_dims = ts.get("Dims")
        break

    if bin_px is not None:
        spot_diam = bin_px
    elif mpp:
        spot_diam = 10.0 / mpp  # 10 um default cell diameter for cell-seg data
    else:
        spot_diam = 1.0
    # tissue_hires_scalef/tissue_lowres_scalef are intentionally BOTH 1.0, not
    # values derived from a real SpaceRanger scalefactors_json.json. This
    # package only ever extracts and stores ONE raster (the highest-zoom tile
    # in the pyramid below), and that raster is in the same pixel space as
    # pxl_row/pxl_col (i.e. genuinely scale 1.0). Setting tissue_lowres_scalef
    # to a real SpaceRanger value (e.g. ~0.19) while only the hires raster
    # exists would break default plotting that assumes lowres scaling --
    # every spot would be mis-registered into the image's top-left corner.
    # Do not "fix" these to look more authentic.
    sf = {
        "tissue_hires_scalef": 1.0,
        "tissue_lowres_scalef": 1.0,
        "fiducial_diameter_fullres": spot_diam,
        "spot_diameter_fullres": spot_diam,
        "microns_per_pixel": mpp,
        "bin_size_um": (mpp * bin_px) if (mpp and bin_px) else None,
    }
    if img_dims:
        sf["image_width_px"]  = img_dims[0]
        sf["image_height_px"] = img_dims[1]

    with open(os.path.join(outdir, "scalefactors_json.json"), "w") as f:
        json.dump(sf, f, indent=2)

    # Tissue image
    if include_image:
        print("[cloupe] Stitching image tiles ...", flush=True)
        img = stitch_tiles(cl)
        if img is not None:
            img.save(os.path.join(outdir, "tissue_hires_image.png"))
            w, h = img.size
            print(f"[cloupe] Image saved: {w}x{h} px", flush=True)
        else:
            print("[cloupe] No image tiles found or Pillow unavailable", flush=True)

    # Non-spatial projections (UMAP, tSNE, ...)
    other_projs = [
        p for p in cl.projections
        if p["name"].lower() not in ("spatial", "fiducials")
    ]
    if other_projs:
        with open(os.path.join(outdir, "projections.csv"), "w") as f:
            header = ["barcode"]
            for p in other_projs:
                for d in range(p["ndim"]):
                    header.append(f"{p['name']}_{d + 1}")
            f.write(",".join(header) + "\n")
            for i, bc in enumerate(barcodes):
                row = [bc]
                for p in other_projs:
                    for d in range(p["ndim"]):
                        row.append(f"{p['coordinates'][d][i]:.6f}")
                f.write(",".join(row) + "\n")

    # Spaceranger clusterings
    clusterings = read_clusterings(cl)
    if clusterings:
        cl_names = list(clusterings.keys())
        with open(os.path.join(outdir, "clusterings.csv"), "w") as f:
            f.write("barcode," + ",".join(cl_names) + "\n")
            for i, bc in enumerate(barcodes):
                row = [bc] + [clusterings[n][i] for n in cl_names]
                f.write(",".join(row) + "\n")

    # User cell tracks
    if cl.has_celltracks:
        ct_names = [t["Name"] for t in cl.celltracks]
        ct_values = [t["Values"] for t in cl.celltracks]
        with open(os.path.join(outdir, "celltracks.csv"), "w") as f:
            f.write("barcode," + ",".join(ct_names) + "\n")
            for i, bc in enumerate(barcodes):
                row = [bc] + [str(ct_values[j][i]) for j in range(len(ct_names))]
                f.write(",".join(row) + "\n")

    print(f"[cloupe] Done -> {outdir}", flush=True)
    return outdir


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python -m cloupe_extract.extract <cloupe_path> <outdir>")
        sys.exit(1)
    extract_cloupe(sys.argv[1], sys.argv[2])
