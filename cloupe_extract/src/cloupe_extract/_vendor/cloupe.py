# ---------------------------------------------------------------------------
# Vendored from cellgeni/cloupe (https://github.com/cellgeni/cloupe),
# commit c0d0564545e163eea95ec24844ace03381786376 (2025-01-13), file
# `cloupe/__init__.py`, unmodified except for this header comment.
#
# Authors: Martin Prete <martin.prete@sanger.ac.uk> and
#          Nithin Mathew Joseph <nj9@sanger.ac.uk>
#          Wellcome Sanger Institute, Cellular Genetics Informatics (cellgeni)
#
# This is the actual reverse-engineering work that makes reading .cloupe
# files possible at all -- everything in cloupe_extract.extract is built on
# top of the Cloupe class below, unchanged. See ../../../THIRD_PARTY_LICENSES/
# for the original AGPL-3.0 license text, preserved as required.
#
# Vendored (rather than kept as an external pip dependency) specifically to
# pin an exact, tested version of this parser and eliminate the class of
# external-path-configuration bugs that affected earlier versions of this
# project. See cloupe_extract's README for the licensing implications of
# vendoring AGPL-3.0 code (this vendoring is why cloupe_extract, and by
# extension loupe2py and Loupe2R which depend on it, are AGPL-3.0-or-later).
# ---------------------------------------------------------------------------

# __init__.py
__version__ = "0.1.0"

import csv
import gzip
import io
import json
import logging
import os
import struct

# Set the logging level and message format
logging.basicConfig(level=logging.DEBUG, format="[%(asctime)s][%(levelname)s] %(message)s")


class Cloupe(object):
    """
    A class to parse, process, and export data from `.cloupe` files.

    Attributes:
        cloupe_path (str): Path to the `.cloupe` file.
        header (dict): Parsed header of the file.
        index_block (dict): Index block containing metadata and sections of the file.
        runs (list): List of runs metadata.
        metrics (list): List of metrics.
        analyses (list): List of analyses.
        projections (list): List of projections with coordinate data.
        matrices (list): List of matrices with properties (barcodes and feature names).
        celltracks (list): List of cell tracks (user annotations) with metadata and values.

    Methods:
        read_header: Reads the header section of the file.
        get_file_size: Returns the file size in bytes.
        decompress: Decompresses a zlib-compressed chunk.
        read_block: Reads a block from the file.
        to_csv: Exports components of the file to CSV format.
        to_anndata: Converts data into AnnData format.
    """

    def __init__(self, cloupe_path, load_csr=False):
        """
        Initializes the Cloupe object by reading and parsing the contents of `.cloupe` file.

        Args:
            cloupe_path (str): Path to the `.cloupe` file.
        """

        logging.debug(f"Cloupe file path: {cloupe_path}")
        self.cloupe_path = cloupe_path

        # read main header
        logging.debug(f"Reading header...")
        self.header = self.read_header()
        logging.debug(f"Version {self.header['version']}")

        # index block will tell us where the data sections start/end
        logging.debug(f"Reading indexBlock...")
        self.index_block = self.read_block(
            start=self.header["indexBlock"]["Start"],
            end=self.header["indexBlock"]["End"],
            as_json=True,
        )
        logging.debug(f"Format version {self.index_block['FormatVersion']}")

        # read Runs
        self.runs = []
        logging.debug(f"Reading Runs...")
        for run in self.index_block.get("Runs", []):
            r = run.copy()
            logging.debug(f"Format version {r['FormatVersion']}")
            if run.get("Metadata", False):
                r["Metadata"] = self.read_block(
                    start=run["Metadata"]["Start"],
                    end=run["Metadata"]["End"],
                    as_json=True,
                )
            else:
                logging.debug(f"Runs don't have Metadata")
            self.runs.append(r)

        # read Metrics
        self.metrics = []
        logging.debug(f"Reading Metrics...")
        for metric in self.index_block.get("Metrics", []):
            logging.debug(f"Format version {metric['FormatVersion']}")
            self.metrics.append(
                self.read_block(
                    start=metric["Contents"]["Start"],
                    end=metric["Contents"]["End"],
                    as_json=True,
                )
            )

        # read Analyses
        self.analyses = []
        logging.debug(f"Reading Analyses...")
        for analysis in self.index_block.get("Analyses", []):
            logging.debug(f"Format version {analysis['FormatVersion']}")
            self.analyses.append(analysis)

        # read Projections
        self.projections = []
        logging.debug(f"Reading Projections...")
        for projection in self.index_block.get("Projections", []):
            logging.debug(f"Projection name {projection['Name']}")
            logging.debug(f"Format version {projection['FormatVersion']}")
            try:
                pname = projection["Name"]
                cols, rows = tuple(projection["Dims"])
                logging.debug(f"Dims ({cols} x {rows})")

                # read data block
                pblock_bytes = self.read_block(
                    start=projection["Matrix"]["Start"],
                    end=projection["Matrix"]["End"],
                )
                pblock_count = projection["Matrix"]["ArraySize"]
                pblock = struct.unpack(f"{pblock_count}d", pblock_bytes)

                # read index block
                # pindex_bytes = self.read_block(
                #    start=projection["Matrix"]["Index"]["Start"],
                #    end=projection["Matrix"]["Index"]["End"],
                # )
                # pindex_count = projection["Matrix"]["Index"]["ArraySize"]
                # pindex = struct.unpack(f"{pindex_count}Q", pindex_bytes)

                # reashape so it has the dimensions of the matrix
                self.projections.append(
                    {
                        "name": pname,
                        "ndim": cols,
                        "coordinates": [pblock[x : x + rows] for x in range(0, len(pblock), rows)],
                    }
                )
            except Exception as e:
                logging.debug(f"Projection '{projection['Name']}' error: {e}")

        # read Matrices
        self.matrices = []
        logging.debug(f"Reading Matrices...")
        for matrix in self.index_block.get("Matrices", []):
            mtx = matrix.copy()

            for prop in [
                "Barcodes",
                "FeatureIds",
                "FeatureNames",
                "FeatureSecondaryNames",
            ]:
                if prop not in matrix:
                    logging.debug(f"Matrix doesn't have property '{prop}'")
                    continue
                mblock = self.read_block(start=matrix[prop]["Start"], end=matrix[prop]["End"])
                step = matrix[prop]["ArrayWidth"]
                stop = matrix[prop]["ArraySize"] * step
                mtx[prop] = [
                    mblock[i : i + step].decode("utf-8", errors="replace").strip("\x00")
                    for i in range(0, stop, step)
                ]

            # load UMI counts for this matricx
            umicounts_bytes = self.read_block(
                start=matrix["UMICounts"]["Start"],
                end=matrix["UMICounts"]["End"],
            )
            umicounts_count = matrix["UMICounts"]["ArraySize"]
            mtx["UMICounts"] = struct.unpack(f"{umicounts_count}Q", umicounts_bytes)

            # load the sparce matrix (CSR)
            if load_csr:
                import scipy.sparse

                # data
                csr_values_bytes = self.read_block(
                    start=matrix["CSRValues"]["Start"], end=matrix["CSRValues"]["End"]
                )
                csr_values = struct.unpack(f"{matrix['CSRValues']['ArraySize']}d", csr_values_bytes)
                # indptr
                csr_pointers_bytes = self.read_block(
                    start=matrix["CSRPointers"]["Start"], end=matrix["CSRPointers"]["End"]
                )
                csr_pointers = struct.unpack(
                    f"{matrix['CSRPointers']['ArraySize']}Q", csr_pointers_bytes
                )
                # indices
                csr_indices_bytes = self.read_block(
                    start=matrix["CSRIndices"]["Start"], end=matrix["CSRIndices"]["End"]
                )
                csr_indices = struct.unpack(
                    f"{matrix['CSRIndices']['ArraySize']}Q", csr_indices_bytes
                )

                # csr matrix is shaped Features x Barcodes
                mtx["CSR"] = scipy.sparse.csr_matrix(
                    (csr_values, csr_indices, csr_pointers),
                    shape=(matrix["FeatureCount"], matrix["BarcodeCount"]),
                )
            self.matrices.append(mtx)

        # read the next header (user content is stored here)
        self.next_index_block = {}
        self.has_celltracks = False
        self.celltracks = []
        if self.header["nextHeaderOffset"] != self.get_file_size():
            logging.debug(f"Reading next header...")
            self.next_header = self.read_header(self.header["nextHeaderOffset"])
            self.next_index_block = self.read_block(
                start=self.next_header["indexBlock"]["Start"],
                end=self.next_header["indexBlock"]["End"],
                as_json=True,
            )
            # now that we have the index block we get the 'userCreated' celltracks
            logging.debug(f"Reading CellTracks")
            for celltrack in self.next_index_block["CellTracks"]:
                self.has_celltracks = True
                ct = celltrack.copy()
                logging.debug(f'Celltrack {ct["Name"]}')
                ct["Metadata"] = self.read_block(
                    start=celltrack["Metadata"]["Start"],
                    end=celltrack["Metadata"]["End"],
                    as_json=True,
                )
                vblock_bytes = self.read_block(
                    start=celltrack["Values"]["Start"],
                    end=celltrack["Values"]["End"],
                )
                vblock_count = celltrack["Values"]["ArraySize"]
                try:
                    group_keys = struct.unpack(f"{vblock_count}h", vblock_bytes)
                    ct["Values"] = [
                        ct["Metadata"]["groups"][gk] if gk != -1 else "UNASSIGNED"
                        for gk in group_keys
                    ]
                except Exception as e:
                    logging.warning(f"Fail to get celltrack groups/values: {e}")
                    ct["Values"] = []
                self.celltracks.append(ct)

    def read_header(self, position=0):
        """
        Reads the header section of the `.cloupe` file.

        Args:
            position (int): Byte offset to start reading the header.

        Returns:
            dict: Parsed JSON object from the header.
        """
        with open(self.cloupe_path, "rb") as f:
            f.seek(position)
            header = f.read(4096)
            return json.loads(header.decode("utf-8", errors="replace").strip("\x00"))

    def get_file_size(self):
        """
        Returns the size of the `.cloupe` file in bytes.

        Returns:
            int: File size in bytes.
        """
        with open(self.cloupe_path, "rb") as f:
            f.seek(0, 2)
            size = f.tell()
            return size

    def __repr__(self):
        return str(self.header) + str(self.index_block["Runs"])

    def decompress(self, chunk):
        """
        Decompresses a gzip-compressed chunk.

        Args:
            chunk (bytes): Compressed data.

        Returns:
            bytes: Decompressed data.
        """
        with gzip.GzipFile(fileobj=io.BytesIO(chunk)) as gz:
            return gz.read()

    def read_block(self, start, end, as_json=False, verbose=False):
        """
        Reads a block of data from the file.

        Args:
            start (int): Start byte offset.
            end (int): End byte offset.
            as_json (bool): Whether to decode the block as JSON.

        Returns:
            bytes or dict: Raw bytes or decoded JSON object.
        """
        with open(self.cloupe_path, "rb") as f:
            if verbose:
                logging.debug(f"reading block [{start} : {end}]")
            f.seek(start)
            read_size = end - start
            block = f.read(read_size)
            if block[:2] == b"\x1f\x8b":
                if verbose:
                    logging.debug("compress block")
                block = self.decompress(block)
            if as_json:
                return json.loads(block.decode("utf-8", errors="replace"))
            return block

    def debug_blocks(self):
        return {
            "index": self.index_block,
            "runs": bool(self.runs),
            "metrics": bool(self.metrics),
            "analyses": bool(self.analyses),
            "projections": bool(self.projections),
            "matrices": bool(self.matrices),
            "has_celltracks": self.has_celltracks,
            "celltracks": bool(self.celltracks),
        }

    def to_csv(self, outdir=".", barcodes=True, features=True, celltracks=True, projections=True):
        """
        Exports various components of the `.cloupe` file to CSV format.

        Args:
            outdir (str): Directory to save the CSV files.
            barcodes (bool): Whether to export barcodes to a file named `barcodes.csv`.
            features (bool): Whether to export features to a file named `features.csv`.
            celltracks (bool): Whether to export cell tracks to a file named `celltracks.csv`.
            projections (bool): Whether to export projections to a file named `projections.csv`.
        """
        if barcodes:
            with open(os.path.join(outdir, "barcodes.csv"), "wt") as f:
                header = ["Barcode"]
                writer = csv.writer(f)
                writer.writerow(header)
                for barcode in self.matrices[0]["Barcodes"]:
                    writer.writerow([barcode])

        if features:
            with open(os.path.join(outdir, "features.csv"), "wt") as f:
                header = ["FeatureId", "FeatureName"]
                writer = csv.writer(f)
                writer.writerow(header)
                features = list(
                    zip(self.matrices[0]["FeatureIds"], self.matrices[0]["FeatureNames"])
                )
                writer.writerows(features)

        if celltracks:
            with open(os.path.join(outdir, "celltracks.csv"), "wt") as f:
                header = ["Barcode"] + [track["Name"] for track in self.celltracks]
                writer = csv.writer(f)
                writer.writerow(header)
                tracks = zip(*[track["Values"] for track in self.celltracks])
                writer.writerows(list(zip(self.matrices[0]["Barcodes"], *tracks)))

        if projections:
            with open(f"projections.csv", "wt") as f:
                header = ["Barcode"]
                # "Barcode" is the first column the rest of the columns are trackA_0,trackA_1,trackB_0,trackB_1,etc.
                for p in self.projections:
                    for i in range(p["ndim"]):
                        header.append(f'{p["name"]}_{i}')
                writer = csv.writer(f)
                writer.writerow(header)
                rows = []
                # first melt all dimensions to a common array (rows)
                # then zip barcodes and rows (all columns) to write rows
                for p in self.projections:
                    for i in range(p["ndim"]):
                        rows.append(p["coordinates"][i])
                writer.writerows(list(zip(self.matrices[0]["Barcodes"], *rows)))

    def to_anndata(self, filename=""):
        """
        Converts the `.cloupe` data into an AnnData object and saves it as an H5AD file.

        Args:
            filename (str): Name of the output H5AD file. Defaults to replacing the `.cloupe` file extension with `.h5ad`.

        Raises:
            Exception: If the `.cloupe` file wasn't parsed with `load_csr=True`.
        """
        import anndata as ad
        import numpy as np

        matrix = self.matrices[0]
        if "CSR" not in matrix:
            raise Exception(
                "Cloupe file wasn't parsed with `load_csr=True`. Reload the file with that option enabled to create an Anndata file."
            )

        # feature x barcode so it needs to be transposed for Annadata
        adata = ad.AnnData(X=matrix["CSR"]).T
        adata.obs_names = matrix["Barcodes"]
        adata.var_names = matrix["FeatureIds"]
        adata.var["FeatureNames"] = matrix["FeatureNames"]

        for projection in self.projections:
            # validate that the projection matches obs_names length
            if len(projection["coordinates"][0]) == len(adata.obs_names):
                adata.obsm[projection["name"]] = np.array(projection["coordinates"]).T

        for celltrack in self.celltracks:
            # validate that the celltrack matches obs_names length
            if len(celltrack["Values"]) == len(adata.obs_names):
                adata.obs[celltrack["Name"]] = celltrack["Values"]

        if not filename:
            filename = os.path.basename(self.cloupe_path).replace(".cloupe", ".h5ad")

        adata.write(filename=filename)
