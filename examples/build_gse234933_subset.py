"""Rebuild the small tutorial fixture from GEO's original RDS count files.

Download GSE234933_MGH_HNSCC_gex_raw_counts.tar.gz and the cell annotation
table from GEO GSE234933. Extract the archive into SOURCE_DIR. This optional
preparation step needs rdata and is not part of the ten-minute tutorial.

python examples/build_gse234933_subset.py SOURCE_DIR --output NEW_OUTPUT

The checked-in manifest supplies expected input SHA-256 checksums. Extraction
fails if the deposited input snapshot has changed. No clinical attributes are
inferred. The source data are not covered by the software's MIT License.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy import sparse

SEED = 20261105
MAX_CELLS = 100
SPECIMENS = ["HN1", "HN13", "HN14", "HN17", "HN2", "HN20", "HN21", "HN22"]
SPP1_TAM_GENES = [
    "SPP1", "PCSK5", "SLC11A1", "VCAN", "SLC25A37", "FLNA", "UPP1", "BCL6",
    "AQP9", "TIMP1", "VEGFA", "ADM", "MARCO", "FN1", "IL1RN",
]
GENES = list(dict.fromkeys(["LILRB4", "MMP14", "MMP9", "SPP1", "APOE"] + SPP1_TAM_GENES))
DATA_DIR = Path(__file__).resolve().parents[1] / "src" / "pase" / "data"
GEO = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234933"
FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE234nnn/GSE234933/suppl/"
ANNOTATION = "GSE234933_MGH_HNSCC_cell_annotation.txt.gz"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def select_barcodes(barcodes: list[str], specimen: str) -> list[str]:
    """Hash selection depends only on seed and barcode, never gene expression."""
    def key(barcode: str) -> tuple[str, str]:
        digest = hashlib.sha256(f"{SEED}|{specimen}|{barcode}".encode("utf-8")).hexdigest()
        return digest, barcode
    return sorted(sorted(barcodes, key=key)[:MAX_CELLS])


def build(source: Path, output: Path, expected: dict[str, str]) -> dict:
    import rdata

    required = [ANNOTATION] + [f"{s}_gex_raw_counts.rds" for s in SPECIMENS]
    actual = {name: sha256(source / name) for name in required}
    for name in required:
        if actual[name] != expected[name]:
            raise ValueError(f"Source checksum mismatch: {name}")
    annotation = pd.read_csv(source / ANNOTATION, sep="\t", dtype=str)
    if not annotation.sample_barcode.is_unique:
        raise ValueError("Repeated sample barcode in source annotation")
    matrices, observations, source_dimensions = [], [], []
    for specimen in SPECIMENS:
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message='Missing constructor for R class "dgCMatrix".*')
            obj = rdata.read_rds(source / f"{specimen}_gex_raw_counts.rds")
        full = sparse.csc_matrix((obj.x, obj.i, obj.p), shape=tuple(obj.Dim))
        genes, barcodes = (np.asarray(values, dtype=str) for values in obj.Dimnames)
        if not pd.Index(genes).is_unique or not pd.Index(barcodes).is_unique:
            raise ValueError(f"Duplicated axis labels in {specimen}")
        if (not np.isfinite(full.data).all() or (full.data < 0).any()
                or not np.equal(full.data, np.floor(full.data)).all()):
            raise ValueError(f"Counts must be nonnegative finite integers: {specimen}")
        myeloid = annotation.loc[
            (annotation["sample"] == specimen) & (annotation.subclustering == "Myeloid cells"),
            "sample_barcode",
        ].tolist()
        if not set(myeloid) <= set(barcodes):
            raise ValueError(f"Unmatched annotated barcode: {specimen}")
        positions = pd.Index(barcodes).get_indexer(myeloid)
        libraries = np.asarray(full[:, positions].sum(axis=0)).ravel().astype(np.int64)
        positive = [bc for bc, total in zip(myeloid, libraries, strict=True) if total > 0]
        chosen = select_barcodes(positive, specimen)
        rows = pd.Index(barcodes).get_indexer(chosen)
        gene_positions = pd.Index(genes).get_indexer(GENES)
        if (gene_positions < 0).any():
            raise ValueError(f"Missing tutorial gene: {specimen}")
        # Crucial: full library sizes are computed BEFORE reducing the gene axis.
        selected_full = full[:, rows]
        full_library = np.asarray(selected_full.sum(axis=0)).ravel().astype(np.int64)
        selected_counts = selected_full[gene_positions, :].T.toarray().astype(np.int32)
        matrices.append(selected_counts)
        observations.append(pd.DataFrame({
            "barcode": chosen, "donor": specimen, "specimen": specimen,
            "myeloid": True, "full_library_size": full_library,
            "log_depth": np.log1p(full_library),
        }))
        source_dimensions.append({
            "specimen": specimen, "full_gene_count": int(full.shape[0]),
            "raw_matrix_barcode_count": int(full.shape[1]), "annotated_myeloid_cells": len(myeloid),
            "positive_library_myeloid_cells": len(positive), "demo_cells": len(chosen),
        })
        print(f"{specimen}: {len(chosen)} cells from {len(positive)} myeloid cells", flush=True)
    obs = pd.concat(observations, ignore_index=True)
    raw = np.concatenate(matrices)
    x = np.log1p(raw.astype(np.float64) * (10000.0 / obs.full_library_size.to_numpy())[:, None])
    output.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output / "gse234933_demo.npz", X=x, counts=raw, genes=np.asarray(GENES))
    obs.to_csv(output / "gse234933_demo_obs.tsv", sep="\t", index=False, float_format="%.17g", lineterminator="\n")
    manifest = {
        "schema_version": 1, "accession": "GSE234933", "source_url": GEO,
        "publication": {"doi": "10.1126/science.ade2292", "pmid": "37535729",
                        "citation": "Bill R, Wirapati P, Messemaker M, et al. CXCL9:SPP1 macrophage polarity identifies a network of cellular programs that control human cancers. Science. 2023;381:515-524."},
        "source_downloads": [FTP + "GSE234933_MGH_HNSCC_gex_raw_counts.tar.gz", FTP + ANNOTATION],
        "source_file_sha256": actual, "source_dimensions": source_dimensions,
        "selection": {
            "seed": SEED, "max_cells_per_specimen": MAX_CELLS,
            "specimens": SPECIMENS,
            "specimen_rule": "First eight lexicographically sorted specimen names, after excluding HN45; no effect-size selection.",
            "cell_rule": "Author annotation subclustering == Myeloid cells; positive full library; lowest SHA256(seed|specimen|barcode) values; output sorted by barcode.",
            "duplicate_donor_rule": "HN45 excluded because HN23 and HN45 are repeated specimens from one donor. Neither specimen is among the first eight selected specimens.",
        },
        "normalization": "log1p(10000 * gene_count / sum(counts across all original genes)); no feature-selection denominator",
        "log_depth": "log1p(sum(counts across all original genes))",
        "n_cells": len(obs), "n_donors": int(obs.donor.nunique()), "n_genes": len(GENES),
        "genes": GENES, "spp1_tam_genes": SPP1_TAM_GENES,
        "signature_source_doi": "10.3389/fimmu.2021.694801",
        "interpretation": "Small real-data software tutorial; no claim to reproduce full-cohort manuscript results or to provide independent biological validation.",
        "data_rights": "Third-party GEO-derived data; not relicensed under the package MIT License. Acknowledge GSE234933 and Bill et al.; see DATA_NOTICE.md and GEO data-use terms.",
        "geo_use_terms": "https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html",
        "outputs_sha256": {name: sha256(output / name) for name in ["gse234933_demo.npz", "gse234933_demo_obs.tsv"]},
        "builder_sha256": sha256(Path(__file__)),
    }
    (output / "gse234933_demo_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, help="Directory containing extracted original RDS files and annotation")
    parser.add_argument("--output", type=Path, required=True, help="New destination for rebuilt tutorial files")
    parser.add_argument("--manifest", type=Path, default=DATA_DIR / "gse234933_demo_manifest.json")
    args = parser.parse_args()
    lock = json.loads(args.manifest.read_text(encoding="utf-8"))
    build(args.source, args.output, lock["source_file_sha256"])
