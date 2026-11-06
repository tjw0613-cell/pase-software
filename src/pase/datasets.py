"""Small real-data examples supplied with PASE."""

from __future__ import annotations

from copy import deepcopy
from importlib.resources import files
from io import BytesIO, StringIO
import hashlib
import json

from anndata import AnnData
import numpy as np
import pandas as pd


def load_gse234933_demo() -> AnnData:
    """Load a small real GSE234933 myeloid-cell subset without network access.

    Returns
    -------
    AnnData
        ``X`` is log1p(CP10K) based on each cell's FULL original gene library.
        ``layers['counts']`` contains untransformed counts for the included genes;
        these genes alone MUST NOT be used to recalculate library sizes.
        ``obs`` retains original barcodes, donor/specimen labels, ``myeloid``,
        ``full_library_size`` and ``log_depth``. ``uns['dataset']`` supplies the
        source, extraction rules, source checksums and third-party data notice.

    Notes
    -----
    This teaching subset does not reproduce the full-cohort manuscript analysis.
    It is selected without consulting gene correlations. Each call returns an
    independent in-memory object. Files are checked against their SHA-256 manifest.
    """
    root = files("pase").joinpath("data")
    manifest = json.loads(root.joinpath("gse234933_demo_manifest.json").read_text(encoding="utf-8"))
    content = {}
    for name, expected in manifest["outputs_sha256"].items():
        value = root.joinpath(name).read_bytes()
        if hashlib.sha256(value).hexdigest() != expected:
            raise RuntimeError(f"Packaged tutorial data checksum mismatch: {name}")
        content[name] = value
    with np.load(BytesIO(content["gse234933_demo.npz"]), allow_pickle=False) as stored:
        x, counts, genes = stored["X"].copy(), stored["counts"].copy(), stored["genes"].astype(str)
    obs = pd.read_csv(StringIO(content["gse234933_demo_obs.tsv"].decode("utf-8")), sep="\t", index_col="barcode")
    obs["myeloid"] = obs["myeloid"].astype(bool)
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    var["spp1_tam_signature"] = var.index.isin(manifest["spp1_tam_genes"])
    return AnnData(X=x, obs=obs, var=var, layers={"counts": counts}, uns={"dataset": deepcopy(manifest)})
