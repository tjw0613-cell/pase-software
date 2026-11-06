"""Checks that the public tutorial preserves genuine counts and normalization."""

import numpy as np

from pase.datasets import load_gse234933_demo


def test_demo_keeps_full_library_denominator():
    adata = load_gse234933_demo()
    total = adata.obs.full_library_size.to_numpy()
    raw = adata.layers["counts"]
    assert np.issubdtype(raw.dtype, np.integer)
    assert np.all(total > raw.sum(axis=1))
    np.testing.assert_allclose(adata.X, np.log1p(raw * (10000.0 / total)[:, None]), rtol=1e-14)
    np.testing.assert_allclose(adata.obs.log_depth.to_numpy(), np.log1p(total), rtol=1e-14)


def test_demo_identity_selection_and_sources():
    adata = load_gse234933_demo()
    assert adata.shape == (776, 19)
    assert adata.obs.donor.nunique() == 8
    assert adata.obs.index.is_unique
    assert adata.obs.myeloid.all()
    assert "HN45" not in set(adata.obs.specimen)
    assert adata.obs.groupby("donor").size().max() == 100
    assert adata.obs.groupby("donor").size().min() == 76
    assert adata.var.spp1_tam_signature.sum() == 15
    for specimen, barcode in zip(adata.obs.specimen, adata.obs_names, strict=True):
        assert barcode.startswith(specimen + "_")
    info = adata.uns["dataset"]
    assert info["accession"] == "GSE234933"
    assert len(info["source_file_sha256"]) == 9
    assert all(row["full_gene_count"] > adata.n_vars for row in info["source_dimensions"])
    assert "not relicensed" in info["data_rights"]


def test_demo_returns_independent_objects():
    first = load_gse234933_demo()
    second = load_gse234933_demo()
    first.X[0, 0] = -1
    first.obs.iloc[0, first.obs.columns.get_loc("donor")] = "altered"
    first.uns["dataset"]["selection"]["seed"] = -1
    assert second.X[0, 0] >= 0
    assert second.obs.iloc[0].donor != "altered"
    assert second.uns["dataset"]["selection"]["seed"] == 20261105
