import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from scipy import sparse, stats

import pase


def fixture_data(n=40, patients=5):
    rng = np.random.default_rng(2305)
    count = n * patients
    depth = rng.normal(size=count)
    x = rng.normal(size=count) + 0.3 * depth
    effect = np.repeat(np.linspace(-0.2, 0.6, patients), n)
    y = effect * x + rng.normal(size=count) + 0.2 * depth
    matrix = np.column_stack([x, y, y + rng.normal(size=count), np.ones(count)])
    obs = pd.DataFrame({"donor": np.repeat([f"p{i}" for i in range(patients)], n),
                        "log_depth": depth, "myeloid": True},
                       index=[f"cell{i}" for i in range(count)])
    return AnnData(matrix, obs=obs,
                   var=pd.DataFrame(index=["A", "B", "C", "constant"]))


def test_equal_matches_manual_partial_correlation_and_t():
    adata = fixture_data()
    out = pase.associate(adata, "A", "B", "donor", covariates=["log_depth"], min_patients=2)
    expected = []
    for pid in adata.obs.donor.unique():
        selected = adata.obs.donor == pid
        x, y = adata.X[selected, :2].T
        z = np.column_stack([np.ones(selected.sum()), adata.obs.loc[selected, "log_depth"]])
        x = x - z @ np.linalg.lstsq(z, x, rcond=None)[0]
        y = y - z @ np.linalg.lstsq(z, y, rcond=None)[0]
        expected.append(np.dot(x, y) / np.sqrt(np.dot(x, x) * np.dot(y, y)))
    expected = np.array(expected)
    summary = out.summary().iloc[0]
    se = expected.std(ddof=1) / np.sqrt(len(expected))
    np.testing.assert_allclose(out.patients.r, expected, atol=1e-13)
    assert summary.estimate == pytest.approx(expected.mean())
    assert summary.ci_low == pytest.approx(expected.mean() - stats.t.ppf(.975, 4) * se)
    assert summary.p_value == pytest.approx(2 * stats.t.sf(abs(expected.mean() / se), 4))
    assert summary.n_patients_eligible == 5
    assert summary.status == "reported"


def test_module_standardized_before_covariate_missingness():
    adata = fixture_data()
    adata.obs.loc[adata.obs.index[:10], "log_depth"] = np.nan
    out = pase.associate(adata, "A", ["B", "C", "constant"], "donor",
                         covariates=["log_depth"], min_cells=10)
    a = adata.X[:40]
    standardized = (a[:, 1:3] - a[:, 1:3].mean(0)) / a[:, 1:3].std(0, ddof=1)
    score = standardized.sum(1) / 3
    good = np.arange(10, 40)
    z = np.column_stack([np.ones(len(good)), adata.obs.log_depth.to_numpy()[good]])
    x = a[good, 0] - z @ np.linalg.lstsq(z, a[good, 0], rcond=None)[0]
    y = score[good] - z @ np.linalg.lstsq(z, score[good], rcond=None)[0]
    assert out.patients.iloc[0].r == pytest.approx(np.corrcoef(x, y)[0, 1])
    assert out.patients.iloc[0].n_expression_complete == 40
    assert out.patients.iloc[0].n_complete == 30


def test_equal_invariant_to_exact_donor_cell_duplication():
    adata = fixture_data()
    one = pase.associate(adata, "A", ["B", "C"], "donor", covariates=["log_depth"], min_patients=2)
    selected = np.concatenate([np.arange(adata.n_obs), np.tile(np.arange(40), 9)])
    with pytest.warns(UserWarning, match="Observation names are not unique"):
        duplicated = adata[selected].copy()
    duplicated.obs_names = [f"dup{i}" for i in range(duplicated.n_obs)]
    two = pase.associate(duplicated, "A", ["B", "C"], "donor", covariates=["log_depth"], min_patients=2)
    np.testing.assert_allclose(one.patients.r, two.patients.r, atol=1e-13)
    np.testing.assert_allclose(one.summary()[["estimate", "ci_low", "ci_high"]],
                               two.summary()[["estimate", "ci_low", "ci_high"]], atol=1e-13)


def test_sparse_dense_parity_and_no_mutation():
    adata = fixture_data()
    original = adata.copy()
    out = pase.associate(adata, "A", ["B", "C"], "donor")
    sparse_adata = adata.copy()
    sparse_adata.X = sparse.csr_matrix(adata.X)
    other = pase.associate(sparse_adata, "A", ["B", "C"], "donor")
    pd.testing.assert_frame_equal(out.patients, other.patients)
    np.testing.assert_array_equal(adata.X, original.X)
    pd.testing.assert_frame_equal(adata.obs, original.obs)
    assert sparse.issparse(sparse_adata.X)


def test_retains_masked_out_absent_and_below_threshold_patients():
    adata = fixture_data()
    adata.obs.loc[adata.obs.donor == "p1", "myeloid"] = False
    out = pase.associate(adata, "A", "B", "donor", cell_mask="myeloid",
                         patient_ids=["p0", "p1", "p2", "p3", "p4", "absent"],
                         min_patients=5)
    table = out.patients.set_index("patient_id")
    assert table.loc["p1", "reason"] == "no_selected_cells"
    assert table.loc["absent", "reason"] == "no_observed_cells"
    assert len(table) == 6
    assert out.summary().iloc[0].status == "withheld_policy"
    assert pd.isna(out.summary().iloc[0].estimate)
    assert pd.isna(out.summary().iloc[0].p_value)
    assert table.loc["p0", "eligible"]
    assert np.isfinite(table.loc["p0", "r"])
    small = pase.associate(adata, "A", "B", "donor", min_cells=50)
    assert not small.patients.eligible.any()
    assert small.patients.r.notna().all()


def test_re_hksj_matches_manual_two_patient_fit():
    adata = fixture_data(patients=2)
    out = pase.associate(adata, "A", "B", "donor", aggregate="re_hksj", min_patients=2)
    z = np.arctanh(out.patients.r.to_numpy())
    v = 1 / 37
    tau = max(0., ((z[0] - z[1]) ** 2 - 2 * v) / 2)
    weights = np.full(2, 1 / (v + tau))
    mu = np.average(z, weights=weights)
    q = np.sum(weights * (z - mu) ** 2)
    se = np.sqrt(q / weights.sum())
    summary = out.summary().iloc[0]
    assert summary.tau2 == pytest.approx(tau)
    assert summary.estimate == pytest.approx(np.tanh(mu))
    assert summary.ci_low == pytest.approx(np.tanh(mu - stats.t.ppf(.975, 1) * se))
    assert out.patients.weight.sum() == pytest.approx(1)


def test_fisher_boundary_kept_for_equal_only():
    adata = fixture_data()
    adata.X[:40, 1] = adata.X[:40, 0]
    out = pase.associate(adata, "A", "B", "donor", aggregate="re_hksj")
    row = out.patients.iloc[0]
    assert row.r == 1
    assert row.equal_eligible
    assert not row.fisher_eligible
    assert row.fisher_reason == "boundary_r"
    assert out.summary().iloc[0].n_patients_eligible == 4


def test_eb_has_no_group_inference_and_policy_suppresses_fit():
    adata = fixture_data()
    out = pase.associate(adata, "A", "B", "donor", aggregate="eb", min_patients=2)
    summary = out.summary().iloc[0]
    assert pd.isna(summary.estimate)
    assert pd.isna(summary.ci_low)
    assert pd.isna(summary.p_value)
    assert out.patients.eb_estimate.notna().all()
    assert "conditional" in out.metadata["eb_interval_interpretation"]
    withheld = pase.associate(adata, "A", "B", "donor", aggregate="eb", min_patients=6)
    assert withheld.patients.eb_estimate.isna().all()
    assert withheld.summary().iloc[0].status == "withheld_policy"
    assert pd.isna(withheld.summary().iloc[0].tau2)


def test_zero_heterogeneity_eb_has_no_false_intervals():
    adata = fixture_data()
    for start in range(40, adata.n_obs, 40):
        adata.X[start:start + 40] = adata.X[:40]
    out = pase.associate(adata, "A", "B", "donor", aggregate="eb", min_patients=2)
    assert out.patients.eb_estimate.notna().all()
    assert out.patients.eb_ci_low.isna().all()
    assert (out.patients.eb_status == "zero_tau2_complete_pooling").all()


def test_collinear_covariates_use_design_rank():
    adata = fixture_data()
    adata.obs["duplicate"] = adata.obs.log_depth * 2
    one = pase.associate(adata, "A", "B", "donor", covariates=["log_depth"])
    two = pase.associate(adata, "A", "B", "donor", covariates=["log_depth", "duplicate"])
    np.testing.assert_allclose(one.patients.r, two.patients.r)
    assert (two.patients.design_rank == 2).all()


def test_covariate_unit_changes_preserve_projection_and_design_rank():
    adata = fixture_data()
    original = pase.associate(adata, "A", ["B", "C"], "donor", covariates=["log_depth"])
    adata.obs["log_depth"] = 1e15 * adata.obs.log_depth
    scaled = pase.associate(adata, "A", ["B", "C"], "donor", covariates=["log_depth"])
    np.testing.assert_array_equal(original.patients.design_rank, scaled.patients.design_rank)
    np.testing.assert_allclose(original.patients.r, scaled.patients.r, atol=1e-13)


@pytest.mark.parametrize("mask", [["yes"] * 200, np.ones(200), [True] * 199,
                                   [True] * 199 + [None]])
def test_invalid_masks_fail(mask):
    with pytest.raises((ValueError, TypeError)):
        pase.associate(fixture_data(), "A", "B", "donor", cell_mask=mask)


def test_series_mask_index_must_match():
    adata = fixture_data()
    mask = pd.Series(True, index=adata.obs_names[::-1])
    with pytest.raises(ValueError, match="index"):
        pase.associate(adata, "A", "B", "donor", cell_mask=mask)


@pytest.mark.parametrize("kwargs", [{"target": ["B", "B"]}, {"target": "missing"},
                                    {"target": ["A", "B"]}, {"min_cells": 1.5},
                                    {"min_patients": 1}, {"alpha": True},
                                    {"patient_ids": ["p0", "p0"]},
                                    {"aggregate": "unknown"}])
def test_invalid_arguments_fail(kwargs):
    args = dict(anchor="A", target="B", patient_key="donor")
    args.update(kwargs)
    with pytest.raises((TypeError, ValueError, KeyError)):
        pase.associate(fixture_data(), **args)


def test_missing_patient_ids_and_nonnumeric_covariates_fail():
    adata = fixture_data()
    adata.obs.loc[adata.obs.index[0], "donor"] = None
    with pytest.raises(ValueError, match="patient"):
        pase.associate(adata, "A", "B", "donor")
    adata = fixture_data()
    adata.obs["text"] = "10"
    with pytest.raises(TypeError, match="numeric"):
        pase.associate(adata, "A", "B", "donor", covariates=["text"])


def test_nonfinite_expression_complete_cases_are_counted():
    adata = fixture_data()
    adata.X[0, 1] = np.nan
    adata.X[1, 0] = np.inf
    out = pase.associate(adata, "A", "B", "donor")
    assert out.patients.iloc[0].n_complete == 38
    assert out.patients.iloc[0].n_nonfinite_expression == 2
    assert np.isfinite(out.patients.iloc[0].r)


def test_layer_selection_and_constant_residual_exclusion():
    adata = fixture_data()
    adata.layers["transformed"] = adata.X.copy()
    original = pase.associate(adata, "A", "B", "donor")
    adata.X[:, 1] = 1
    constant = pase.associate(adata, "A", "B", "donor")
    assert (constant.patients.reason == "zero_residual_variance").all()
    from_layer = pase.associate(adata, "A", "B", "donor", layer="transformed")
    pd.testing.assert_frame_equal(original.patients, from_layer.patients)


def test_residual_degrees_of_freedom_gate_survives_low_cell_threshold():
    adata = fixture_data(n=4)
    out = pase.associate(adata, "A", "B", "donor", covariates=["log_depth"], min_cells=1)
    assert (out.patients.residual_df == 2).all()
    assert (out.patients.reason == "insufficient_residual_df").all()
    assert not out.patients.eligible.any()


def test_duplicate_gene_names_and_stringified_patient_collisions_fail():
    adata = fixture_data()
    adata.var_names = ["A", "B", "C", "B"]
    with pytest.raises(ValueError, match="duplicate"):
        pase.associate(adata, "A", "B", "donor")
    adata = fixture_data()
    ids = adata.obs.donor.to_numpy(dtype=object)
    ids[:40] = 1
    ids[40:80] = "1"
    adata.obs["donor"] = ids
    with pytest.raises(ValueError, match="collide"):
        pase.associate(adata, "A", "B", "donor")


def test_empty_anndata_keeps_declared_patients():
    adata = fixture_data()[:0].copy()
    out = pase.associate(adata, "A", "B", "donor", patient_ids=["not_sampled"])
    assert out.patients.iloc[0].reason == "no_observed_cells"
    assert out.summary().iloc[0].n_patients_total == 1


def test_plot_returns_axes():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out = pase.associate(fixture_data(), "A", "B", "donor")
    ax = out.plot_forest()
    assert isinstance(ax, plt.Axes)
    plt.close(ax.figure)
