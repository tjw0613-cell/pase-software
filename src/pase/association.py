"""AnnData interface for patient-level molecular associations."""
from __future__ import annotations

from collections.abc import Sequence
from numbers import Integral, Real

import numpy as np
import pandas as pd
from scipy import sparse

from ._reporting import DEFAULT_MIN_PATIENTS, add_reporting, validate_planning
from ._statistics import summarize_patients
from .result import AssociationResult


def _integer(value, name, minimum):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")
    return int(value)


def _names(value, name, allow_empty=False):
    if isinstance(value, str):
        result = [value]
    elif isinstance(value, Sequence) or isinstance(value, (np.ndarray, pd.Index)):
        result = list(value)
    else:
        raise TypeError(f"{name} must be a column/gene name or a sequence of names.")
    if (not result and not allow_empty) or any(not isinstance(x, str) or not x.strip() for x in result):
        raise ValueError(f"{name} must contain nonempty string names.")
    if len(set(result)) != len(result):
        raise ValueError(f"{name} contains duplicate names.")
    return result


def _patient_names(values, name):
    result = []
    for value in values:
        if not np.isscalar(value) or pd.isna(value) or not str(value).strip():
            raise ValueError(f"{name} contains missing or empty patient identifiers.")
        result.append(str(value))
    return result


def _mask(obs, value):
    if value is None:
        return np.ones(len(obs), dtype=bool)
    if isinstance(value, str):
        if value not in obs:
            raise KeyError(f"cell_mask column {value!r} is absent from adata.obs.")
        value = obs[value]
    if isinstance(value, pd.Series) and not value.index.equals(obs.index):
        raise ValueError("A Series cell_mask must have exactly the same index and order as adata.obs.")
    arr = np.asarray(value)
    if arr.ndim != 1 or len(arr) != len(obs):
        raise ValueError("cell_mask must be a one-dimensional boolean vector with one value per cell.")
    if pd.isna(arr).any():
        raise ValueError("cell_mask must not contain missing values.")
    if not all(isinstance(v, (bool, np.bool_)) for v in arr):
        raise TypeError("cell_mask must contain booleans, not strings or numeric indicators.")
    return arr.astype(bool, copy=True)


def _partial_correlation(x, y, design):
    """Residualize against a rank-aware least-squares nuisance model."""
    rank = int(np.linalg.matrix_rank(design))
    residual_df = len(x) - rank
    if residual_df < 3:
        return np.nan, rank, residual_df, "insufficient_residual_df"
    # Center and scale outcomes first for numerical stability. The design has
    # an intercept, so these operations leave its partial correlation unchanged.
    x = np.asarray(x, dtype=float) - np.mean(x)
    y = np.asarray(y, dtype=float) - np.mean(y)
    x_scale = np.max(np.abs(x))
    y_scale = np.max(np.abs(y))
    if x_scale == 0 or y_scale == 0:
        return np.nan, rank, residual_df, "zero_residual_variance"
    x, y = x / x_scale, y / y_scale
    residuals = np.column_stack([x, y])
    residuals -= design @ np.linalg.lstsq(design, residuals, rcond=None)[0]
    ux, uy = residuals.T
    nx, ny = np.linalg.norm(ux), np.linalg.norm(uy)
    tol = 64 * np.finfo(float).eps
    if nx <= tol * np.linalg.norm(x) or ny <= tol * np.linalg.norm(y):
        return np.nan, rank, residual_df, "zero_residual_variance"
    correlation = float(np.dot(ux / nx, uy / ny))
    if not np.isfinite(correlation):
        return np.nan, rank, residual_df, "numerical_failure"
    correlation = float(np.clip(correlation, -1.0, 1.0))
    if 1 - abs(correlation) <= 8 * np.finfo(float).eps:
        correlation = float(np.sign(correlation))
    return correlation, rank, residual_df, None


def _patient_effect(pid, indices, mask, expression, covariate_values, module, min_cells):
    selected = indices[mask[indices]]
    result = dict(patient_id=pid, n_total=len(indices), n_selected=len(selected),
                  n_expression_complete=0, n_complete=0, n_cells=0,
                  n_nonfinite_expression=0, n_nonfinite_covariates=0,
                  design_rank=0, residual_df=0, r=np.nan, eligible=False,
                  failure_reason=None)
    if not len(indices):
        result["failure_reason"] = "no_observed_cells"
        return result
    if not len(selected):
        result["failure_reason"] = "no_selected_cells"
        return result
    values = expression[selected]
    finite_expression = np.isfinite(values).all(axis=1)
    result["n_nonfinite_expression"] = int((~finite_expression).sum())
    selected = selected[finite_expression]
    values = values[finite_expression]
    result["n_expression_complete"] = len(selected)
    if not len(selected):
        result["failure_reason"] = "no_finite_expression"
        return result
    x = values[:, 0]
    if module:
        outcomes = values[:, 1:]
        centered = outcomes - outcomes.mean(axis=0)
        scale = outcomes.std(axis=0, ddof=1) if len(outcomes) > 1 else np.zeros(outcomes.shape[1])
        varying = np.ptp(outcomes, axis=0) != 0
        standardized = np.divide(centered, scale, out=np.zeros_like(centered),
                                 where=varying & (scale > 0))
        y = standardized.mean(axis=1)
    else:
        y = values[:, 1]
    nuisance = covariate_values[selected]
    complete = np.isfinite(nuisance).all(axis=1)
    result["n_nonfinite_covariates"] = int((~complete).sum())
    n = int(complete.sum())
    result["n_complete"] = result["n_cells"] = n
    if not n:
        result["failure_reason"] = "no_complete_cases"
        return result
    # Standardization preserves the projection with an intercept while keeping
    # numerical rank invariant to covariate units. Constant columns add no rank.
    design_columns = [np.ones(n)]
    for values in nuisance[complete].T:
        if n > 1 and np.ptp(values) != 0:
            centered = values - values.mean()
            size = np.max(np.abs(centered))
            scaled = centered / size
            design_columns.append(scaled / scaled.std(ddof=1))
    design = np.column_stack(design_columns)
    try:
        r, rank, df, failure = _partial_correlation(x[complete], y[complete], design)
    except (np.linalg.LinAlgError, FloatingPointError, ValueError):
        r, rank, df, failure = np.nan, 0, 0, "numerical_failure"
    result.update(r=r, design_rank=rank, residual_df=df)
    result["failure_reason"] = failure or ("below_min_cells" if n < min_cells else None)
    result["eligible"] = result["failure_reason"] is None
    return result


def associate(adata, anchor, target, patient_key, covariates=None, cell_mask=None,
              aggregate="equal", min_cells=30, min_patients=DEFAULT_MIN_PATIENTS, layer=None,
              alpha=0.05, patient_ids=None, *, max_ci_width=0.2,
              minimum_relevant_effect=None, planning_sd=None, planning_power=0.8,
              planning_sd_source=None):
    """Estimate an association within each patient and summarize across patients.

    Parameters
    ----------
    adata : anndata.AnnData
        Cells by genes. X (or ``layer``) must already contain transformed
        expression, such as log1p(CP10K). This function does not normalize counts.
    anchor : str
        Anchor gene name in ``adata.var_names``.
    target : str or sequence of str
        A gene, or a module computed as the mean of within-patient gene z scores
        (sample SD, ddof=1). Constant measured genes contribute zero. Module
        standardization precedes covariate complete-case filtering.
    patient_key : str
        Column in ``adata.obs`` identifying independent patients.
    covariates : sequence of str, optional
        Numeric observation columns centered and scaled within each patient's
        complete cells, then included with an intercept in least-squares
        residualization. Constant columns add no rank; collinear covariates use
        the design rank. Scaling does not change the regression projection.
    cell_mask : str or boolean vector, optional
        Boolean observation column or cell-aligned selection. Patients with no
        selected cells remain in the output. Missing mask values are an error.
    aggregate : {'equal', 're_hksj', 'eb'}, default 'equal'
        Arithmetic mean of patient correlations with a t interval; Fisher-z
        random effects fitted by REML with classical HKSJ; or conditional
        empirical Bayes patient estimates, respectively. Their targets differ.
    min_cells : int, default 30
        Minimum complete cells for patient eligibility. Regardless of this
        setting, at least three residual degrees of freedom are required.
    min_patients : int, default 6
        Minimum method-specific eligible patients before inference is returned.
        Below this floor, only raw effects and descriptive summaries are returned.
        Values from two to five are expert overrides and generate a warning.
        Six is a reporting policy, not a validated adequacy threshold.
    layer : str, optional
        Transformed-expression layer instead of X. ``adata.raw`` is not used.
    alpha : float, default 0.05
        Two-sided interval error probability.
    patient_ids : sequence, optional
        Additional patient IDs to retain even if they have no observed cells.
        Observed IDs are always retained. Output IDs are strings.
    max_ci_width : float or None, keyword-only, default 0.2
        Warn when an available cohort interval exceeds this correlation-scale
        width. Must be positive and finite; None disables this warning. The
        interval is never changed or suppressed by this preference.
    minimum_relevant_effect : float, optional, keyword-only
        Relevant absolute mean correlation, in (0, 1]. Requires planning_sd.
        Warn when it is smaller than the model-based planning MDE.
    planning_sd : float, optional, keyword-only
        Externally assumed SD of patient correlations, in (0, 1]. Enables a
        two-sided noncentral-t MDE calculation for equal aggregation only.
        It is never estimated from the input cohort. Below min_patients the
        calculation is withheld; re_hksj and eb return not_supported.
    planning_power : float, keyword-only, default 0.8
        Desired planning power, strictly between zero and one. When planning_sd
        is provided, it must also exceed alpha.
    planning_sd_source : str, optional, keyword-only
        User-declared source of planning_sd, such as an independent pilot cohort.
        The software records but cannot verify it. Providing an SD without a
        source generates a warning. A source requires planning_sd.

    Returns
    -------
    AssociationResult
        Per-patient effects and exclusion reasons, a summary table, and metadata.

    Notes
    -----
    Cells with nonfinite anchor or target expression are excluded before module
    construction. Nonfinite covariates are then excluded. There is no outcome
    detection gate. Fisher variance is approximated as 1/(n-k-3), with k equal
    to design rank minus one; its independent-cell assumptions may not hold for
    all single-cell data. Raw patient correlations remain visible when cohort
    inference is withheld. Inputs are not modified.

    Reported inference below 20 eligible patients carries a risk warning;
    reaching 20 does not guarantee calibration. MDE is a planning quantity under
    independent, identically distributed normal patient effects with the supplied
    SD. It is not a retrospective power estimate or a count-data power guarantee.
    See result.metadata for the scope of the simulation-informed reporting policy.
    """
    min_cells = _integer(min_cells, "min_cells", 1)
    min_patients = _integer(min_patients, "min_patients", 2)
    if isinstance(alpha, (bool, np.bool_)) or not isinstance(alpha, Real) or not 0 < alpha < 1:
        raise ValueError("alpha must be a finite number strictly between zero and one.")
    (max_ci_width, minimum_relevant_effect, planning_sd,
     planning_power, planning_sd_source) = validate_planning(
        max_ci_width, minimum_relevant_effect, planning_sd, planning_power, float(alpha),
        planning_sd_source)
    if aggregate not in ("equal", "re_hksj", "eb"):
        raise ValueError("aggregate must be 'equal', 're_hksj', or 'eb'.")
    if not all(hasattr(adata, attr) for attr in ("obs", "var_names", "X", "layers")):
        raise TypeError("adata must provide the AnnData cells-by-genes interface.")
    obs = adata.obs
    if not obs.columns.is_unique:
        raise ValueError("adata.obs has duplicate column names.")
    if not isinstance(patient_key, str) or patient_key not in obs:
        raise KeyError("patient_key must name a column in adata.obs.")
    anchor = _names(anchor, "anchor")
    if len(anchor) != 1:
        raise ValueError("anchor must be a single gene name.")
    anchor = anchor[0]
    module = not isinstance(target, str)
    targets = _names(target, "target")
    if anchor in targets:
        raise ValueError("anchor must not also occur in target.")
    if not adata.var_names.is_unique:
        raise ValueError("adata.var_names contains duplicate gene names.")
    genes = [anchor, *targets]
    gene_indices = adata.var_names.get_indexer(genes)
    missing = [gene for gene, index in zip(genes, gene_indices) if index < 0]
    if missing:
        raise KeyError(f"Genes missing from adata.var_names: {missing}.")
    if layer is not None and (not isinstance(layer, str) or layer not in adata.layers):
        raise KeyError(f"Requested layer {layer!r} is absent from adata.layers.")
    matrix = adata.X if layer is None else adata.layers[layer]
    selected_matrix = matrix[:, gene_indices]
    if sparse.issparse(selected_matrix):
        selected_matrix = selected_matrix.toarray()
    expression = np.asarray(selected_matrix)
    if not np.issubdtype(expression.dtype, np.number) or np.iscomplexobj(expression):
        raise TypeError("Selected expression values must be real numeric data.")
    expression = expression.astype(float, copy=False)
    selected_mask = _mask(obs, cell_mask)
    covariates = [] if covariates is None else _names(covariates, "covariates", allow_empty=True)
    columns = []
    for name in covariates:
        if name not in obs:
            raise KeyError(f"Covariate {name!r} is absent from adata.obs.")
        if not pd.api.types.is_numeric_dtype(obs[name].dtype) or pd.api.types.is_complex_dtype(obs[name].dtype):
            raise TypeError(f"Covariate {name!r} must be numeric; encode categories explicitly.")
        columns.append(obs[name].to_numpy(dtype=float, na_value=np.nan))
    nuisance = np.column_stack(columns) if columns else np.empty((len(obs), 0))
    raw_ids = obs[patient_key].tolist()
    ids = _patient_names(raw_ids, "patient_key")
    # Do not merge distinct identifiers merely because their string forms match.
    original_unique = pd.unique(obs[patient_key])
    unique_ids = _patient_names(original_unique, "patient_key")
    if len(set(unique_ids)) != len(unique_ids):
        raise ValueError("Distinct patient identifiers collide after conversion to strings.")
    if patient_ids is not None:
        if isinstance(patient_ids, (str, bytes)):
            raise TypeError("patient_ids must be a sequence, not a single string.")
        additional_ids = _patient_names(list(patient_ids), "patient_ids")
        if len(set(additional_ids)) != len(additional_ids):
            raise ValueError("patient_ids contains duplicates.")
        unique_ids.extend(pid for pid in additional_ids if pid not in set(unique_ids))
    groups = {pid: [] for pid in unique_ids}
    for index, pid in enumerate(ids):
        groups[pid].append(index)
    records = [_patient_effect(pid, np.asarray(indices, dtype=int), selected_mask,
                               expression, nuisance, module, min_cells)
               for pid, indices in groups.items()]
    # The statistical kernel retains the verified rank-aware Fisher variance
    # and equal-r / REML / HK calculations used by the standalone implementation.
    fits = summarize_patients(records, alpha=float(alpha), shrinkage=aggregate == "eb")
    patient_columns = ["patient_id", "n_total", "n_selected", "n_expression_complete",
                       "n_complete", "n_cells", "n_nonfinite_expression", "n_nonfinite_covariates",
                       "design_rank", "residual_df", "r", "eligible", "failure_reason"]
    patients = pd.DataFrame(records, columns=patient_columns)
    flags = {row["patient_id"]: row for row in fits["patients"]}
    for key in ("equal_eligible", "equal_reason", "fisher_eligible", "fisher_reason", "z", "v"):
        patients[key] = [flags[pid][key] for pid in patients.patient_id]
    fisher = aggregate != "equal"
    patients["eligible"] = patients["fisher_eligible" if fisher else "equal_eligible"].astype(bool)
    patients["reason"] = patients["fisher_reason" if fisher else "equal_reason"]
    patients["status"] = np.where(patients.eligible, "eligible", "excluded")
    eligible = int(patients.eligible.sum())
    target_definition = ("Mean of eligible patients' observed within-patient correlations"
                         if aggregate == "equal" else
                         "Back-transformed mean of the random-effects distribution on the Fisher-z scale"
                         if aggregate == "re_hksj" else
                         "Individual patient associations shrunk toward the fitted Fisher-z population mean")
    summary = dict(method=aggregate, target=target_definition, n_patients_total=len(patients),
                   n_patients_eligible=eligible, n_patients_excluded=len(patients)-eligible,
                   estimate=np.nan, ci_low=np.nan, ci_high=np.nan, p_value=np.nan,
                   standard_error=np.nan, se_scale="r" if aggregate == "equal" else "Fisher z",
                   df=np.nan, tau2=np.nan, confidence_level=1-float(alpha),
                   min_cells=min_cells, min_patients=min_patients, status="withheld_technical", reason="")
    patients["weight"] = np.nan
    for key in ("eb_estimate", "eb_estimate_z", "eb_ci_low", "eb_ci_high", "eb_data_weight", "eb_conditional_se_z"):
        patients[key] = np.nan
    patients["eb_status"] = "not_requested" if aggregate != "eb" else "not_eligible"
    if aggregate == "equal":
        fitted = fits["equal_r"]
        interval = fitted
        if eligible:
            patients.loc[patients.eligible, "weight"] = 1 / eligible
    else:
        fitted = fits["reml"]
        interval = fitted["hk"]
    if eligible < min_patients:
        summary.update(status="withheld_policy", reason="below_min_eligible_patients")
        if aggregate == "eb":
            patients.loc[patients.eligible, "eb_status"] = "withheld_policy"
    elif aggregate == "eb":
        if fits["shrinkage_status"] in ("conditional_plugin", "zero_tau2_complete_pooling"):
            summary.update(status="individual_estimates_only", reason=fits["shrinkage_status"],
                           tau2=fitted["tau2"])
            eb_rows = {p["patient_id"]: p for p in fits["shrinkage"]}
            mapping = {"eb_estimate": "estimate", "eb_estimate_z": "estimate_z",
                       "eb_ci_low": "ci_low", "eb_ci_high": "ci_high",
                       "eb_data_weight": "observed_data_weight", "eb_conditional_se_z": "conditional_se_z",
                       "eb_status": "status"}
            for index, row in patients.loc[patients.eligible].iterrows():
                for public, private in mapping.items():
                    patients.loc[index, public] = eb_rows[row.patient_id][private]
        else:
            summary["reason"] = fits["shrinkage_status"]
    elif interval["status"] == "ok":
        summary.update(status="reported", reason="", estimate=fitted["estimate"],
                       ci_low=interval["ci_low"], ci_high=interval["ci_high"],
                       p_value=interval["p_value"], standard_error=interval["se"], df=interval["df"],
                       tau2=fitted.get("tau2", np.nan))
    else:
        summary["reason"] = interval["status"]
    if fisher and eligible >= min_patients and fitted["status"] == "ok":
        patients["weight"] = patients.patient_id.map(fitted["weights"])
    warnings = list(fits["warnings"])
    warnings.append("min_patients is a configurable reporting policy, not a validated adequacy threshold.")
    if fisher:
        warnings.append("Fisher sampling variances are approximations and assume independent cell-level information.")
    if aggregate == "eb":
        warnings.append("EB intervals condition on fitted population parameters and omit their estimation uncertainty.")
    if patients.n_nonfinite_expression.sum() or patients.n_nonfinite_covariates.sum():
        warnings.append("Nonfinite expression or covariates excluded cells; see per-patient complete-case counts.")
    metadata = dict(anchor=anchor, target_genes=targets,
                    target_label=targets[0] if not module else "mean z score (" + ", ".join(targets) + ")",
                    target_definition=target_definition,
                    target_score="within-patient mean z score; ddof=1; constant genes contribute zero" if module else "single transformed gene",
                    expression="Pretransformed input; no normalization performed",
                    expression_source="X" if layer is None else f"layers[{layer!r}]",
                    score_standardization="Selected cells with finite anchor and all target genes; before covariate complete cases",
                    patient_key=patient_key, covariates=covariates, cell_mask=cell_mask if isinstance(cell_mask, str) else "provided vector" if cell_mask is not None else None,
                    covariate_standardization="Centered and scaled within patient complete cases; constant columns omitted",
                    aggregate=aggregate, min_cells=min_cells, min_patients=min_patients,
                    minimum_residual_df=3, alpha=float(alpha), warnings=warnings,
                    eb_interval_interpretation="conditional plug-in intervals, holding fitted population parameters fixed; no cohort inference",
                    interval_method="Student t" if aggregate == "equal" else "REML with classical HKSJ" if aggregate == "re_hksj" else "conditional empirical Bayes",
                    patient_id_order="first observed, then additional patient_ids")
    add_reporting(summary, patients, metadata, max_ci_width=max_ci_width,
                  minimum_relevant_effect=minimum_relevant_effect,
                  planning_sd=planning_sd, planning_power=planning_power,
                  planning_sd_source=planning_sd_source)
    return AssociationResult(patients, summary, metadata)
