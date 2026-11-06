"""Reporting safeguards and explicitly model-based planning calculations.

This module does not modify the association estimates or interval kernels.
"""
from __future__ import annotations

from numbers import Real

import numpy as np
from scipy import optimize, stats


DEFAULT_MIN_PATIENTS = 6
SMALL_COHORT_WARNING_BELOW = 20
VALIDATION_SCOPE = (
    "A simulation-informed reporting policy, not a calibration guarantee. "
    "The default floor of six was a registered reporting candidate, not an "
    "independently validated minimum for reliable inference. The warning below "
    "20 eligible patients is a post-result risk flag; reaching 20 does not certify "
    "coverage. The formal main simulation used a fixed 15-gene score and "
    "depth-only adjustment. Arbitrary modules, covariates and masks in this API "
    "are not covered by that validation. This API standardizes modules after "
    "cell_mask selection; the formal outcome-gating experiment standardized "
    "before its gate. The configurable interval-width tolerance is a reporting "
    "preference, not a simulation-calibrated cutoff."
)
MDE_SCOPE = (
    "Planning calculation for a two-sided one-sample t test of a zero equal-patient "
    "mean correlation, using independent, identically distributed normal patient "
    "effects and an externally specified r-scale standard deviation. Solves "
    "noncentral-t power at the eligible patient count and requested alpha. "
    "It does not estimate power from the observed effects, guarantee power for "
    "heterogeneous/count-data distributions, or apply to Fisher-scale or EB targets. "
    "Attainability checks only whether the required magnitude exceeds one; it "
    "does not validate joint mean/variance feasibility for bounded correlations."
)


def _positive_real(value, name, *, upper=None, optional=True):
    if value is None and optional:
        return None
    if (isinstance(value, (bool, np.bool_)) or not isinstance(value, Real)
            or not np.isfinite(value) or value <= 0
            or (upper is not None and value > upper)):
        bound = f" and <= {upper}" if upper is not None else ""
        raise ValueError(f"{name} must be a finite number > 0{bound}.")
    return float(value)


def validate_planning(max_ci_width, minimum_relevant_effect, planning_sd,
                      planning_power, alpha, planning_sd_source):
    """Reject ambiguous or impossible reporting inputs before fitting any data."""
    width = _positive_real(max_ci_width, "max_ci_width")
    effect = _positive_real(minimum_relevant_effect, "minimum_relevant_effect", upper=1)
    sd = _positive_real(planning_sd, "planning_sd", upper=1)
    power = _positive_real(planning_power, "planning_power", optional=False)
    if power >= 1:
        raise ValueError("planning_power must be strictly between zero and one.")
    if sd is not None and power <= alpha:
        raise ValueError("planning_power must exceed alpha when planning_sd is provided.")
    if effect is not None and sd is None:
        raise ValueError("minimum_relevant_effect requires an externally specified planning_sd.")
    if planning_sd_source is not None:
        if not isinstance(planning_sd_source, str) or not planning_sd_source.strip():
            raise ValueError("planning_sd_source must be nonempty text or None.")
        if sd is None:
            raise ValueError("planning_sd_source requires planning_sd.")
        planning_sd_source = planning_sd_source.strip()
    return width, effect, sd, power, planning_sd_source


def _planning_mde(n_patients, sd, alpha, power):
    """Invert two-sided noncentral-t power on the noncentrality scale."""
    df = n_patients - 1
    critical = stats.t.isf(alpha / 2, df)
    if not np.isfinite(critical):
        raise ValueError("Nonfinite t critical value in planning calculation.")

    def difference(noncentrality):
        # The reflected survival function equals the left-tail probability.
        # It avoids NaN from some SciPy negative-tail CDF evaluations at small
        # df and large noncentrality during root bracketing.
        attained = (stats.nct.sf(critical, df, -noncentrality)
                    + stats.nct.sf(critical, df, noncentrality))
        if not np.isfinite(attained):
            raise ValueError("Nonfinite noncentral-t power in planning calculation.")
        return float(attained) - power

    upper = 1.0
    while difference(upper) < 0:
        upper *= 2
        if upper > 1e8:
            raise ValueError("Could not bracket the planning noncentrality.")
    noncentrality = optimize.brentq(difference, 0.0, upper, xtol=1e-12, rtol=1e-12)
    mde = float(noncentrality * (sd / np.sqrt(n_patients)))
    if not np.isfinite(mde) or mde <= 0:
        raise ValueError("Planning MDE could not be represented as a positive finite number.")
    return mde


def add_reporting(summary, patients, metadata, *, max_ci_width,
                  minimum_relevant_effect, planning_sd, planning_power,
                  planning_sd_source):
    """Annotate, without changing, already fitted method-specific results."""
    count = summary["n_patients_eligible"]
    method = summary["method"]
    floor = summary["min_patients"]
    reported = summary["status"] == "reported"
    individual = summary["status"] == "individual_estimates_only"
    effects = patients.loc[patients.eligible, "r"].to_numpy(dtype=float)
    codes, messages = [], []

    def warn(code, message):
        codes.append(code)
        messages.append(message)

    summary.update(
        inference_available=reported,
        descriptive_mean_r=float(np.mean(effects)) if len(effects) else np.nan,
        descriptive_median_r=float(np.median(effects)) if len(effects) else np.nan,
        n_positive=int(np.count_nonzero(effects > 0)),
        n_negative=int(np.count_nonzero(effects < 0)),
        n_zero=int(np.count_nonzero(effects == 0)),
        ci_width=np.nan, planning_mde=np.nan, mde_status="not_requested",
        mde_attainable=None,
    )
    if floor < DEFAULT_MIN_PATIENTS:
        warn("below_default_min_patients",
             "Expert override: min_patients is below the default reporting floor of six.")
    if count < floor:
        warn("below_min_patients",
             "Only descriptive results are returned below the chosen eligible-patient floor.")
    elif (reported or individual) and count < SMALL_COHORT_WARNING_BELOW:
        warn("small_cohort_calibration_warning",
             "Fewer than 20 eligible patients: coverage failures occurred in some small-cohort "
             "simulation settings. This is a risk flag, not a calibrated cutoff.")
    if method != "equal":
        warn("fisher_variance_approximation",
             "Fisher sampling variances approximate independent cell-level information.")
        warn("weighted_aggregation_assumptions",
             "If cell yield and patient association strength are related, precision weighting "
             "can change the represented patient population. Independence is a model assumption "
             "and is not checked by this fit.")
    if method == "eb" and individual:
        warn("eb_conditional_uncertainty",
             "EB intervals hold fitted population parameters fixed and omit their uncertainty.")
    if summary["status"] == "withheld_technical":
        warn("inference_unavailable", "Inference unavailable: " + summary["reason"].replace("_", " ") + ".")
    if reported:
        width = float(summary["ci_high"] - summary["ci_low"])
        summary["ci_width"] = width
        if max_ci_width is not None and width > max_ci_width:
            warn("wide_interval", "The correlation-scale interval exceeds the chosen width tolerance.")

    if planning_sd is not None:
        if planning_sd_source is None:
            warn("planning_sd_source_unspecified",
                 "The planning SD has no user-declared source. Supply an independent external assumption; "
                 "the software cannot verify its origin.")
        if count < floor:
            summary["mde_status"] = "withheld_policy"
        elif method != "equal":
            summary["mde_status"] = "not_supported"
            warn("mde_not_supported", "MDE planning is implemented only for equal-patient mean correlations.")
        else:
            try:
                mde = _planning_mde(count, planning_sd, metadata["alpha"], planning_power)
            except (ValueError, RuntimeError, FloatingPointError, OverflowError):
                summary["mde_status"] = "numerical_failure"
                warn("mde_numerical_failure", "The specified noncentral-t planning calculation could not be evaluated.")
            else:
                summary.update(planning_mde=mde, mde_status="model_based_planning",
                               mde_attainable=bool(mde <= 1))
                if mde > 1:
                    warn("mde_exceeds_correlation_range", "The planning MDE exceeds the maximum correlation magnitude of one.")
                if minimum_relevant_effect is not None and mde > minimum_relevant_effect:
                    warn("model_based_power_warning",
                         "Under the supplied planning model, the relevant effect is smaller than the MDE.")

    summary["warning_codes"] = tuple(codes)
    if reported:
        summary["reporting_status"] = "inferential_with_warnings" if codes else "inferential"
    elif individual:
        summary["reporting_status"] = "individual_estimates_with_warnings" if codes else "individual_estimates_only"
    else:
        summary["reporting_status"] = "descriptive_only"
    metadata.update(
        calibration_guaranteed=False, validation_scope=VALIDATION_SCOPE,
        warning_codes=list(codes), max_ci_width=max_ci_width,
        minimum_relevant_effect=minimum_relevant_effect,
        planning_sd=planning_sd, planning_power=planning_power,
        planning_sd_source=planning_sd_source, planning_sd_source_verified=False,
        mde_scope=MDE_SCOPE,
        mde_null_mean_r=0.0,
        mde_sd_source=("user-declared: " + planning_sd_source if planning_sd_source is not None
                       else "source unspecified" if planning_sd is not None else "not provided"),
        mde_patient_count=count, mde_alpha=metadata["alpha"],
        descriptive_population="Patients eligible for the selected aggregation method; raw correlations, before EB shrinkage",
        reporting_policy={
            "default_min_patients": DEFAULT_MIN_PATIENTS,
            "min_patients": floor,
            "custom_min_patients": floor != DEFAULT_MIN_PATIENTS,
            "small_cohort_warning_below": SMALL_COHORT_WARNING_BELOW,
            "width_tolerance": max_ci_width,
            "protocol_doi": "10.5281/zenodo.24010816",
            "policy_basis": "registered candidate floor; post-result risk warning; configurable width preference",
        },
    )
    metadata["warnings"].extend(messages)
