# Reading the output

## Summary table

`result.summary()` returns a new one-row DataFrame. Its principal fields are:

| Field | Interpretation |
|---|---|
| `method`, `target` | Selected method and the population quantity it addresses. |
| `n_patients_total` | All observed patients plus any additional IDs provided. |
| `n_patients_eligible`, `n_patients_excluded` | Counts for the selected method. |
| `estimate`, `ci_low`, `ci_high` | Correlation-scale group result, or NaN when withheld or in EB mode. |
| `p_value` | Two-sided test of zero group effect, where returned. |
| `standard_error`, `se_scale` | Standard error and its scale; Fisher summaries use z-scale uncertainty. |
| `df` | Degrees of freedom used by the interval. |
| `tau2` | Estimated between-patient variance on the Fisher-z scale; not applicable to equal-r summaries. |
| `status`, `reason` | Whether inference was returned and why it was withheld. |
| `reporting_status`, `inference_available` | Descriptive, cohort-inferential or individual-EB reporting; availability refers only to cohort inference. |
| `descriptive_mean_r`, `descriptive_median_r` | Unshrunk correlations among patients eligible for the selected method, even when inference is withheld. |
| `n_positive`, `n_negative`, `n_zero` | Direction counts in that same eligible set. All other observed correlations remain in the patient table. |
| `warning_codes` | Machine-readable cautions; the complete explanations are in `metadata["warnings"]`. |
| `ci_width` | Full displayed correlation-scale cohort interval width; NaN if no cohort interval is returned. |
| `planning_mde`, `mde_status`, `mde_attainable` | External-SD equal-r planning calculation, its availability, and whether its magnitude is at most one. This is not achieved power or a full bounded-distribution feasibility check. |

Status `reported` denotes a computed interval, not demonstrated reliability.
`withheld_policy` means the selected patient threshold was not met.
`withheld_technical` means an interval was not mathematically estimable.
`individual_estimates_only` marks EB mode.

`reporting_status` uses `descriptive_only`, `inferential`,
`inferential_with_warnings`, `individual_estimates_only` or
`individual_estimates_with_warnings`. Even `inferential` does not certify
calibration; `metadata["calibration_guaranteed"]` is always false.

MDE status is `not_requested`, `withheld_policy`, `not_supported`,
`numerical_failure` or `model_based_planning`. A planning result can be
available with `inference_available=False` when the external planning inputs
are usable but the observed-data interval is technically undefined.

The principal warning codes are `below_min_patients`,
`below_default_min_patients`, `small_cohort_calibration_warning`,
`wide_interval`, `weighted_aggregation_assumptions`,
`fisher_variance_approximation`, `eb_conditional_uncertainty`,
`inference_unavailable`, `planning_sd_source_unspecified`,
`model_based_power_warning`, `mde_exceeds_correlation_range`,
`mde_not_supported` and `mde_numerical_failure`.
See [Reporting policy](reporting-policy.md) for their scope.

## Patient table

`result.patients` includes all patients, including those with no selected cells.
`n_total`, `n_selected`, `n_expression_complete` and `n_complete` follow the
processing steps. `n_cells` equals `n_complete`. `r` is the observed partial
correlation; it can be shown descriptively even when the patient falls below
the minimum cell count. The table includes design rank and residual degrees
of freedom.

`eligible`, `reason` and `status` refer to the selected summary. The
`equal_eligible`/`equal_reason` and `fisher_eligible`/`fisher_reason` columns
also reveal differences between the method-specific eligible sets.
`weight` is the normalized group weight when available. In EB mode,
`eb_data_weight` is the patient's shrinkage coefficient, not a group weight.

EB fields include `eb_estimate`, `eb_ci_low`, `eb_ci_high`,
`eb_conditional_se_z` and `eb_status`. A missing EB interval at zero estimated
heterogeneity is intentional. Raw patient correlations do not receive
individual intervals in the equal or random-effects forest plot.

## Saving tables and metadata

```python
import json
result.summary().to_csv("summary.tsv", sep="\t", index=False)
result.patients.to_csv("patients.tsv", sep="\t", index=False)
with open("settings.json", "w", encoding="utf-8") as handle:
    json.dump(result.metadata, handle, indent=2)
```

