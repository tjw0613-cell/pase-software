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

Status `reported` denotes a computed interval, not demonstrated reliability.
`withheld_policy` means the selected patient threshold was not met.
`withheld_technical` means an interval was not mathematically estimable.
`individual_estimates_only` marks EB mode.

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

