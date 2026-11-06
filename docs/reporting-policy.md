# Reporting policy

Version 1.1.0 uses the number of **eligible patients for the selected method**
to decide what to report. The attempted cohort size cannot satisfy this check.
These defaults are result-informed reporting choices, not validated guarantees.

| Condition | Default behavior |
|---|---|
| Fewer than 6 eligible patients | Retain observed patient correlations, exclusion reasons, descriptive mean/median and direction counts. Withhold inferential estimates, intervals and P values, including EB shrinkage intervals. |
| 6–19 eligible patients | Permit mathematically estimable inference and add a small-cohort warning. |
| 20 or more eligible patients | Permit mathematically estimable inference; calibration remains conditional on the method and data assumptions. |
| Full correlation-scale interval width exceeds 0.2 | Retain the interval and add a precision warning. |
| Planning MDE exceeds the minimum relevant effect | Add a model-based planning-power warning. |

`min_patients` is configurable, with a mathematical minimum of two. An explicit
value below six is flagged as an override. `max_ci_width` sets a practical
precision tolerance; `None` disables that check. Selecting intervals on their
observed width would change the reported subset, so a width warning does not
remove an interval. All results retain `calibration_guaranteed=False` in metadata.

## Evidence and timing

The completed count simulation contains 293 generating conditions and
3,302,000 study cohorts. The design and executable protocol were archived
before those formal runs as [protocol 2.0.0](https://doi.org/10.5281/zenodo.24010816).
Reporting floors of 2, 4 and 6 were registered illustrative comparisons.
Selecting six as the software default, adding a warning below twenty, and
choosing a width tolerance of 0.2 were development decisions made after
examining the results. They have not undergone independent confirmation.

With the primary cell cutoff of 30 and reporting floor of 6, 234 conditions
produced at least one reported interval. Across those conditions, equal-t
conditional coverage ranged from 81.76% to 96.78%; REML–HKSJ coverage ranged
from 71.68% to 96.01%. The method targets differ, so these ranges are not a
ranking of two estimators of one quantity. A stress condition with forty
attempted patients still had poor weighted-interval coverage. Increasing the
patient count cannot by itself certify inference.

The audit checked all 13,208 archived interval-shard hashes and reproduced the
registered floor-2/4/6 counts. This verifies the extraction of the policy
evidence, not the statistical adequacy of the chosen reporting rules.

The formal main workflow used a fixed 15-gene outcome and depth-only
adjustment. Arbitrary API gene sets, covariates, normalization choices or
outcome-dependent masks are outside a universal calibration claim. In
particular, `cell_mask` is applied before within-patient score standardization;
it does not reproduce the formal outcome gate applied after score construction.

## Planning MDE

For `aggregate="equal"`, specify an external `planning_sd` on the correlation
scale and record its origin in `planning_sd_source`. With eligible count
$J$, PASE solves for the noncentrality parameter $\lambda$ such that

$$
\Pr\bigl(\lvert T_{J-1,\lambda}\rvert > t_{J-1,1-\alpha/2}\bigr)
= \text{planning_power},\qquad
\mathrm{MDE}=\lambda\,\text{planning_sd}/\sqrt{J}.
$$

The default planning power is 0.8. `minimum_relevant_effect` is a scientifically
specified absolute effect against which to compare the MDE; supplying it
requires a planning SD. PASE cannot verify that a supplied SD came from an
independent source. If its source is omitted, the output records a warning.
Do not substitute the current observed effect for the scientific effect of
interest or label an observed sample SD as independent planning information.

The calculation assumes independent normal patient effects and tests a zero
arithmetic mean. It is not achieved power and is not a guarantee under the
mixture scenarios that caused undercoverage. An MDE above one is retained
and flagged as unattainable on the correlation scale. Fisher and EB MDEs
are not computed by applying the arithmetic-r formula to a different target.

Below the reporting floor, the MDE is suppressed. Above that floor, an external
planning calculation can remain available even if the observed-data interval
is technically undefined; `inference_available` stays false in that case.

## Weighted and individual summaries

Fisher-weighted inference always carries a model-assumption caution: its
variance approximation and relationship between cell yield and effect matter
even in larger cohorts. EB intervals condition on fitted population parameters
and are not cohort confidence intervals. No absence of a small-cohort or width
warning should be read as a declaration that either method is calibrated.
