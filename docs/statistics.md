# Statistical meaning

## Within-patient association

For patient $i$, PASE fits separate least-squares nuisance models to the anchor
and target scores, using the same design matrix including an intercept.
The Pearson correlation of the two residual vectors is $r_i$.
The inferential units are eligible patients; cells define each patient's
observed association.

## Equal-patient mean

For $J$ eligible patients, `aggregate="equal"` reports

$$
\bar r = J^{-1}\sum_{i=1}^J r_i,\qquad
\mathrm{SE}(\bar r) = s_r/\sqrt{J},\qquad
\bar r \pm t_{J-1,1-\alpha/2}\,s_r/\sqrt{J}.
$$

The target is the mean observed correlation among patients meeting the chosen
eligibility rule. It is not a cell-weighted correlation or the correlation of
pseudobulk patient means. The t interval assumes independent patient effects;
small-sample calibration depends on their distribution. Intervals are not
silently clipped to [-1, 1]. If the across-patient standard error is zero,
the interval and P value are withheld.

## Random-effects Fisher summary

`aggregate="re_hksj"` transforms $z_i=\operatorname{atanh}(r_i)$ and uses the
approximate sampling variance

$$
v_i = (n_i-k_i-3)^{-1},\qquad
w_i = (v_i+\hat\tau^2)^{-1},\qquad
\hat\mu = \frac{\sum_i w_i z_i}{\sum_i w_i}.
$$

$n_i$ is the complete-cell count; $k_i$ is the nuisance design rank minus one.
The between-patient variance $\tau^2$ is fitted by restricted maximum likelihood
(REML), with its boundary at zero included. Under the random-effects model,
the target is the population mean on the Fisher-z scale. The displayed point
estimate is $\tanh(\hat\mu)$, which is not the arithmetic mean of $r_i$.

The **classical HKSJ** interval uses

$$
q = \frac{\sum_i w_i(z_i-\hat\mu)^2}{J-1},\qquad
\hat\mu \pm t_{J-1,1-\alpha/2}\sqrt{q/\sum_i w_i}.
$$

Both endpoints are transformed back with tanh. The two-sided P value tests
zero mean on the z scale. This API does not substitute `max(q, 1)` for $q$;
it is not the modified HKSJ procedure. The interval can be narrower than the
corresponding Wald interval when $q<1$. A zero standard error is withheld.

The Fisher variance approximation assumes independent cell-level information
and a suitable correlation model. Count noise, dependence within cell
subpopulations and selection can violate these assumptions. HKSJ does not
remove those violations or validate inference at an arbitrary patient count.
If cell yields are related to effects, weighted and equal-patient summaries
can diverge for substantive as well as statistical reasons.

## Individual empirical Bayes estimates

`aggregate="eb"` fits the same REML model and computes

$$
B_i=\frac{\hat\tau^2}{\hat\tau^2+v_i},\qquad
\tilde z_i=\hat\mu+B_i(z_i-\hat\mu),\qquad
\operatorname{Var}(\theta_i\mid z_i,\hat\mu,\hat\tau^2)=B_i v_i.
$$

Here $\theta_i$ is the latent patient effect on the Fisher-z scale, and
$\tilde z_i$ is its conditional posterior mean. Normal conditional intervals
are formed on the z scale and transformed back.
They treat estimated population parameters as fixed, omit their uncertainty,
and should not be interpreted as fully calibrated small-sample intervals.
When $\hat\tau^2=0$, point estimates pool completely and individual intervals
are withheld. EB estimates are in `res.patients`; `res.summary()` deliberately
contains no cohort estimate, cohort confidence interval or P value for this mode.

## Reporting policy

The default `min_patients=6` is a reporting floor applied after method-specific
eligibility checks. Below it, patient correlations and descriptive summaries
remain available, while inferential quantities are withheld. The result-informed
policy and its limitations are detailed in [Reporting policy](reporting-policy.md).
Small-cohort and wide-interval warnings do not remove otherwise computable
intervals. Passing either check does not establish calibration.

For the equal-patient mean only, a supplied external `planning_sd` enables a
two-sided noncentral-t MDE calculation at the requested `alpha` and
`planning_power`. This assumes independent normally distributed patient effects
with the specified planning SD. It does not estimate achieved power from the
observed effect, guarantee power under mixture distributions, or support the
Fisher-z or EB target. The API does not provide a prediction interval.

The package does not adjust across multiple calls: when testing many
associations, define the testing family and apply the selected multiplicity
correction separately.
