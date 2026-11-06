# Preparing an AnnData object

## Expression and library size

`adata.X` (or `adata.layers[layer]`) contains transformed expression, with cells
as rows and genes as columns. Gene names must be unique. PASE extracts the
requested genes and densifies only that submatrix if the input is sparse.
It does not normalize, select variable genes, or modify the input.

For count data, calculate library sizes from all measured genes, then apply
CP10K and log1p. Do this before restricting the object to the anchor and targets:

```python
import numpy as np
from scipy import sparse

counts = adata.layers["counts"].astype(float).copy()
totals = np.asarray(counts.sum(axis=1)).ravel()
if np.any(totals <= 0):
    raise ValueError("Handle zero-library cells before normalization.")
if sparse.issparse(counts):
    transformed = counts.multiply((10000 / totals)[:, None]).tocsr()
    transformed.data = np.log1p(transformed.data)
else:
    transformed = np.log1p(counts * (10000 / totals)[:, None])
adata.layers["log1p_cp10k"] = transformed
adata.obs["log_depth"] = np.log1p(totals)
```

Use `layer="log1p_cp10k"` in the analysis. If a shared sample or library has
other normalization requirements, choose and document preprocessing before
calling PASE.

## Patients and cells

`patient_key` identifies the independent patient, not a specimen, sequencing
lane or individual cell. Resolve repeated specimens in the study design and
retain one inference unit per independent patient. Missing IDs are an error.

`cell_mask` is either a boolean `obs` column or a boolean vector in the same
order as the cells. Strings such as `"False"`, numeric 0/1 values and missing
mask values are rejected. Do not pre-subset the object if you want patients
with zero selected cells to remain visible. Use `patient_ids` to add patients
absent from the object entirely. All observed IDs remain in the output.

## Scores, covariates and missing values

A string target selects one transformed gene. A list creates a score by
standardizing each component within the patient's selected cells (sample
standard deviation, `ddof=1`) and averaging the z scores. A measured constant
gene contributes zero. A gene absent from `var_names` is an error.

The order is fixed:

1. Select cells using the mask.
2. Remove cells with nonfinite anchor or target expression.
3. Construct the patient's target score.
4. Remove cells with nonfinite covariates.
5. Residualize anchor and score using the same numeric covariates and intercept.
6. Correlate the residual vectors and apply patient eligibility rules.

No outcome-detection gate is applied. The input function accepts numeric
covariates, including explicit dummy variables, but does not encode categories.
Within each patient's complete cases, nuisance columns are centered and scaled
for numerical stability, and constant columns are omitted alongside the
intercept. Rank-deficient designs use their numerical rank. Nonfinite values are counted
in the patient output instead of being silently replaced with zero.

## Eligibility

The selected method requires at least `min_cells` complete cells and at least
three residual degrees of freedom. Constant residual vectors are not estimable.
Fisher methods additionally require `abs(r) < 1` and `n - k - 3 > 0`.
Here `k` is the nuisance design rank minus one (the intercept).

Equal and Fisher summaries can therefore have different eligible patients.
Inspect `equal_eligible`, `fisher_eligible` and both exclusion-reason columns
before comparing summaries. Small-cell patients are not automatically rescued
by weighting; the caller must lower `min_cells` explicitly and justify doing so.

