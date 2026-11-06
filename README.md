# PASE

PASE estimates molecular associations within individual patients and summarizes them across a single-cell cohort.

## Installation

Use Python 3.10 or later in a fresh environment. Install the published release
from [PyPI](https://pypi.org/project/pase-framework/1.0.0/):

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install pase-framework==1.0.0
```

The distribution is named **pase-framework**, but the Python import is `pase`.
The unrelated audio package named `pase` on PyPI uses the same import name;
do not install both in one environment. The
[v1.0.0 release](https://github.com/tjw0613-cell/pase-software/releases/tag/v1.0.0)
is also available on GitHub. For development, clone this repository and use
`python -m pip install -e ".[test]"`.

## Five-line example

```python
import pase
adata = pase.datasets.load_gse234933_demo()
res = pase.associate(adata, anchor="LILRB4", target=["MMP14", "MMP9", "SPP1", "APOE"], patient_key="donor", covariates=["log_depth"], cell_mask="myeloid", aggregate="equal", min_cells=30)
print(res.summary())
res.plot_forest()
```

The bundled GSE234933 subset contains real cells and works without a download.
The four-gene target above illustrates the interface; it is not the current
manuscript's published 15-gene signature. The complete tutorial is
[examples/gse234933_quickstart.ipynb](examples/gse234933_quickstart.ipynb).
Install the notebook tools with
`python -m pip install "pase-framework[tutorial]==1.0.0"`. Download the notebook
or clone this repository, then run
`jupyter lab examples/gse234933_quickstart.ipynb` from the checkout.

## Input format

Supply an AnnData object with cells in rows and unique gene names in
`adata.var_names`. `adata.X`, or the requested `layer`, must already contain
transformed expression. For count data, the tutorial uses
`log1p(10000 * count / total_counts)` with totals calculated **before** selecting
genes. `associate` does not normalize raw counts.

| Input | Meaning |
|---|---|
| `adata.obs["donor"]` | Independent patient ID, with no missing values. Repeated samples from one patient must share an ID. |
| `adata.obs["log_depth"]` | Numeric covariate; the tutorial uses `log1p(total_counts)` from the full matrix. |
| `adata.obs["myeloid"]` | Boolean cell-selection column. Patients with no selected cells remain in the output. |
| `anchor` | One gene name. |
| `target` | One gene name, or a list whose within-patient gene z scores are averaged. |

Covariates enter a separate least-squares model for each patient, with an
intercept. Numeric categorical encodings must be supplied explicitly.
Patients absent from the object can be retained using `patient_ids`.
See [the input guide](docs/inputs.md) for missing values and eligibility rules.

## Outputs and choice of summary

`res.patients` reports every patient, cell counts, observed correlation,
eligibility, exclusion reasons and, where applicable, weights or shrunken effects.
`res.summary()` returns one row with the selected method, target, estimate,
interval, P value and status. `res.metadata` records the analysis settings and
warnings. `res.plot_forest()` returns a Matplotlib Axes.

| `aggregate` | What it reports |
|---|---|
| `"equal"` | Arithmetic mean of eligible patients' observed correlations, with a t interval. Each eligible patient has equal weight. |
| `"re_hksj"` | REML random-effects mean on the Fisher-z scale, back-transformed for display, with a classical HKSJ interval. It has a different target from the arithmetic mean. |
| `"eb"` | Individual patient effects shrunk toward the fitted Fisher-z population mean. Conditional intervals exclude uncertainty in fitted population parameters; no cohort P value is returned. |

`min_cells=30` is adjustable. `min_patients=2` is a technical minimum, **not a
validated recommendation for reliable inference**. Set a stricter reporting
threshold explicitly, for example `min_patients=6`. When that threshold is not
met, cohort estimates, intervals and P values are withheld, and the patient
records remain available. A configurable threshold does not establish its
statistical adequacy. See [statistical details](docs/statistics.md).

## Citation

For this code, cite the authors, PASE software version **1.0.0**, and its
[release URL](https://github.com/tjw0613-cell/pase-software/releases/tag/v1.0.0).
Machine-readable author
information is in [CITATION.cff](CITATION.cff). An earlier analysis archive is
available at [Zenodo v3](https://doi.org/10.5281/zenodo.23032874); that DOI does
not identify this 1.0.0 software release.

Code is distributed under the MIT License. The real-data example retains its
source attribution and third-party terms; see [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md).

