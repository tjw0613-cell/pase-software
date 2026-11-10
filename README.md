# PASE

PASE estimates molecular associations within individual patients and summarizes them across a single-cell cohort.

Read the [online documentation](https://pase-software.readthedocs.io/en/latest/)
for installation, the API and the statistical meaning of each parameter.

## Installation

Use Python 3.10 or later in a fresh environment. Install the published release
from [PyPI](https://pypi.org/project/pase-framework/1.1.1/):

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install pase-framework==1.1.1
```

macOS or Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install pase-framework==1.1.1
```

These commands use the new environment directly; activation is not required.
Version 1.1.1 updates the documentation and packaging only. Statistical code
and bundled data are unchanged from 1.1.0, which was used for the manuscript.

The distribution is named **pase-framework**, but the Python import is `pase`.
The unrelated audio package named `pase` on PyPI uses the same import name;
do not install both in one environment. The
[v1.1.1 release](https://github.com/tjw0613-cell/pase-software/releases/tag/v1.1.1)
is also available on GitHub. For development, clone this repository and use
`python -m pip install -e ".[test]"`.

## Quick example

Save this as `quickstart.py`:

```python
import pase
adata = pase.datasets.load_gse234933_demo()
res = pase.associate(
    adata, anchor="LILRB4", target=["MMP14", "MMP9", "SPP1", "APOE"],
    patient_key="donor", covariates=["log_depth"], cell_mask="myeloid",
    aggregate="equal", min_cells=30,
)
print(res.summary().to_string(index=False))
print(res.patients[["patient_id", "n_cells", "r", "eligible", "reason"]])
ax = res.plot_forest()
ax.figure.savefig("association.pdf", bbox_inches="tight")
```

Run it with `.\.venv\Scripts\python.exe quickstart.py` on Windows or
`.venv/bin/python quickstart.py` on macOS/Linux. It prints the cohort summary
and patient table and saves a forest plot as `association.pdf`.

The bundled GSE234933 subset contains real cells and works without a download.
The four-gene target above illustrates the interface; it is not the current
manuscript's 15-gene signature. The
[tutorial notebook](https://github.com/tjw0613-cell/pase-software/blob/v1.1.1/examples/gse234933_quickstart.ipynb)
uses the 15-gene list on this small teaching subset. For the complete cohort
analyses, use the [RNA workflow](https://github.com/tjw0613-cell/pase-reproducibility/tree/main/workflows/rna).
See the [quick-start guide](https://pase-software.readthedocs.io/en/latest/quickstart.html)
for notebook installation and execution.

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
intercept. Encode categorical covariates as numeric dummy variables; integer
labels such as 0, 1 and 2 are otherwise treated as a numeric covariate.
Patients absent from the object can be retained using `patient_ids`.
See [the input guide](https://pase-software.readthedocs.io/en/latest/inputs.html)
for missing values and eligibility rules.

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

REML means restricted maximum likelihood; HKSJ is the
Hartung–Knapp–Sidik–Jonkman interval adjustment; EB means empirical Bayes.
See [the output guide](https://pase-software.readthedocs.io/en/latest/outputs.html)
for field definitions and warning codes.

`min_cells=30` is adjustable. Version 1.1.0 changes the default
`min_patients` from 2 to **6**. Below the selected method's eligible-patient
threshold, PASE returns observed patient correlations, descriptive summaries
and exclusion reasons, but withholds inferential estimates, intervals and
P values. An explicit lower threshold is allowed and recorded as an override.

Reported inference receives a small-cohort warning below 20 eligible patients
and a precision warning when the full correlation-scale interval is wider
than `max_ci_width=0.2`. Warnings retain the interval. The patient boundaries
are result-informed reporting choices, and 0.2 is a configurable practical
tolerance; none guarantees calibrated coverage. See the
[reporting policy and its evidence](https://pase-software.readthedocs.io/en/latest/reporting-policy.html).

For the equal-patient mean, an externally specified `planning_sd` enables a
model-based minimum detectable effect (MDE). Record its origin with
`planning_sd_source`; set `minimum_relevant_effect` to compare the MDE with a
scientifically meaningful effect. This is a planning calculation, not achieved
power, and is not provided for Fisher or EB targets.

## Reproduce the manuscript

The [reproduction repository](https://github.com/tjw0613-cell/pase-reproducibility)
contains figure commands and the full 15-gene RNA workflow. The software
tutorial above is a smaller example, not a reproduction of all paper results.

- [Simulation protocol, version 2.0.0](https://doi.org/10.5281/zenodo.24010816):
  the design archived before the formal simulations.
- [Core results and reproducibility materials](https://doi.org/10.5281/zenodo.24015713):
  the subsequent code and results archive, including PASE 1.1.0.
- Supplementary Data S1 and S2 accompany the revised manuscript. They contain
  the detection-matched gene analysis and unequal-yield native count benchmark,
  respectively, and are not included in the core-results DOI above.

See [reproduction instructions](https://pase-software.readthedocs.io/en/latest/reproducibility.html)
for the scope of each workflow and the software test commands.

## Citation

For this code, cite the authors, PASE software version **1.1.1**, and its
[release URL](https://github.com/tjw0613-cell/pase-software/releases/tag/v1.1.1).
Machine-readable author information is in
[CITATION.cff](https://github.com/tjw0613-cell/pase-software/blob/v1.1.1/CITATION.cff).
An earlier analysis archive is
available at [Zenodo v3](https://doi.org/10.5281/zenodo.23032874); that DOI does
not identify this 1.1.1 software release.

Code is distributed under the MIT License. The real-data example retains its
source attribution and third-party terms; see
[THIRD_PARTY_DATA.md](https://github.com/tjw0613-cell/pase-software/blob/v1.1.1/THIRD_PARTY_DATA.md).
Report software problems through [GitHub Issues](https://github.com/tjw0613-cell/pase-software/issues).

