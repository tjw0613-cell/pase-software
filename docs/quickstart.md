# Quick start

Create and activate an environment, then install the published PASE release:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install "pase-framework[tutorial]==1.1.0"
```

The package distribution is `pase-framework`; the Python module is `pase`.
Use a fresh environment to avoid the unrelated `pase` audio package.

```python
import pase
adata = pase.datasets.load_gse234933_demo()
result = pase.associate(
    adata, anchor="LILRB4", target=["MMP14", "MMP9", "SPP1", "APOE"],
    patient_key="donor", covariates=["log_depth"], cell_mask="myeloid",
    aggregate="equal", min_cells=30,
)
print(result.summary().to_string(index=False))
print(result.patients[["patient_id", "n_cells", "r", "eligible", "reason"]])
ax = result.plot_forest()
ax.figure.savefig("association.pdf", bbox_inches="tight")
```

This is an interface example using real cells, not a full-cohort analysis.
The four-gene target demonstrates module input; the revised manuscript's main
biological analysis uses a separately sourced 15-gene signature.

Download [the tutorial notebook](https://github.com/tjw0613-cell/pase-software/blob/v1.1.0/examples/gse234933_quickstart.ipynb)
or clone the software repository. Open `examples/gse234933_quickstart.ipynb`
with JupyterLab for a complete
walkthrough of normalization metadata, exclusions, the three summaries,
threshold sensitivity and saving results. After installation, that notebook
needs neither network access nor a separate raw-data download.


## Inspect reporting warnings

The teaching subset has eight donors. Its small-cohort warning is expected.
Check `result.metadata["warnings"]` and the eligible patient count rather than
treating a computed interval as proof of reliability. The default reporting
floor is six eligible patients; the warning boundary is twenty.

For an externally specified planning scenario, use:

```python
planned = pase.associate(
    adata, anchor="LILRB4", target=["MMP14", "MMP9", "SPP1", "APOE"],
    patient_key="donor", covariates=["log_depth"], cell_mask="myeloid",
    aggregate="equal", planning_sd=0.15,
    planning_sd_source="Illustrative planning assumption, not an empirical estimate",
    minimum_relevant_effect=0.10, planning_power=0.80,
)
print(planned.summary())
```

The SD and effect above are illustrative user choices, not values established
for this biological association. Replace them with justified study-specific
planning values. Do not use an observed effect as the minimum relevant effect.
