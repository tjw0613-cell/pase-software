# Quick start

From the source checkout, create and activate an environment and install PASE:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install ".[tutorial]"
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

Open `examples/gse234933_quickstart.ipynb` with JupyterLab for a complete
walkthrough of normalization metadata, exclusions, the three summaries,
threshold sensitivity and saving results. After installation, that notebook
needs neither network access nor a separate raw-data download.

