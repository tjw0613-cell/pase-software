# Quick start

Install the published release in a new environment. These commands call that
environment's Python directly, so no activation step is needed.

Windows PowerShell:

```powershell
python -m venv .venv
.\.venv\Scripts\python.exe -m pip install "pase-framework[tutorial]==1.1.1"
```

macOS or Linux:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install "pase-framework[tutorial]==1.1.1"
```

The package distribution is `pase-framework`; the Python module is `pase`.
Use a fresh environment to avoid the unrelated `pase` audio package.

Save the following example as `quickstart.py`. Run it with
`.\.venv\Scripts\python.exe quickstart.py` on Windows or
`.venv/bin/python quickstart.py` on macOS/Linux.

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

Download [the tutorial notebook](https://github.com/tjw0613-cell/pase-software/blob/v1.1.1/examples/gse234933_quickstart.ipynb)
into the working directory, then open it with JupyterLab:

```powershell
.\.venv\Scripts\python.exe -m jupyter lab gse234933_quickstart.ipynb
```

On macOS/Linux, replace `.\.venv\Scripts\python.exe` with `.venv/bin/python`.
If you cloned the repository instead, use
`examples/gse234933_quickstart.ipynb` as the notebook path. The notebook uses
the 15-gene target and walks through normalization metadata, exclusions, the
three summaries, threshold sensitivity and saving results. After installation,
it needs neither network access nor a separate raw-data download. The
[full-cohort RNA workflow](https://github.com/tjw0613-cell/pase-reproducibility/tree/main/workflows/rna)
is separate; its inputs and computing requirements are larger.


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
