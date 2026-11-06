# Software tests and manuscript reproduction

This repository supplies the software API, tests and a teaching dataset.
The separate [pase-reproducibility repository](https://github.com/tjw0613-cell/pase-reproducibility) supplies manuscript figure
scripts and the inputs used by those scripts.

## Test the software

```bash
python -m pip install ".[test]"
python -m pytest
```

Tests cover independent calculation examples, sparse/dense parity, patient
exclusions, cell duplication, model boundaries, input validation and lack of
input mutation. The CI configuration requests Linux, macOS and Windows runs
across Python versions. A workflow file is not evidence that those remote runs
have passed; use the actual run logs for that claim.

## Run the tutorial

```bash
python -m pip install ".[tutorial]"
python scripts/run_tutorial.py --output tutorial-output
```

The runner launches a kernel using the current Python executable and writes an
executed notebook and measured runtime. Installation time and reconstruction
of the teaching subset from original study files are separate from notebook
execution time. The runner fails if execution exceeds ten minutes or any cell
raises an error.

## Rebuild paper figures

In the separate repository, `python reproduce.py fig2` is a portable equivalent
of `make fig2`. Each figure command documents whether it recreates a plot from
archived summary inputs or repeats the upstream analysis. Missing formal
simulation results are reported as missing, not replaced by synthetic charts.
The bundled notebook does not reproduce all manuscript results.

