# Maintainer release instructions

## What must be configured outside this checkout

The source is prepared as version 1.0.0. A local wheel, CI configuration or
documentation build is not a public software release. The maintainer must
connect actual GitHub, PyPI and Read the Docs accounts. Never commit API tokens.

The software repository is `tjw0613-cell/pase-software`; the manuscript
repository is `tjw0613-cell/pase-reproducibility`.
Push each checkout to its own repository. Add the actual documentation URL
to `pyproject.toml` once the documentation site exists.
Do not reuse the earlier manuscript archive DOI as this version's software DOI.

## Build and check

```bash
python -m pip install '.[test,docs,tutorial,build]'
python -m pytest
python scripts/run_tutorial.py
python -m sphinx -W --keep-going -b html docs docs/_build/html
python -m build
python -m twine check dist/*
```

Inspect the actual GitHub Actions results on Linux, macOS and Windows. Fix
failures before tagging the final release. Version 1.0.0 is already set in
both package metadata and `pase.__version__`; the GitHub release tag must be
`v1.0.0`.

## PyPI

The distribution name is **pase-framework**. `pase` is already used by an
unrelated project. A name lookup does not reserve a name; registration remains
subject to PyPI's checks.

Configure a pending Trusted Publisher for the new project in PyPI using the
GitHub owner `tjw0613-cell`, repository `pase-software`, workflow `publish.yml` and
environment `pypi`. Configure the corresponding TestPyPI publisher with
environment `testpypi`. Add those environments to the GitHub repository.

The manual **Publish release** workflow uploads to TestPyPI. A **published
GitHub release** triggers the same tests and then uploads to PyPI. The publisher
uses OIDC, not a password stored in the repository. After public upload,
verify the project page, download the published distribution and repeat the
five-line example in a fresh environment using:

```bash
python -m pip install pase-framework==1.0.0
```

Update the README installation section only after that command has succeeded
against the intended public index. See the official
[Trusted Publisher guide](https://docs.pypi.org/trusted-publishers/using-a-publisher/).

## Read the Docs

Import the real software repository into Read the Docs. The root
`.readthedocs.yaml` selects Python, installs the documentation extra and builds
Sphinx with warnings treated as errors. Activate the version tag and verify the
public API pages after the build succeeds. The project URL is assigned by
Read the Docs; no live URL is claimed in this source checkout.
See the [configuration reference](https://docs.readthedocs.com/platform/stable/config-file/v2.html).

## Bioconda

After the PyPI source distribution is public, run
`python scripts/prepare_bioconda.py --source-url URL --maintainer GITHUB_NAME --homepage REPOSITORY_URL`.
The script downloads the exact source artifact, verifies its package/version,
and writes a recipe with its SHA-256. Supply a maintainer who has agreed to the
role. Review and submit the recipe to the Bioconda recipes repository and run
its required checks. Recipe generation is not a channel publication.
See [Bioconda contributor guidelines](https://bioconda.github.io/contributor/guidelines.html).

