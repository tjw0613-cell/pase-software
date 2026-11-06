# Maintainer release instructions

## Release status

Version 1.0.0 is published on
[PyPI](https://pypi.org/project/pase-framework/1.0.0/) and
[GitHub](https://github.com/tjw0613-cell/pase-software/releases/tag/v1.0.0).
The [publication workflow](https://github.com/tjw0613-cell/pase-software/actions/runs/42162951531)
passed its cross-platform checks and PyPI upload. A separate Windows environment
installed the public PyPI wheel and ran the equal, RE-HKSJ and EB examples.
Read the Docs has not yet been registered or published, and Bioconda publication
has not been requested. Never commit API tokens.

The software repository is `tjw0613-cell/pase-software`; the manuscript
repository is `tjw0613-cell/pase-reproducibility`.
Keep each checkout in its own repository. Add the actual documentation URL
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
failures before tagging a release. For a future version, update both package
metadata and `pase.__version__`, then use the matching GitHub version tag.
Do not move the published `v1.0.0` tag or replace its distributions.

## PyPI

The published distribution name is **pase-framework**. `pase` is used by an
unrelated project and must not be substituted in installation commands.

The production Trusted Publisher uses GitHub owner `tjw0613-cell`, repository
`pase-software`, workflow `publish.yml` and environment `pypi`. For TestPyPI,
configure a separate publisher with environment `testpypi` before using the
manual testing workflow.

The manual **Publish release** workflow uploads to TestPyPI. A **published
GitHub release** triggers the same tests and then uploads to PyPI. The publisher
uses OIDC, not a password stored in the repository. After public upload,
verify the project page, download the published distribution and repeat the
five-line example in a fresh environment using:

```bash
python -m pip install pase-framework==1.0.0
```

Repeat the clean installation check for each future version before updating
its installation instructions. See the official
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

