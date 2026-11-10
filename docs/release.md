# Maintainer release instructions

## Release records

Version 1.1.1 corrects documentation and packaging. Statistical code, reporting
defaults and bundled data are unchanged from 1.1.0, the version used for the
manuscript analyses.

Version 1.1.0 changes reporting defaults and adds precision warnings and a
model-based planning MDE. See [the changelog](https://github.com/tjw0613-cell/pase-software/blob/v1.1.1/CHANGELOG.md)
and [reporting policy](reporting-policy.md) before upgrading an existing workflow.
Release artifacts are indexed on [PyPI](https://pypi.org/project/pase-framework/1.1.1/)
and [GitHub](https://github.com/tjw0613-cell/pase-software/releases/tag/v1.1.1).
Check the actual [release workflow results](https://github.com/tjw0613-cell/pase-software/actions/workflows/publish.yml)
for cross-platform tests and upload status; the configuration alone is not proof
that a run passed.

The software repository is `tjw0613-cell/pase-software`; the manuscript
repository is `tjw0613-cell/pase-reproducibility`. A Bioconda package has not
been published. The earlier manuscript archive DOI and the simulation protocol
DOI do not identify this software version. Never commit API tokens or move an
existing published release tag.

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
failures before tagging a release. For a future version, update `pyproject.toml`,
`pase.__version__` and `CITATION.cff`, then use the matching GitHub version tag.
Keep installation commands, release links and the changelog consistent with
the version being published. Do not move an existing published tag or replace
its distributions.

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
quick example in a fresh environment using:

```bash
python -m pip install pase-framework==1.1.1
```

Repeat the clean installation check for each future version before updating
its installation instructions. See the official
[Trusted Publisher guide](https://docs.pypi.org/trusted-publishers/using-a-publisher/).

PyPI displays the README packaged with a release. Editing the GitHub README
does not update an already published distribution. A documentation-only patch
release can carry corrected instructions without changing statistical code;
state that scope in its changelog and retain the earlier release. PyPI does
not allow an existing distribution filename to be reused. See
[PyPI upload troubleshooting](https://docs.pypi.org/publishing/troubleshooting/).

## Read the Docs

The project is [pase-software](https://app.readthedocs.org/projects/pase-software/).
Its public documentation is available at
https://pase-software.readthedocs.io/en/latest/.
The root `.readthedocs.yaml` selects Python, installs the documentation extra
and builds Sphinx with warnings treated as errors. The initial successful
build was [35689123](https://app.readthedocs.org/projects/pase-software/builds/35689123/),
from commit `f15b4ad5748de1279d0b5b640f07b1db4f51a530`.
The `latest` documentation follows `main`; it is not an immutable version archive.
For future tagged documentation, activate the tag and verify its own public URL
after a successful build. Do not move an existing release tag.
See the [configuration reference](https://docs.readthedocs.com/platform/stable/config-file/v2.html).

## Bioconda

After the PyPI source distribution is public, run
`python scripts/prepare_bioconda.py --source-url URL --maintainer GITHUB_NAME --homepage REPOSITORY_URL`.
The script downloads the exact source artifact, verifies its package/version,
and writes a recipe with its SHA-256. Supply a maintainer who has agreed to the
role. Review and submit the recipe to the Bioconda recipes repository and run
its required checks. Recipe generation is not a channel publication.
See [Bioconda contributor guidelines](https://bioconda.github.io/contributor/guidelines.html).

