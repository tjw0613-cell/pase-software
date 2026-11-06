"""Prepare a recipe from an existing PyPI source artifact; does not publish it."""
from __future__ import annotations

import argparse
from email.parser import BytesParser
import hashlib
import io
import json
from pathlib import Path
import re
import tarfile
from urllib.parse import urlparse
from urllib.request import urlopen

VERSION = "1.1.0"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-url", required=True, help="Public PyPI sdist URL")
    parser.add_argument("--maintainer", required=True, help="Consenting Bioconda recipe maintainer's GitHub username")
    parser.add_argument("--homepage", required=True, help="Actual GitHub software repository URL")
    parser.add_argument("--output", type=Path, default=Path("bioconda-recipe"))
    args = parser.parse_args()
    source = urlparse(args.source_url)
    home = urlparse(args.homepage)
    if source.scheme != "https" or source.hostname != "files.pythonhosted.org":
        parser.error("Use the exact HTTPS sdist URL hosted on files.pythonhosted.org.")
    if home.scheme != "https" or home.hostname != "github.com" or len(home.path.strip("/").split("/")) != 2:
        parser.error("--homepage must be the actual HTTPS GitHub owner/repository URL.")
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9-]{0,38}", args.maintainer):
        parser.error("Invalid GitHub maintainer name.")
    with urlopen(args.source_url, timeout=60) as response:
        payload = response.read(20 * 1024 * 1024 + 1)
    if len(payload) > 20 * 1024 * 1024:
        raise ValueError("Unexpected source archive size (>20 MB).")
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r:gz") as archive:
        matches = [m for m in archive.getmembers() if m.name.count("/") == 1 and m.name.endswith("/PKG-INFO")]
        if len(matches) != 1:
            raise ValueError("Expected one top-level PKG-INFO in a source distribution.")
        metadata = BytesParser().parsebytes(archive.extractfile(matches[0]).read())
    if metadata["Name"].replace("_", "-").lower() != "pase-framework" or metadata["Version"] != VERSION:
        raise ValueError("Source artifact has an unexpected package name/version.")
    sha = hashlib.sha256(payload).hexdigest()
    q = json.dumps
    recipe = f'''package:
  name: pase-framework
  version: {q(VERSION)}
source:
  url: {q(args.source_url)}
  sha256: {sha}
build:
  noarch: python
  number: 0
  script: "{{{{ PYTHON }}}} -m pip install . --no-deps --no-build-isolation -vv"
requirements:
  host:
    - python >=3.10
    - pip
    - setuptools >=77
    - wheel
  run:
    - python >=3.10
    - numpy >=1.24
    - scipy >=1.10
    - pandas >=2.0
    - anndata >=0.10
    - matplotlib-base >=3.7
test:
  imports:
    - pase
  commands:
    - python -c "import pase; a = pase.datasets.load_gse234933_demo(); assert a.n_obs > 0"
about:
  home: {q(args.homepage)}
  license: MIT
  license_file:
    - LICENSE
    - src/pase/data/DATA_NOTICE.md
  summary: "Patient-level association estimates for multi-patient single-cell studies"
extra:
  recipe-maintainers:
    - {args.maintainer}
'''
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "meta.yaml").write_text(recipe, encoding="utf-8")
    (args.output / "source-check.json").write_text(json.dumps({
        "source_url": args.source_url, "sha256": sha, "bytes": len(payload),
        "name": metadata["Name"], "version": metadata["Version"],
        "status": "recipe_generated; not built or submitted to Bioconda"
    }, indent=2), encoding="utf-8")
    print(args.output / "meta.yaml")


if __name__ == "__main__":
    main()

