"""Check that the built API page contains rendered Python-domain objects."""
from __future__ import annotations

import argparse
from html.parser import HTMLParser
from pathlib import Path


class ApiPage(HTMLParser):
    def __init__(self):
        super().__init__()
        self.ids = set()
        self.text = []

    def handle_starttag(self, tag, attrs):
        attributes = dict(attrs)
        if "id" in attributes:
            self.ids.add(attributes["id"])

    def handle_data(self, data):
        self.text.append(data)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("html_dir", type=Path)
    args = parser.parse_args()
    page = ApiPage()
    page.feed((args.html_dir / "api.html").read_text(encoding="utf-8"))
    required = {
        "pase.associate",
        "pase.AssociationResult",
        "pase.AssociationResult.summary",
        "pase.AssociationResult.plot_forest",
        "pase.datasets.load_gse234933_demo",
    }
    missing = required - page.ids
    if missing:
        raise SystemExit("Missing rendered API objects: " + ", ".join(sorted(missing)))
    text = " ".join(page.text)
    if any(marker in text for marker in (".. py:function::", ".. py:class::", ":param adata:")):
        raise SystemExit("API page exposes unrendered reStructuredText directives.")
    print("API_HTML_OK: five public API objects rendered with addressable signatures")


if __name__ == "__main__":
    main()
