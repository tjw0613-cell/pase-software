from importlib.metadata import version as package_version

project = "PASE"
author = "Jinwei Tang, Yunfan Lin, Linrui Cai, Yali Wen, Guangyang Chen, Li Huang"
copyright = "2026, PASE authors"
release = package_version("pase-framework")
version = release
extensions = ["myst_parser", "sphinx.ext.autodoc", "sphinx.ext.napoleon", "sphinx.ext.mathjax"]
source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
master_doc = "index"
exclude_patterns = ["_build"]
html_theme = "sphinx_rtd_theme"
html_title = f"PASE {release}"
myst_enable_extensions = ["dollarmath", "colon_fence"]
autodoc_member_order = "bysource"
napoleon_numpy_docstring = True

