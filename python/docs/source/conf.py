# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

# Import the local package regardless of the directory used to run Sphinx.
srcdir = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(srcdir))

# -- Project information -----------------------------------------------------

project = "rotlib"
copyright = "2026, Sarit Dutta"
author = "Sarit Dutta"
version = "1.0.1"
release = "1.0.1"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "myst_parser",
    "numpydoc",
]

exclude_patterns = []

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

autodoc_member_order = "bysource"
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "inherited-members": True,
    "show-inheritance": True,
    "member-order": "bysource",
}
toc_object_entries_show_parents = "hide"

intersphinx_mapping = {"numpy": ("https://numpy.org/doc/stable/", None)}

# MyST syntax extensions.
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "substitution",
    "tasklist",
]

# Shared TeX macros for HTML math rendering.
macros = {
    "bs": [r"\boldsymbol{#1}", 1],
    "bm": [r"\mathbf{#1}", 1],
    "omag": [r"\lVert \boldsymbol{\upomega} \rVert", 0],
}
mathjax4_config = {
    "loader": {"load": ["[tex]/upgreek"]},
    "tex": {
        "packages": {"[+]": ["upgreek"]},
        "macros": macros,
    },
}

# Provide the same macros when generating LaTeX/PDF documentation.
latex_elements = {
    "preamble": r"""
\usepackage{upgreek}
\newcommand{\bs}[1]{\boldsymbol{#1}}
\newcommand{\bm}[1]{\mathbf{#1}}
\newcommand{\omag}{\lVert \boldsymbol{\upomega} \rVert}
""",
}

# -- Options for HTML output -------------------------------------------------
html_theme = "pydata_sphinx_theme"

# SmartyPants will be used to convert quotes and dashes to typographically
# correct entities.
smartquotes = True
html_context = {"default_mode": "dark"}
html_theme_options = {
    "primary_sidebar_end": ["indices.html"],
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/saridut/rotlib",
            "icon": "fa-brands fa-github",
            "type": "fontawesome",
        },
    ],
}
