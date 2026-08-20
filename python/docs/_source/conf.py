# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path
import numpydoc

srcdir = Path('../../rotlib').resolve()
sys.path.insert(0, str(srcdir))


# -- Project information -----------------------------------------------------

project = 'rotlib'
copyright = '2026, Sarit Dutta'
author = 'Sarit Dutta'
release = '1.0.0'

# -- General configuration ---------------------------------------------------

extensions = [
        'sphinx.ext.duration',
        'sphinx.ext.autodoc',
        'sphinx.ext.autosummary',
        'sphinx.ext.viewcode',
        'sphinx.ext.intersphinx',
        'numpydoc',
        ]

exclude_patterns = []

autodoc_member_order = 'bysource'
autosummary_generate = True
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "inherited-members": True,
    "show-inheritance": True,
    "member-order": 'bysource',
}
toc_object_entries_show_parents = "hide"

intersphinx_mapping = {'numpy': ('https://numpy.org/doc/stable/', None)}

# -- Options for HTML output -------------------------------------------------
html_theme = 'pydata_sphinx_theme'

# SmartyPants will be used to convert quotes and dashes to typographically
# correct entities.
html_use_smartypants = True
html_theme_options = {
    "primary_sidebar_end": ["indices.html"]
}
#html_theme_options = {
#    "use_edit_page_button": True,
#}
