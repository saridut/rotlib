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
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'rotlib'
copyright = '2026, Sarit Dutta'
author = 'Sarit Dutta'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

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
    "private-members": False,
}
toc_object_entries_show_parents = "hide"

intersphinx_mapping = {'numpy': ('https://numpy.org/doc/stable/', None)}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output
 
# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
html_theme = 'pydata_sphinx_theme'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True

html_theme_options = {
    #[...]
    "primary_sidebar_end": ["indices.html", "sidebar-ethical-ads.html"]
    #[...]
}
#html_theme_options = {
#    "use_edit_page_button": True,
#}
