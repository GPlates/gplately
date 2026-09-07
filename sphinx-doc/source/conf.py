# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import gplately

version = gplately.__version__
release = ".".join(version.split(".")[:2])

project = "gplately"
copyright = "2023-2026, The University of Sydney"
author = "EarthByte Group"
release = "2.0.0"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_favicon = "favicon.ico"

autosummary_generate = True

# the following autodoc_default_options is here because a bug of autosummary
# some links in the autosummary tables are not clickable without these autodoc_default_options settings
# this is a workaround for a autosummary bug. If the bug is fixed, the autodoc_default_options can be removed.
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "private-members": False,
    "show-inheritance": True,
}


def autodoc_skip_member(app, what, name, obj, skip, options):
    """
    app     - the Sphinx application object
    what    - the type of the parent: 'module', 'class', 'exception', 'function', 'method', 'attribute'
    name    - the name of the member (e.g. 'generate', '__init__')
    obj     - the member object itself
    skip    - bool: whether Sphinx's default logic decided to skip it
    options - the options given to the directive (autodoc_default_options etc.)

    Return True to skip (exclude) the member, False to include it,
    or None to let Sphinx's default behavior stand.
    """
    if name in {
        "CURRENT_LATITUDES_KEY",
        "CURRENT_LONGITUDES_KEY",
    }:
        return True
    return None


def setup(app):
    app.connect("autodoc-skip-member", autodoc_skip_member)
