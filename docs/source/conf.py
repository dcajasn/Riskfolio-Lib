# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import os
import sys
import sphinx_rtd_theme

import riskfolio as rp

sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(1, os.path.abspath('../../riskfolio'))
sys.path.insert(2, os.path.abspath('../../riskfolio/src'))
sys.path.insert(3, os.path.abspath('../../riskfolio/external'))

# -- Project information -----------------------------------------------------

project = 'Riskfolio-Lib'
copyright = '2020-2024, Dany Cajas'
author = 'Dany Cajas'


# The short X.Y version
release = '.'.join(rp.__version__.split('.')[:2])

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
# #    'numpydoc',
#     'sphinx_rtd_theme',
     'sphinx.ext.autodoc',
#     'sphinx.ext.doctest',
     'sphinx.ext.autosummary',
#     'sphinx_autopackagesummary',
     'sphinx.ext.intersphinx',
     'sphinx.ext.todo',
#     'sphinx.ext.coverage',
     'sphinx.ext.mathjax',
#     'sphinx.ext.ifconfig',
     'sphinx.ext.viewcode',
# #    'sphinx.ext.githubpages',
     'sphinx.ext.napoleon',
     'sphinxcontrib.bibtex',
     'sphinxemoji.sphinxemoji',
     'sphinx_immaterial',
]

autodoc_mock_imports = ["riskfolio.external.functions"]
autodoc_member_order = 'bysource'
autosummary_generate = True
keep_warnings = False

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
#source_suffix = ['.rst', '.bib']
source_suffix = '.rst'
bibtex_bibfiles = ['biblio.bib']


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_title = f"{project} {release}"
html_theme = 'sphinx_immaterial'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_search_language = 'en'
html_theme_options = {
    "palette": {"scheme": "default"},
    "icon": {"repo": "fontawesome/brands/github"},
    "site_url": "https://riskfolio-lib.readthedocs.io/en/latest",
    "repo_url": "https://github.com/dcajasn/Riskfolio-Lib",
    "repo_name": "Riskfolio-Lib",
    "globaltoc_collapse": True,
    "toc_title": "Contents",
    "toc_title_is_page_title": True,
    "social": [
        {
            "icon": "fontawesome/brands/github",
            "link": "https://github.com/dcajasn/Riskfolio-Lib",
            "name": "Source on github.com",
        },
        {
            "icon": "fontawesome/brands/python",
            "link": "https://pypi.org/project/Riskfolio-Lib/",
        },
    ],}
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_favicon = "_static/Riskfolio.ico"
html_logo = "_static/Riskfolio.png"

# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
# intersphinx_mapping = {'https://docs.python.org/': None}
python_version = '.'.join(map(str, sys.version_info[0:2]))
intersphinx_mapping = {
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
    'python': ('https://docs.python.org/' + python_version, None),
    'matplotlib': ('https://matplotlib.org', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'sklearn': ('https://scikit-learn.org/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/reference/', None),
    'arch': ('https://bashtage.github.io/arch/', None),
    'xlsxwriter': ('https://xlsxwriter.readthedocs.io', None),
    'networkx': ('https://networkx.org/documentation/stable/', None),
    'astropy': ('https://docs.astropy.org/en/stable/', None),
    'pybind11': ('https://pybind11.readthedocs.io/en/stable/',None),
}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True