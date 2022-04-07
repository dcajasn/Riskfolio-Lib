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

sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(1, os.path.abspath('../../riskfolio'))


# -- Project information -----------------------------------------------------

project = 'Riskfolio-Lib'
copyright = '2020-2022, Dany Cajas'
author = 'Dany Cajas'

# The short X.Y version
version = 'latest'
# The full version, including alpha/beta/rc tags
release = '3.1.0'


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
]

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
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_search_language = 'en'

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
    'networkx': ('https://networkx.org', None),
    'astropy': ('https://www.astropy.org', None),
}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True