# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html
import os

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'RETIMAT'
author = 'David Romero-Bascones'
version = '1.0.0'
release = version

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['sphinx.ext.viewcode', 'sphinxcontrib.matlab','sphinx.ext.autodoc']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

add_module_names = False # show only function name

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'classic'
html_static_path = ['_static']
html_sidebars = {'**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html']}
html_show_sphinx = False
html_title = f'{project} {version}'
html_domain_indices = False
html_use_index = False
html_show_sourcelink = False
html_show_copyright = False

# -- MATLAB configuration ----------------------------------------------------
primary_domain = 'mat'
parent_dir = os.path.abspath(os.path.join(os.getcwd(), os.pardir))
matlab_src_dir = os.path.join(parent_dir, 'src')

