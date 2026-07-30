# Configuration file for the Sphinx documentation builder.

# -- Project information -----------------------------------------------------
project = 'Computing Nexus docs'
copyright = '2026, Shahid Maqbool'
author = 'Shahid Maqbool'
release = '1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    "myst_parser",              # Markdown support
    "sphinx.ext.napoleon",      # Google/NumPy style docstrings
    "sphinx_fortran_domain",    #  Fortran support
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

# Default highlighting language
highlight_language = 'cpp'