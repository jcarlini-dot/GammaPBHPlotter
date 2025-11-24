# docs/source/conf.py
import os, sys
sys.path.insert(0, os.path.abspath('../../src'))  # make 'gammapbh' importable

project = 'GammaPBHPlotter'
author = 'John Carlini and Ilias Cholis'
copyright = '2025'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
]
autosummary_generate = True
autodoc_default_options = {
    'members': True,
    'undoc-members': True,
    'show-inheritance': True,
}
napoleon_google_docstring = True
napoleon_numpy_docstring = True

templates_path = ['_templates']
exclude_patterns = []
html_theme = 'furo'  # or 'alabaster'
