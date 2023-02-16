# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import sys, os
sys.path.insert(0, os.path.abspath('C:/Users/ahenry01/Desktop/Python_code/SL_C0de_lib_0_4_0/src/SL_C0de'))

project = 'SL_C0de'
copyright = '2023, HENRY Adrien'
author = 'HENRY Adrien'
release = '0.4.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'

html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]
#html_logo = 'img/logo2.png'
# html_theme_options = {
#     'logo_only': True,
#     'display_version': False,
# }

# html_theme = 'alabaster'
# html_static_path = ['_static']

# -- Options for LaTeX output ---------------------------------------------

master_doc = 'index'

latex_elements = {
    'fncychap': '\\usepackage[Sonny]{fncychap}',
    'figure_align': 'tbh'
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}
latex_engine = 'pdflatex'
#latex_logo = 'img/logo_latex.pdf'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    (master_doc, 'badlands.tex', 'badlands Documentation',
     'Tristan Salles', 'manual'),
]


