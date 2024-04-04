# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

#import sys, os
# sys.path.insert(0, 'slcode/')

project = 'SL_C0de'
copyright = '2024, HENRY Adrien'
author = 'HENRY Adrien'
release = '1.0.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinxcontrib.bibtex',
]

bibtex_encoding = 'latin'
bibtex_bibfiles = ['reference.bib']
bibtex_default_style = 'plain'
bibtex_reference_style= 'author_year'
templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# import sphinx_rtd_theme

html_theme = 'sphinx_rtd_theme'

#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme_options = {
    # 'analytics_id': 'G-XXXXXXXXXX',  #  Provided by Google in your dashboard
    # 'analytics_anonymize_ip': False,
    # 'logo_only': False,
    # 'display_version': True,
    # 'prev_next_buttons_location': 'bottom',
    # 'style_external_links': False,
    # 'vcs_pageview_mode': '',
    # 'style_nav_header_background': 'white',
    # # Toc options
    'collapse_navigation': True,
    # 'sticky_navigation': True,
    'navigation_depth': -1
    # 'includehidden': True,
    # 'titles_only': False
}

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
    'figure_align': 'tbh',
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    'preamble': r'''\renewcommand{\hyperref}[2][]{(#2 p.\pageref{#1})}'''
    #'preamble': r'''\renewcommand{\hyperref}[2][]{}'''

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
    (master_doc, 'SLC0de.tex', 'SL$_{C0de}$ Documentation',
     'Adrien Henry', 'manual'),
]




