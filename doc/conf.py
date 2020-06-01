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
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))
import sys
import os
import shutil
sys.path.insert(1,'/Users/simon/GPlatesBuilds/pygplates_rev28_python36_MacOS64/')

# Based on https://github.com/spatialaudio/nbsphinx/issues/170,
# with some modifications to remove checkpoint files
print("Copy example notebooks into docs/_examples")

def all_but_ipynb(dir, contents):
    result = []
    for c in contents:
        if ("checkpoint" in c):
            result += [c]
        elif os.path.isfile(os.path.join(dir,c)) and (not c.endswith(".ipynb")):
            result += [c]
    return result

shutil.rmtree(os.path.join(os.path.abspath('..'), "doc/_examples"), ignore_errors=True)
shutil.copytree(os.path.join(os.path.abspath('..'), "test_notebooks"),
                os.path.join(os.path.abspath('..'), "doc/_examples"),
                ignore=all_but_ipynb)

# -- Project information -----------------------------------------------------

project = 'GPlatesReconstructionModel'
copyright = '2020, Simon Williams'
author = 'Simon Williams'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'nbsphinx']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']