#!/bin/bash

# PLEASE READ!!!
# this script was only tested on macOS. it may not work on Windows.
# run this script in the root directory of gplately repository.
# you need to have micromamba installed and an env called gplately.
# if something strange happened, make sure the executables of these commands were from the micromamba env, such as pip-compile.
# ask Michael Chin if you need help.

eval "$(micromamba shell hook --shell bash)"
micromamba activate gplately

pip-compile pyproject.toml
pip3 install .
pip3 install -U sphinx sphinx_rtd_theme
sphinx-apidoc -f -o sphinx-doc/source gplately/
cd sphinx-doc
make html