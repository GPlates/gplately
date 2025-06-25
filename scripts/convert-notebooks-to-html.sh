#!/bin/bash

# run this script in the root directory of GPlately repository

jupyter-nbconvert --to=html --output-dir=notebook-html/ Notebooks/*.ipynb

jupyter-nbconvert --to=html --output-dir=notebook-html/ Notebooks/Examples/*.ipynb