#!/bin/bash

BASEDIR=$(dirname "$0")
echo "$BASEDIR"

rm -rf $BASEDIR/../dist

pip3 install -U build twine packaging
python3 -m build
twine check dist/*
twine upload --non-interactive -u __token__ -p $1 dist/*
