#!/bin/bash
python -m pip install --user --upgrade setuptools wheel build --ignore-installed -v .
python -m build
twine upload dist/*
