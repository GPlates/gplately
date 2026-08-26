#!/usr/bin/env bash
#
# Execute one or more Jupyter notebooks and convert each to HTML.
#
# Usage:
#   ./scripts/run_notebooks.sh <notebook> [<notebook> ...]
#
# Each <notebook> is the path of a notebook under Notebooks/, without the
# .ipynb extension, e.g.:
#
#   ./scripts/run_notebooks.sh 01-GettingStarted 10-SeafloorGrids
#   ./scripts/run_notebooks.sh Examples/hello_world
#
# Run this from the repository root (the directory containing Notebooks/).
# It assumes `jupyter nbconvert` is already available on PATH — activate
# your conda/micromamba environment first if running locally, e.g.:
#
#   micromamba activate <env-name>
#   ./scripts/run_notebooks.sh 10-SeafloorGrids
#
# Environment variables:
#   NOTEBOOKS_SRC_DIR   Directory containing the source .ipynb files
#                        (default: Notebooks)
#   NOTEBOOK_RUN_DIR     Scratch directory for executed notebooks
#                        (default: /tmp/notebook-run)
#   NOTEBOOK_HTML_DIR    Output directory for the generated HTML
#                        (default: notebook-html)

set -euo pipefail

if [ "$#" -eq 0 ]; then
  echo "Usage: $0 <notebook> [<notebook> ...]" >&2
  echo "Example: $0 01-GettingStarted 10-SeafloorGrids" >&2
  exit 1
fi

src_dir="${NOTEBOOKS_SRC_DIR:-Notebooks}"
run_dir="${NOTEBOOK_RUN_DIR:-/tmp/notebook-run}"
html_dir="${NOTEBOOK_HTML_DIR:-notebook-html}"

export DISABLE_GPLATELY_DEV_WARNING=true

if ! command -v jupyter >/dev/null 2>&1; then
  echo "Error: 'jupyter' not found on PATH. Activate your environment first." >&2
  exit 1
fi

if [ ! -d "${src_dir}" ]; then
  echo "Error: source notebook directory '${src_dir}' not found." >&2
  echo "Run this script from the repository root, or set NOTEBOOKS_SRC_DIR." >&2
  exit 1
fi

for notebook in "$@"; do
  source_path="${src_dir}/${notebook}.ipynb"
  if [ ! -f "${source_path}" ]; then
    echo "Error: notebook not found: ${source_path}" >&2
    exit 1
  fi

  notebook_file="$(basename "${notebook}").ipynb"

  echo "==> Executing ${source_path}"
  mkdir -p "${run_dir}/${notebook}"
  jupyter nbconvert \
    --to notebook \
    --execute \
    --ExecutePreprocessor.force_raise_errors=True \
    --output "${notebook_file}" \
    --output-dir "${run_dir}/${notebook}" \
    "${source_path}"

  echo "==> Converting to HTML"
  mkdir -p "${html_dir}/$(dirname "${notebook}")"
  jupyter nbconvert \
    --to html \
    --ExecutePreprocessor.force_raise_errors=True \
    --output "${notebook}.html" \
    --output-dir "${html_dir}" \
    "${run_dir}/${notebook}/${notebook_file}"

  echo "==> Done: ${html_dir}/${notebook}.html"
done
