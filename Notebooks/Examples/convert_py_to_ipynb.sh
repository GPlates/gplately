#!/usr/bin/env bash

set -euo pipefail

if ! command -v jupytext >/dev/null 2>&1; then
	echo "jupytext is required but not found in PATH." >&2
	exit 1
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
	echo "This script must be run inside a git repository." >&2
	exit 1
fi

files=(
	hello_world
	PNG_reconstruction_copper_deposits
	icosahedron_mesh
	reconstruct_files
	introducing_plate_model_manager
	save_reconstructed_data
	plot_map_with_cartopy
	use_auxiliary_functions
	plot_map_with_pygmt
	use_your_own_plate_model
)

for base in "${files[@]}"; do
	py_file="${base}.py"
	ipynb_file="${base}.ipynb"

	if [ ! -f "${py_file}" ]; then
		echo "Skipping ${py_file}: file not found"
		continue
	fi

	# Skip tracked files that are unchanged vs HEAD.
	if git ls-files --error-unmatch "${py_file}" >/dev/null 2>&1 && \
		git diff --quiet HEAD -- "${py_file}"; then
		echo "Skipping ${py_file}: no git changes"
		continue
	fi

	echo "Converting ${py_file} -> ${ipynb_file}"
	jupytext --to notebook "${py_file}" -o "${ipynb_file}"
done