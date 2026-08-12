#!/usr/bin/env python3
"""
Fail if any given Jupyter notebook has code cells with saved outputs or a
non-null execution_count. Used to keep notebooks committed "clean" (no run
artifacts) in version control.

Usage:
    python scripts/check_notebook_outputs.py [notebook.ipynb ...]

With no arguments, it checks every *.ipynb file tracked by git in the repo.
"""

import json
import subprocess
import sys


def tracked_notebooks():
    result = subprocess.run(
        ["git", "ls-files", "*.ipynb"],
        capture_output=True,
        text=True,
        check=True,
    )
    return [line for line in result.stdout.splitlines() if line]


def notebook_has_output(path):
    with open(path, "r", encoding="utf-8") as f:
        nb = json.load(f)

    for cell in nb.get("cells", []):
        if cell.get("cell_type") != "code":
            continue
        if cell.get("outputs"):
            return True
        if cell.get("execution_count") is not None:
            return True
    return False


def main():
    paths = sys.argv[1:] or tracked_notebooks()

    if not paths:
        print("No .ipynb files found.")
        return 0

    dirty = []
    for path in paths:
        try:
            if notebook_has_output(path):
                dirty.append(path)
        except (json.JSONDecodeError, OSError) as exc:
            print(f"::error file={path}::Could not read notebook: {exc}")
            dirty.append(path)

    if dirty:
        print("The following notebook(s) contain output cells or execution counts:")
        for path in dirty:
            print(f"  - {path}")
        print(
            "\nClear outputs before committing, e.g.:\n"
            "  jupyter nbconvert --clear-output --inplace <notebook>.ipynb"
        )
        return 1

    print(f"Checked {len(paths)} notebook(s); no output cells found.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
