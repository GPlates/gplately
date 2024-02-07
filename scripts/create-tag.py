#!/usr/bin/env python3
import os
import subprocess
import sys

if len(sys.argv) != 2:
    print("Usage: ./create-tag.py VERSION(such as 1.1.0)")
    sys.exit(1)

new_version = sys.argv[1]

this_file_path = os.path.dirname(__file__)

# 1. update version number in gplately/__init__.py
# 2. commit change, create a tag and push changes and the new tag

lines = []
with open(f"{this_file_path}/../gplately/__init__.py", "r") as f:
    for line in f:
        if line.startswith("__version__"):
            lines.append(f'__version__ = "{new_version}"\n')
        else:
            lines.append(line)

with open(f"{this_file_path}/../gplately/__init__.py", "w") as of:
    of.writelines(lines)

subprocess.call(
    ["git", "add", f"{this_file_path}/../gplately/__init__.py"]
)
subprocess.call(["git", "commit", "-m", f"update version to {new_version}"])
subprocess.call(["git", "push"])

subprocess.call(["git", "tag", "v" + new_version])
subprocess.call(["git", "push", "origin", "v" + new_version])

