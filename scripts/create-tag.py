#!/usr/bin/env python3
import subprocess
import sys

if len(sys.argv) != 2:
    print("Usage: ./create-tag.py 2.0.0")
    sys.exit(1)

new_version = sys.argv[1]

subprocess.call(["git", "tag", "v" + new_version])
subprocess.call(["git", "push", "origin", "v" + new_version])
