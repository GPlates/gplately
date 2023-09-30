import argparse
import os
import sys, importlib


def main():
    print("execute command ....")
    print(importlib.resources.files("gplately"))


if __name__ == "__main__":
    main()
