## Developer Note

This folder includes local Pylance type stubs for `pygplates` in
[`typings/pygplates/__init__.pyi`](typings/pygplates/__init__.pyi).
These stubs exist to avoid false-positive attribute diagnostics from static
analysis on the compiled `pygplates` module.

Below is what Claude told me.

Add a stub-of-nothing to explicitly tell Pylance "this module has no types, stop complaining":

Placed in your stub path, this single-line trick tells the type checker "anything accessed on this module is Any" — silences every attribute/import warning in one shot, without writing any real signatures. This is actually a common, legitimate pattern for untyped C extensions, not just a hack.