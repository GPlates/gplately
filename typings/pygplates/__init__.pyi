# pygplates.pyi (empty or near-empty file)
from typing import Any

def __getattr__(name: str) -> Any: ...

# Add a stub-of-nothing to explicitly tell Pylance "this module has no types, stop complaining":

# Placed in your stub path, this single-line trick tells the type checker "anything accessed on this module is Any"
# — silences every attribute/import warning in one shot, without writing any real signatures. This is actually a common,
# legitimate pattern for untyped C extensions, not just a hack.
