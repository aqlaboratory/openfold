import os
import glob
import importlib as importlib

_files = glob.glob(os.path.join(os.path.dirname(__file__), "*.py"))
__all__ = [
    os.path.basename(f)[:-3]
    for f in _files
    if os.path.isfile(f) and not f.endswith("__init__.py")
]
_modules = [(m, importlib.import_module("." + m, __name__)) for m in __all__]
for _m in _modules:
    globals()[_m[0]] = _m[1]

# Avoid needlessly cluttering the global namespace
del _files, _m, _modules
