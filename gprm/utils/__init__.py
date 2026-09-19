"""Utility submodules for gprm.

Submodules are imported lazily (PEP 562) on first attribute access, so that importing one of
them does not import all fourteen. ``utils.spatial`` and friends behave exactly as before;
they are simply loaded at the point of use.
"""
from importlib import import_module as _import_module

_SUBMODULES = (
    'sphere',
    'create_gpml',
    'inpaint',
    'paleogeography',
    'pmag',
    'platetree',
    'geometry',
    'raster',
    'vector',
    'fileio',
    'spatial',
    'rotation',
    'proximity',
    'molchan',
)

__all__ = list(_SUBMODULES)


def __getattr__(name):
    if name in _SUBMODULES:
        module = _import_module('.' + name, __name__)
        globals()[name] = module     # cache, so this runs once per submodule
        return module
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))


def __dir__():
    return sorted(set(list(globals()) + __all__))
