"""Loaders for published palaeo-geoscience datasets and reconstruction models.

Submodules are imported lazily (PEP 562), so that reaching for one dataset does not import
the others. This matters because the modules have very different dependencies: Paleogeography,
Geology, Rocks, Strat and Zircons need only pooch and pandas/xarray, while Reconstructions
needs pygplates.
"""
from importlib import import_module as _import_module

_SUBMODULES = (
    'Seafloor',
    'Rocks',
    'Reconstructions',
    'Paleogeography',
    'Strat',
    'Geology',
    'Zircons',
)

__all__ = list(_SUBMODULES) + ['cache_path', 'DatasetFetchError', 'age_description']


def cache_path(*parts):
    """Return the directory gprm downloads datasets into, optionally joined with sub-paths.

    The location is platform-dependent (``~/Library/Caches/gprm`` on macOS,
    ``~/.cache/gprm`` on Linux, ``%LOCALAPPDATA%\\gprm`` on Windows), so it should always be
    asked for rather than written out by hand.

    Use this to reach a file that a fetcher downloads but does not itself return, for example
    an auxiliary layer inside a bundle::

        from gprm.datasets import cache_path
        hotspots = cache_path('TorsvikCocks2017', 'Hotspot_Surface_Motion_PD2012.shp')

    :param parts: Optional path components to append to the cache directory.
    :returns: pathlib.Path. Note that it is only populated once the relevant fetch_ function
        has been run at least once; this does not download anything.
    """
    from pathlib import Path
    from pooch import os_cache

    return Path(os_cache('gprm')).joinpath(*parts)


def __getattr__(name):
    if name in _SUBMODULES:
        module = _import_module('.' + name, __name__)
        globals()[name] = module     # cache, so this runs once per submodule
        return module
    if name == 'age_description':
        # What "age" means in each dataset; see _ages.py
        from ._ages import age_description
        globals()[name] = age_description
        return age_description
    if name == 'DatasetFetchError':
        from ._fetch import DatasetFetchError
        globals()[name] = DatasetFetchError
        return DatasetFetchError
    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))


def __dir__():
    return sorted(set(list(globals()) + __all__))
