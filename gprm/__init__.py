"""Tools for building and working with GPlates plate tectonic reconstruction models.

The public classes are imported lazily (PEP 562), so that importing a part of gprm does not
drag in the whole stack. In particular ``import gprm.datasets.Paleogeography`` pulls in only
pooch and xarray -- it does not need pygplates, pygmt or ptt, which between them account for
most of the import cost.
"""
import logging as _logging
from importlib import import_module as _import_module

# Progress messages from the slower routines go through logging rather than print, so that
# importing and calling gprm produces no output of its own. A library that prints cannot be
# quietened by its caller, and the messages were unusable anyway inside a loop over hundreds
# of reconstruction times. To see them::
#
#     import logging
#     logging.basicConfig(level=logging.INFO)
#
# The NullHandler is the standard way for a library to avoid "No handlers could be found"
# while leaving all configuration to the application.
_logging.getLogger(__name__).addHandler(_logging.NullHandler())

# Public name -> the submodule that defines it.
_LAZY_ATTRIBUTES = {
    'ReconstructionModel': '.GPlatesReconstructionModel',
    'ReconstructedPolygonSnapshot': '.GPlatesReconstructionModel',
    'PlateTree': '.GPlatesReconstructionModel',
    'GPlatesRaster': '.GPlatesReconstructionModel',
    'PlateSnapshot': '.GPlatesReconstructionModel',
    'MotionPathFeature': '.GPlatesReconstructionModel',
    'FlowlineFeature': '.GPlatesReconstructionModel',
    'VelocityField': '.GPlatesReconstructionModel',
    'SubductionConvergence': '.GPlatesReconstructionModel',
    'PointDistributionOnSphere': '.GPlatesReconstructionModel',
    'CrossSection': '.GPlatesReconstructionModel',
}

_SUBPACKAGES = ('utils', 'datasets')

__all__ = sorted(_LAZY_ATTRIBUTES) + list(_SUBPACKAGES)


def __getattr__(name):
    if name in _LAZY_ATTRIBUTES:
        value = getattr(_import_module(_LAZY_ATTRIBUTES[name], __name__), name)
        globals()[name] = value      # cache, so this runs once per name
        return value

    if name in _SUBPACKAGES:
        module = _import_module('.' + name, __name__)
        globals()[name] = module
        return module

    raise AttributeError("module {!r} has no attribute {!r}".format(__name__, name))


def __dir__():
    return sorted(set(list(globals()) + __all__))
