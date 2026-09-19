"""The import contract.

gprm's package __init__ used to import the whole stack eagerly, which cost ~21 s cold and
forced at least one downstream project to load a submodule by file path to get around it.
These tests pin the laziness down, because it is the kind of property that regresses silently
the moment someone adds a convenient top-level import.
"""
import subprocess
import sys

import pytest


PUBLIC_CLASSES = [
    'ReconstructionModel',
    'ReconstructedPolygonSnapshot',
    'PlateTree',
    'GPlatesRaster',
    'PlateSnapshot',
    'MotionPathFeature',
    'FlowlineFeature',
    'VelocityField',
    'SubductionConvergence',
    'PointDistributionOnSphere',
    'CrossSection',
]


def run_in_subprocess(code):
    """Run code in a clean interpreter and return stdout; import state must not leak between tests."""
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    return result.stdout.strip()


def test_importing_gprm_does_not_import_the_heavy_stack():
    """import gprm must not pull in pygplates, pygmt, ptt or matplotlib."""
    loaded = run_in_subprocess(
        'import sys, gprm; '
        "print(','.join(m for m in ('pygplates','pygmt','ptt','matplotlib','xrspatial','datashader') "
        'if m in sys.modules))')
    assert loaded == '', 'import gprm pulled in: {}'.format(loaded)


def test_dataset_modules_import_without_pygplates_or_pygmt():
    """The pooch-only dataset loaders must work with pygplates and pygmt absent entirely.

    This is the property a downstream caller needs in order to read a palaeogeography grid
    without installing the reconstruction stack.
    """
    code = '''
import sys
BLOCKED = {'pygplates', 'pygmt', 'ptt'}
class Blocker:
    def find_spec(self, name, path=None, target=None):
        if name.split('.')[0] in BLOCKED:
            raise ImportError(name)
        return None
sys.meta_path.insert(0, Blocker())

import gprm
import gprm.datasets.Paleogeography
import gprm.datasets.Strat
import gprm.datasets.Rocks
import gprm.datasets.Zircons
import gprm.datasets.Geology
print('ok')
'''
    assert run_in_subprocess(code) == 'ok'


@pytest.mark.parametrize('name', PUBLIC_CLASSES)
def test_public_class_is_reachable(name):
    """Every documented class must still resolve through the lazy __getattr__."""
    import gprm

    assert isinstance(getattr(gprm, name), type)


def test_public_classes_appear_in_dir():
    """dir(gprm) must advertise the public API even before anything has been accessed."""
    listed = run_in_subprocess(
        'import gprm; print(",".join(sorted(n for n in dir(gprm) if n[0].isupper())))')
    assert sorted(listed.split(',')) == sorted(PUBLIC_CLASSES)


def test_utils_submodules_load_on_attribute_access():
    import gprm.utils as utils

    assert hasattr(utils.spatial, 'plate_boundary_intersections')
    assert hasattr(utils.proximity, 'boundary_proximity')
    assert hasattr(utils.geometry, 'nearest_feature')


def test_unknown_attribute_still_raises_attribute_error():
    import gprm

    with pytest.raises(AttributeError):
        gprm.NoSuchThing


def test_cache_path_is_platform_correct():
    """Callers must be able to ask for the cache directory rather than hardcoding it."""
    from pooch import os_cache

    from gprm.datasets import cache_path

    assert str(cache_path()) == str(os_cache('gprm'))
    assert cache_path('Model', 'file.shp').name == 'file.shp'
