"""The per-model reference frames in gprm/datasets/_frames.py, checked without a network.

Whether each entry matches the rotation files a fetcher actually loads is checked in
test_reconstruction_loaders.py, which needs the data.
"""
import re
from pathlib import Path

import pytest

from gprm.datasets import _frames

DATASETS_DIR = Path(__file__).parent.parent / 'gprm' / 'datasets'


@pytest.mark.parametrize('name', _frames.models())
def test_every_entry_is_well_formed(name):
    frames = _frames.reference_frames(name)

    assert frames, 'a registered model with no frames should not be registered'
    plate_ids = [frame['plate_id'] for frame in frames]
    assert plate_ids == sorted(set(plate_ids)), 'plate ids must be unique and in order'
    for frame in frames:
        assert set(frame) == {'plate_id', 'reference', 'description', 'valid', 'source', 'note'}
        assert isinstance(frame['plate_id'], int) and frame['plate_id'] >= 0
        assert frame['reference'] in _frames.REFERENCES
        assert isinstance(frame['description'], str) and frame['description']
        assert isinstance(frame['source'], str) and frame['source']
        assert frame['note'] is None or isinstance(frame['note'], str)
        if frame['valid'] is not None:
            young, old = frame['valid']
            assert 0 <= young < old


@pytest.mark.parametrize('name', _frames.models())
def test_plate_0_is_always_documented(name):
    """Plate 0 is the anchor everyone uses by default, so every registered model says what it is."""
    assert _frames.reference_frames(name)[0]['plate_id'] == 0


def test_there_are_only_two_references():
    """Frames are classified by what they are fixed to, not how they were built: no-net-rotation
    and Pacific hotspot frames are kinds of mantle frame, not categories of their own."""
    assert _frames.REFERENCES == ('mantle', 'spin axis')


def test_frames_with_two_references_are_the_ones_the_sources_describe():
    """The models known to carry both a mantle and a spin-axis frame, and on which plates."""
    both = {}
    for name in _frames.models():
        by_reference = {}
        for frame in _frames.reference_frames(name):
            by_reference.setdefault(frame['reference'], []).append(frame['plate_id'])
        if len(by_reference) == 2:
            both[name] = by_reference['spin axis']

    assert both == {'TorsvikCocks2017': [1], 'DomeierTorsvik2014': [1],
                    'Muller2025:Opt': [5], 'Muller2025:NNR': [5]}


def test_reference_frames_returns_a_copy():
    frames = _frames.reference_frames('TorsvikCocks2017')
    frames[0]['reference'] = 'spin axis'
    frames.append({})
    assert _frames.reference_frames('TorsvikCocks2017')[0]['reference'] == 'mantle'
    assert len(_frames.reference_frames('TorsvikCocks2017')) == 5


def test_unknown_model_is_rejected_with_the_known_names():
    with pytest.raises(KeyError, match='TorsvikCocks2017'):
        _frames.reference_frames('NoSuchModel')


def test_reference_frames_is_exported_from_gprm_datasets():
    from gprm.datasets import reference_frames
    assert reference_frames is _frames.reference_frames


def test_every_fetcher_records_its_frames_under_a_registered_name():
    """Read the fetchers' source rather than running them, which would need the network. Every
    fetch_ function must hand its model to _with_frames, and every registry entry must be
    reachable."""
    source = (DATASETS_DIR / 'Reconstructions.py').read_text()
    functions = re.split(r'\ndef ', source)

    used = set()
    for function in functions:
        name = function.split('(', 1)[0]
        if not name.startswith('fetch_'):
            continue
        calls = re.findall(r'_with_frames\(reconstruction_model, (.+)\)', function)
        assert calls, '{} does not record its reference frames'.format(name)
        for call in calls:
            literals = re.findall(r"'([^']+)'", call)
            if call.endswith('+ model_case'):
                prefix = literals[0]
                used |= {key for key in _frames.models() if key.startswith(prefix)}
            else:
                used |= set(literals)

    assert used <= set(_frames.models()), 'unregistered names: {}'.format(used - set(_frames.models()))
    assert set(_frames.models()) <= used, 'unreachable entries: {}'.format(set(_frames.models()) - used)


# --- ReconstructionModel side ----------------------------------------------------------------

def test_a_new_model_has_no_frames_and_says_so(reconstruction_model):
    assert reconstruction_model.reference_frames == []
    assert 'Reference Frames (anchor_plate_id):\n   - not documented' in repr(reconstruction_model)


def test_add_reference_frame_records_and_prints_it(reconstruction_model):
    reconstruction_model.add_reference_frame(801, 'mantle', 'test hotspots', valid=(0, 80))
    reconstruction_model.add_reference_frame(0, 'spin axis', 'test palaeomag')

    frames = reconstruction_model.reference_frames
    assert [frame['plate_id'] for frame in frames] == [0, 801]
    assert frames[1] == dict(plate_id=801, reference='mantle', description='test hotspots',
                             valid=(0., 80.), source='user', note=None)

    text = repr(reconstruction_model)
    assert '   - 0: spin axis -- test palaeomag\n' in text
    assert '   - 801: mantle -- test hotspots, 0-80 Ma only\n' in text


@pytest.mark.parametrize('kwargs, message', [
    (dict(reference='hotspot'), 'reference must be one of'),
    (dict(reference='no-net-rotation'), 'reference must be one of'),
    (dict(reference='mantle', valid=(80, 0)), 'young < old'),
    (dict(reference='mantle', valid=(-5, 10)), 'young < old'),
])
def test_add_reference_frame_rejects_bad_input(reconstruction_model, kwargs, message):
    with pytest.raises(ValueError, match=message):
        reconstruction_model.add_reference_frame(0, description='x', **kwargs)


def test_add_reference_frame_refuses_a_second_frame_on_the_same_plate(reconstruction_model):
    reconstruction_model.add_reference_frame(0, 'mantle', 'first')
    with pytest.raises(ValueError, match='already recorded for plate 0'):
        reconstruction_model.add_reference_frame(0, 'spin axis', 'second')


def test_a_shallow_copy_does_not_share_the_frame_list(reconstruction_model):
    reconstruction_model.add_reference_frame(0, 'mantle', 'original')
    duplicate = reconstruction_model.copy()
    duplicate.add_reference_frame(801, 'spin axis', 'only on the copy')

    assert len(reconstruction_model.reference_frames) == 1
    assert len(duplicate.reference_frames) == 2
