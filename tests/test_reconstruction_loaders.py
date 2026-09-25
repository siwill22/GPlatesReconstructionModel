"""Correctness checks for the reconstruction-model fetchers. Marked ``network`` and excluded
from the default run, like test_dataset_loaders.py: these download published models, and a
failure means an upstream archive changed rather than that gprm's own logic is broken.
"""
import os

import pytest

pytestmark = pytest.mark.network

pygplates = pytest.importorskip('pygplates')

# Ages spanning every rotation file segment of the Cao et al (2024) model (1000-0 and
# 1800-1000 Ma), including the boundaries between them.
AGES = [0., 50., 100., 250., 410., 600., 850., 1000., 1001., 1200., 1500., 1800.]


def _rotation_difference_degrees(r1, r2):
    difference = r1 * r2.get_inverse()
    if difference.represents_identity_rotation():
        return 0.
    return abs(difference.get_lat_lon_euler_pole_and_angle_degrees()[2])


def _max_difference(model1, anchor1, model2, anchor2, plate_ids):
    """Largest angle (degrees) between the two models' total rotations, over plate_ids x AGES."""
    worst = 0.
    for age in AGES:
        for plate_id in plate_ids:
            worst = max(worst, _rotation_difference_degrees(
                model1.get_rotation(age, plate_id, anchor_plate_id=anchor1),
                model2.get_rotation(age, plate_id, anchor_plate_id=anchor2)))
    return worst


@pytest.fixture(scope='module')
def muller2025():
    from gprm.datasets.Reconstructions import fetch_Muller2025
    return fetch_Muller2025()


@pytest.fixture(scope='module')
def muller2025_nnr():
    from gprm.datasets.Reconstructions import fetch_Muller2025
    return fetch_Muller2025(NNR=True)


@pytest.fixture(scope='module')
def cao2024():
    from gprm.datasets.Reconstructions import fetch_Cao2024
    return fetch_Cao2024()


@pytest.fixture(scope='module')
def shared_plate_ids(muller2025, cao2024):
    """Every real plate in the Cao et al (2024) rotation files. Plate 5 is left out: in Cao2024
    it is an empty placeholder, in Muller2025 it carries the optimised frame."""
    plate_ids = (cao2024.known_plate_ids() & muller2025.known_plate_ids()) - {0, 5}
    assert len(plate_ids) > 500
    return sorted(plate_ids)


def test_muller2025_components_load(muller2025):
    assert muller2025.name == 'Muller++2025_Opt'
    assert [os.path.basename(f) for f in muller2025.rotation_files] == [
        'optimised_rotation_model_20240725.rot']
    assert len(muller2025.static_polygons) > 0
    assert len(muller2025.coastlines) > 0
    assert len(muller2025.continent_polygons) > 0
    assert len(muller2025.dynamic_polygon_files) == 8


def test_muller2025_anchor_5_is_the_cao2024_palaeomag_frame(muller2025, cao2024, shared_plate_ids):
    """The point of loading this rotation file rather than the authors' recommended pair:
    anchoring on 5 skips the optimised 005-000 rotation, leaving the original frame."""
    worst = _max_difference(muller2025.rotation_model, 5,
                            cao2024.rotation_model, 0, shared_plate_ids)
    assert worst < 1e-6


def test_muller2025_anchor_0_is_the_authors_recommended_mantle_frame(muller2025, shared_plate_ids):
    """The combined file and the recommended pair (with plate 5 merged away) are two
    arrangements of the same mantle frame, per the archive README."""
    from pooch import os_cache
    from gprm.datasets._remote_zip import retrieve_zip_member

    recommended = [retrieve_zip_member(
        url='https://zenodo.org/records/17142287/files/Cao_etal_2024_1.8_Ga_mantle_ref_frame.zip',
        member='Cao_etal_2024_1.8_Ga_mantle_ref_frame/optimisation/{}'.format(name),
        known_hash=known_hash,
        path=os.path.join(str(os_cache('gprm')), 'Muller2025', 'recommended'))
        for name, known_hash in [
            ('1000_0_rotfile_20240725.rot',
             'sha256:e7c3b7d201033d51ce75b2545bfab01f904c2db5bf01d0dccc7436004ec1b1c5'),
            ('1800_1000_rotfile_20240725.rot',
             'sha256:fea2ac9d62e31e5fcaf671e8feba01bc018201423ed109f007aa3a06b124014f')]]

    # Not exact: merging 005-000 into the relative poles perturbs a few of them, by at most
    # 0.001 degrees (~100 m) at 1001 Ma and 0.0002 degrees at 1500 Ma; every other age agrees
    # to within 1e-6. Differences between frames are tens of degrees, so 0.01 still discriminates.
    worst = _max_difference(muller2025.rotation_model, 0,
                            pygplates.RotationModel(recommended), 0, shared_plate_ids)
    assert worst < 0.01


def test_muller2025_mantle_frame_is_not_the_palaeomag_frame(muller2025):
    """Guards against the two anchors silently collapsing onto each other: they agree at
    present day and part company in the past."""
    rotations = muller2025.rotation_model
    assert _rotation_difference_degrees(rotations.get_rotation(0., 701, anchor_plate_id=0),
                                        rotations.get_rotation(0., 701, anchor_plate_id=5)) < 1e-6
    for age in [100., 500., 1000., 1500.]:
        assert _rotation_difference_degrees(
            rotations.get_rotation(age, 701, anchor_plate_id=0),
            rotations.get_rotation(age, 701, anchor_plate_id=5)) > 5.


def test_muller2025_nnr_keeps_the_palaeomag_frame_on_plate_5(muller2025_nnr, muller2025,
                                                             cao2024, shared_plate_ids):
    assert muller2025_nnr.name == 'Muller++2025_NNR'
    assert [os.path.basename(f) for f in muller2025_nnr.rotation_files] == [
        'no_net_rotation_model_20240725.rot']

    worst = _max_difference(muller2025_nnr.rotation_model, 5,
                            cao2024.rotation_model, 0, shared_plate_ids)
    assert worst < 1e-6

    # ...while its plate 0 is a different frame from the optimised mantle one
    assert _rotation_difference_degrees(
        muller2025_nnr.rotation_model.get_rotation(500., 701, anchor_plate_id=0),
        muller2025.rotation_model.get_rotation(500., 701, anchor_plate_id=0)) > 1.


def test_muller2025_only_caches_the_members_not_the_335mb_archive(muller2025):
    from pooch import os_cache

    cache = str(os_cache('gprm'))
    assert not [f for f in os.listdir(cache) if f.endswith('Cao_etal_2024_1.8_Ga_mantle_ref_frame.zip')]

    member_dir = os.path.join(cache, 'Muller2025')
    total = sum(os.path.getsize(os.path.join(member_dir, f))
                for f in os.listdir(member_dir) if os.path.isfile(os.path.join(member_dir, f)))
    assert total < 100e6


@pytest.mark.parametrize('age', [0., 300., 700., 1200., 1700.])
def test_muller2025_topologies_match_cao2024_in_the_palaeomag_frame(muller2025, cao2024, age):
    """The archive README says topologies are close to identical to Cao et al (2024), with a few
    subduction zones nudged offshore. So the same plates should resolve, tiling the sphere."""
    ours = muller2025.plate_snapshot(age, anchor_plate_id=5)
    theirs = cao2024.plate_snapshot(age, anchor_plate_id=0)

    assert sorted(ours.plate_ids) == sorted(theirs.plate_ids)
    sphere_area = 4 * 3.141592653589793 * pygplates.Earth.mean_radius_in_kms ** 2
    assert sum(ours.plate_areas) == pytest.approx(sphere_area, rel=0.01)


@pytest.mark.parametrize('age', [0., 700., 1700.])
def test_muller2025_topologies_resolve_in_the_mantle_frame(muller2025, age):
    """Changing frame rotates the whole globe rigidly, so plate count and total area are unchanged."""
    mantle = muller2025.plate_snapshot(age, anchor_plate_id=0)
    palaeomag = muller2025.plate_snapshot(age, anchor_plate_id=5)

    assert sorted(mantle.plate_ids) == sorted(palaeomag.plate_ids)
    assert sum(mantle.plate_areas) == pytest.approx(sum(palaeomag.plate_areas), rel=1e-6)


def test_muller2025_polygons_reconstruct_in_the_palaeomag_frame(muller2025):
    snapshot = muller2025.polygon_snapshot('static_polygons', 1000., anchor_plate_id=5)
    assert len(snapshot.reconstructed_polygons) > 0


# --- Reference frames, checked against the rotation files each fetcher loads ------------------

FETCHERS_WITH_FRAMES = [
    ('fetch_Cao2024', {}, 'Cao2024'),
    ('fetch_CaoToyRodinia', {'model_case': 'NNR'}, 'CaoToyRodinia:NNR'),
    ('fetch_CaoToyRodinia', {'model_case': 'OV'}, 'CaoToyRodinia:OV'),
    ('fetch_CaoToyRodinia', {'model_case': 'SSL'}, 'CaoToyRodinia:SSL'),
    ('fetch_Li2008', {}, 'Li2008'),
    ('fetch_Li2023', {'model_case': 'East'}, 'Li2023:East'),
    ('fetch_Li2023', {'model_case': 'West'}, 'Li2023:West'),
    ('fetch_DomeierTorsvik2014', {}, 'DomeierTorsvik2014'),
    ('fetch_Matthews2016', {}, 'Matthews2016'),
    ('fetch_Merdith2021', {}, 'Merdith2021'),
    ('fetch_Muller2022', {'NNR': False}, 'Muller2022:Opt'),
    ('fetch_Muller2022', {'NNR': True}, 'Muller2022:NNR'),
    ('fetch_Muller2025', {'NNR': False}, 'Muller2025:Opt'),
    ('fetch_Muller2025', {'NNR': True}, 'Muller2025:NNR'),
    ('fetch_Muller2016', {}, 'Muller2016'),
    ('fetch_Muller2019', {}, 'Muller2019'),
    ('fetch_Pehrsson2015', {}, 'Pehrsson2015'),
    ('fetch_Seton2012', {}, 'Seton2012'),
    ('fetch_TorsvikCocks2017', {}, 'TorsvikCocks2017'),
    ('fetch_Young2019', {}, 'Young2019'),
    ('fetch_Scotese', {}, 'Scotese2008'),
    ('fetch_Clennett', {'model_case': 'M2019'}, 'Clennett:M2019'),
    ('fetch_Clennett', {'model_case': 'S2013'}, 'Clennett:S2013'),
]


def _pole_times(model, plate_id):
    """Every time at which the model's rotation files give plate_id a pole (as moving plate)."""
    times = set()
    for rotation_file in model.rotation_files:
        for feature in pygplates.FeatureCollection(rotation_file):
            pole = feature.get_total_reconstruction_pole()
            if pole and pole[1] == plate_id:
                times |= {sample.get_time() for sample in pole[2].get_enabled_time_samples()}
    return times


@pytest.fixture(scope='module', params=FETCHERS_WITH_FRAMES,
                ids=[key for _, _, key in FETCHERS_WITH_FRAMES])
def fetched(request):
    from gprm.datasets import Reconstructions
    fetcher, kwargs, key = request.param
    return getattr(Reconstructions, fetcher)(**kwargs), key


def test_fetcher_attaches_its_registered_frames(fetched):
    from gprm.datasets import reference_frames
    model, key = fetched
    assert model.reference_frames == reference_frames(key)
    assert 'not documented' not in repr(model)


def test_every_frame_plate_is_in_the_rotation_files(fetched):
    model, _ = fetched
    for frame in model.reference_frames:
        assert frame['plate_id'] in model.known_plate_ids(), frame


def test_a_frame_span_ends_at_its_last_pole(fetched):
    """The spans were read off the rotation files by hand; hold them to the files."""
    model, _ = fetched
    for frame in model.reference_frames:
        if frame['valid'] is None:
            continue
        young, old = frame['valid']
        times = _pole_times(model, frame['plate_id'])
        assert young in times and old in times, (frame, sorted(times))


def test_no_two_frames_of_a_model_are_the_same_frame(fetched):
    """A listed frame that coincided with another would be an alias, not a frame: somewhere in
    their common span, some major plate must sit in a different place."""
    model, _ = fetched
    frames = model.reference_frames
    for i, first in enumerate(frames):
        for second in frames[i + 1:]:
            oldest = min((frame['valid'] or (0., 1e9))[1] for frame in (first, second))
            ages = [age for age in (50., 100., 300., 500.) if age < oldest]
            worst = 0.
            for age in ages:
                for plate_id in (701, 101, 901, 801, 501):
                    worst = max(worst, _rotation_difference_degrees(
                        model.rotation_model.get_rotation(age, plate_id, anchor_plate_id=first['plate_id']),
                        model.rotation_model.get_rotation(age, plate_id, anchor_plate_id=second['plate_id'])))
            assert worst > 0.5, (first['plate_id'], second['plate_id'])
