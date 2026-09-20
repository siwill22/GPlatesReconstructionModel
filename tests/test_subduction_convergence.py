"""Tests for gprm.utils.spatial.subduction_convergence, the pygplates-1.0-native replacement
for ptt.subduction_convergence.subduction_convergence().

The replacement exists because the ptt algorithm walks each subduction segment's shared
sub-segments individually and needs exactly one subducting plate directly attached; it warns
and silently drops a segment otherwise (observed on the standard EarthByte 230-0 Ma AREPS
model at several reconstruction times: "Unable to find the subducting plate of the subducting
sub-segment ..."). Sampling pygplates 1.0's plate-boundary statistics instead finds the
subducting plate as long as one exists on that side of the boundary. These tests check the two
implementations still agree wherever ptt's algorithm succeeds, on a small synthetic model that
does not need any external data file.
"""
import numpy as np


def test_matches_ptt_on_a_synthetic_two_plate_model(two_plate_subduction_topology):
    """Cross-check against the algorithm being replaced, on a model simple enough that both
    can be trusted to walk the same one subduction zone the same way."""
    import ptt.subduction_convergence
    from gprm.utils.spatial import subduction_convergence

    rotation_model, topological_features = two_plate_subduction_topology
    threshold = np.radians(5.)

    old = ptt.subduction_convergence.subduction_convergence(
        rotation_model, topological_features, threshold, 0., 1.0, 0)
    new = subduction_convergence(rotation_model, topological_features, threshold, 0.,
                                 velocity_delta_time=1., anchor_plate_id=0)

    assert len(old) == len(new) == 12  # 60 degrees of trench at 5 degree spacing

    old = sorted(old, key=lambda row: row[1])   # sort by latitude
    new = sorted(new, key=lambda row: row[1])

    old_arr = np.array(old, dtype=float)
    new_arr = np.array(new, dtype=float)

    np.testing.assert_allclose(old_arr[:, 0], new_arr[:, 0], atol=1e-6)   # lon
    np.testing.assert_allclose(old_arr[:, 1], new_arr[:, 1], atol=1e-6)   # lat
    np.testing.assert_allclose(old_arr[:, 2], new_arr[:, 2], atol=0.05)   # conv_rate, cm/yr
    np.testing.assert_allclose(old_arr[:, 3], new_arr[:, 3], atol=1.0)    # conv_obliq, degrees
    np.testing.assert_allclose(old_arr[:, 4], new_arr[:, 4], atol=0.05)   # migr_rate, cm/yr
    np.testing.assert_allclose(old_arr[:, 5], new_arr[:, 5], atol=1.0)    # migr_obliq, degrees
    np.testing.assert_allclose(old_arr[:, 6], new_arr[:, 6], atol=0.5)    # arc_length, degrees
    np.testing.assert_allclose(old_arr[:, 7], new_arr[:, 7], atol=1.0)    # arc_azimuth, degrees
    np.testing.assert_array_equal(old_arr[:, 8], new_arr[:, 8])          # subducting_plate
    np.testing.assert_array_equal(old_arr[:, 9], new_arr[:, 9])          # overriding_plate ("trench_plate")


def test_subducting_and_overriding_plates_are_not_swapped(two_plate_subduction_topology):
    """Regression for a sign error: subducting_plate_id was being read from the overriding
    plate's own side (stat.left_plate when overriding was on the left, rather than
    stat.right_plate), so it came out equal to overriding_plate_id whenever the trench
    feature's own reconstruction plate id matched the overriding plate."""
    from gprm.utils.spatial import subduction_convergence

    rotation_model, topological_features = two_plate_subduction_topology
    result = subduction_convergence(rotation_model, topological_features, np.radians(5.), 0.)

    assert len(result) > 0
    for row in result:
        subducting_plate_id, overriding_plate_id = row[8], row[9]
        assert subducting_plate_id == 100     # west, held fixed -- the subducting plate
        assert overriding_plate_id == 200     # east, rotating -- the overriding plate
        assert subducting_plate_id != overriding_plate_id


def test_subduction_convergence_class_still_builds_a_dataframe(two_plate_subduction_topology):
    """Wiring check for GPlatesReconstructionModel.SubductionConvergence, which now calls
    gprm.utils.spatial.subduction_convergence rather than ptt.subduction_convergence directly."""
    from gprm import ReconstructionModel, SubductionConvergence

    rotation_model, topological_features = two_plate_subduction_topology

    model = ReconstructionModel('TwoPlateTest')
    model.rotation_model = rotation_model
    model.dynamic_polygons = topological_features

    result = SubductionConvergence(model, reconstruction_times=0., threshold_sampling_distance_radians=np.radians(5.))

    assert len(result.df) == 12
    assert set(result.df.columns) == {'lon', 'lat', 'conv_rate', 'conv_obliq', 'migr_rate',
                                      'migr_obliq', 'arc_length', 'arc_azimuth',
                                      'subducting_plate', 'overriding_plate', 'time'}
    assert (result.df.subducting_plate == 100).all()
    assert (result.df.overriding_plate == 200).all()
