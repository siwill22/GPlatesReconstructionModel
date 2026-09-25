"""The BRIDGE run list gprm pins checksums against. Offline: no download."""
from gprm.datasets._valdes2021_runs import RUNS


def test_109_runs_oldest_first():
    assert len(RUNS) == 109
    ages = [age for _, age, _ in RUNS]
    assert ages == sorted(ages, reverse=True)
    assert ages[0] == 541 and ages[-1] == 0


def test_every_run_has_its_own_checksum():
    assert len({sha for _, _, sha in RUNS}) == 109
    assert all(len(sha) == 64 for _, _, sha in RUNS)


def test_cached_file_names_do_not_collide_on_a_case_insensitive_filesystem():
    """Run codes are unique only when case matters (teXPb, teXpb and texpb are three runs).
    Geode cached them under bare run codes on macOS and 45 of 109 runs were silently
    overwritten by another run's data. fetch_Valdes2021 prefixes each file with its index."""
    run_codes = [run for run, _, _ in RUNS]
    assert len({r.lower() for r in run_codes}) < len(run_codes)  # the hazard is real
    names = ['{:03d}_{:s}a.pdclann.nc'.format(i, run) for i, run in enumerate(run_codes)]
    assert len({n.lower() for n in names}) == len(names)


def test_ocean_surface_checksums_cover_the_same_runs():
    from gprm.datasets._valdes2021_runs import RUNS, OCEAN_SURFACE_SHA256

    assert set(OCEAN_SURFACE_SHA256) == {run for run, _, _ in RUNS}
    assert len(set(OCEAN_SURFACE_SHA256.values())) == len(RUNS)

