# Bundle the Boucot dataset in the repo rather than fetching it with pooch

Every other loader in `gprm/datasets/` retrieves its data over the network with `pooch`.
The Boucot palaeolithology compilation is bundled in `gprm/Data` instead, as a slimmed
GeoPackage (~1.8 MB), because we already hold and curate the file — it is not fetched from
a canonical upstream, so a pooch URL would mean first publishing our own derived product
somewhere and then depending on the network to read data we ship anyway.

## Consequences

`gprm/Data` becomes a real dependency of `gprm.datasets`, not just of the healpix meshes in
the core module. `setup.py` sets `include_package_data=True` with no `MANIFEST.in`,
`setup.cfg` or `pyproject.toml`, which means package data is silently excluded from built
wheels and sdists — this already breaks the healpix meshes for pip-installed users. A
`MANIFEST.in` is added as part of this change.

The 12.2 MB original shapefile set remains in git history even after removal, so the
repository does not shrink; only the working tree does.

The source shapefile's `PLATEID1`/`PLATEID2`/`L_PLATE`/`R_PLATE` columns are dropped
entirely rather than carried forward: they were assigned by partitioning against some
unrecorded static polygon set, and shipping them would silently mislead anyone who
reconstructs against a different plate model. `PaleoLithology()` therefore returns points
with no plate ID; callers must assign one themselves (e.g. via
`ReconstructionModel`/`pygplates.partition_into_plates` against their chosen static
polygons) before calling `reconstruct()`.
