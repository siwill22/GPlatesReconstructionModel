# GPlatesReconstructionModel

*Prototyping for GPlates plate reconstructions as a python class*

- Usage examples are available in the test_notebooks directory

## Installation

Install directly from GitHub with pip:
```
pip install git+https://github.com/siwill22/GPlatesReconstructionModel
```

This pulls in every core dependency automatically, including `pygplates` and
`PlateTectonicTools`, both of which are ordinary PyPI packages.

To also install the optional dependencies for plotting, geophysics helpers, and additional
spatial algorithms, install the `all` extra:
```
pip install "gprm[all] @ git+https://github.com/siwill22/GPlatesReconstructionModel"
```

A plain install needs **no system packages** — every core dependency has wheels on PyPI for
Linux, macOS and Windows. That means the data-preparation and analysis path (fetching
reconstruction models, assigning plate IDs, reconstructing data, plate snapshots, and all
proximity calculations) runs on a headless machine or in CI with nothing but `pip install`.

`pygmt` is deliberately **not** a core dependency, because it is a wrapper around the GMT
command-line library and does not bundle GMT itself — that is the one thing pip cannot
install for you. It is only needed for plotting and for a few grid-sampling helpers, all of
which import it at the point of use and tell you what to install if it is missing.

## Python Dependencies
### Required (installed automatically, all pip-installable):
- numpy, scipy, pandas
- geopandas, shapely, rasterio
- matplotlib, xarray, cartopy
- pygplates
- PlateTectonicTools (https://github.com/EarthByte/PlateTectonicTools) — this itself
  imports `cartopy` unconditionally, which is why `cartopy` is required rather than optional
- pooch, requests, tqdm, pyyaml

### Optional extras:
- `pip install "gprm[viz]"` — pygmt and basemap (plotting, grid sampling). **pygmt also needs
  GMT >=6 installed separately**: `conda install -c conda-forge gmt`, Homebrew's `gmt`, or
  your Linux package manager.
- `pip install "gprm[geophysics]"` — pyshtools, litho1pt0, pmagpy
- `pip install "gprm[spatial]"` — stripy, astropy-healpix, scikit-image, scikit-learn
  (if astropy-healpix is not installed, some functions fall back on precomputed point
  distributions in the `Data` folder)
- `pip install "gprm[all]"` — everything above

### Where downloaded data goes
Datasets are cached with `pooch`, in a platform-dependent location. Ask for it rather than
writing it out by hand:
```python
from gprm.datasets import cache_path
cache_path()                                    # the cache root
cache_path('TorsvikCocks2017', 'CEED6_LAND.gpml')  # a file inside a fetched bundle
```

**Note on `stripy`**: `stripy` (pulled in directly by `spatial`, and transitively by
`geophysics` via `litho1pt0`) has no prebuilt wheel for Apple Silicon macOS on any Python
version, and none at all for Python 3.13+. On those platforms pip will try to build it from
source, which needs a Fortran compiler (e.g. `brew install gcc` on macOS, or
`conda install -c conda-forge stripy` as an alternative to pip for just that package).


