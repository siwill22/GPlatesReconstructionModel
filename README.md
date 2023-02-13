# GPlatesReconstructionModel

*Prototyping for GPlates plate reconstructions as a python class*

- To use, clone the repo and add the main directory to the pythonpath
- Or, try typing the following in the command line:
> pip install --no-cache-dir --upgrade git+https://github.com/siwill22/GPlatesReconstructionModel
- pygplates will need to be installed separately
- Usage examples are available in the test_notebooks directory

## Python Dependencies
### Required:
- numpy
- scipy
- pandas
- geopandas
- matplotlib
- xarray
- pygplates
- PlateTectonicTools (https://github.com/EarthByte/PlateTectonicTools)
- pooch
- pygmt

### Optional:
- astropy_healpix (if astropy_healpix is not available, the some functions will fall back on precomputed point distributions in the 'data' folder)
- rasterio
- skimage
- stripy
- moviepy
- cartopy


