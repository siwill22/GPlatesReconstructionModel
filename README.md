# GPlatesReconstructionModel

*Prototyping for GPlates plate reconstructions as a python class*

- To use, clone the repo and add the main directory to the pythonpath

- Usage examples are available in the test_notebooks directory

## Python Dependencies
### Required:
- numpy
- scipy
- pandas
- matplotlib
- xarray
- pygplates
- PlateTectonicTools (https://github.com/EarthByte/PlateTectonicTools)

### Optional:
- healpy
- stripy
- gmt (accessed via calls to external command-line processes)
- moviepy
- cartopy and python basemap (for legacxy reasons, map plotting options are currently a mixture of both, in addition to gmt)


Current version is still preliminary with minimal error handling and likely numerous bugs
