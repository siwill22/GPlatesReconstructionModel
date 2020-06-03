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
- healpy (if healpy is not available, the some functions will fall back on precomputed point distributions in the 'data' folder)
- stripy
- pygmt ~~gmt (accessed via calls to external command-line processes)~~
- moviepy
- cartopy and python basemap (but planned for removal with all plotting handled by matplotlib and pygmt)


## Warning
- Current version is still preliminary with minimal error handling and likely numerous bugs
- The current version in the process of being transitoined to python3 compatibility, which will temporarily break some functions 
