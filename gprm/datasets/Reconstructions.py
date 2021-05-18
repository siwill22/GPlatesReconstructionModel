from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
#import geopandas as _gpd
import os as _os
from gprm import ReconstructionModel as _ReconstructionModel



def fetch_Pehrsson2015(load=True):
    '''
    Load Nuna reconstruction from Pehrsson et al, doi:

    '''
    fnames = _retrieve(
        url="https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip",
        known_hash="sha256:12e7ed7f1f736b0421a60c60151fed7b46ce028b3348f8bf39ba6d7916651b6f",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = _os.path.split(fnames[0])[0]

    # process the rotation file to enable continuous reconstruction
    df = _pd.read_fwf('{:s}/T_Rot_Model_Abs_25Ma_20131004.rot'.format(dirname),
                     colspecs=[(0,7),(7,17),(17,28),(28,40),(40,51),(51,56),(56,81)],
                     names=['MovingPlate','Time','Lat','Long','Angle','FixedPlate','Comment'])

    dfs = df.sort_values(['MovingPlate','Time'])

    dfs.to_csv('{:s}/T_Rot_Model_Abs_25Ma_20131004_sort.rot'.format(dirname),
            sep=' ',
            header=False,
            index=False)
    
    reconstruction_model = _ReconstructionModel('Pehrsson++2015')
    reconstruction_model.add_rotation_model('{:s}/T_Rot_Model_Abs_25Ma_20131004_sort.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/PlatePolygons.shp'.format(dirname))

    return reconstruction_model


def fetch_Muller2016(load=True):
    '''
    Load Pangea breakup reconstruction from Muller et al, doi:

    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp_data/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip",
        known_hash="sha256:a671d6f2318b329e6f633065771fe37d29b6932e805e619039c4405dcb0fb91a",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = _os.path.split(fnames[0])[0]

    reconstruction_model = _ReconstructionModel('Muller++2016')
    reconstruction_model.add_rotation_model('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Shapefiles/StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_2015_v1.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_Coastlines.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_Topology_BuildingBlocks.gpml'.format(dirname))

    return reconstruction_model

