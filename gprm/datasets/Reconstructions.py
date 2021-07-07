from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
#import geopandas as _gpd
import os as _os
#from gprm import ReconstructionModel as _ReconstructionModel


# TODO
# Add: Domeier and Torsvik
#      Shephard?? (or a version of Seton with Static Polygons)
#      Scotese
#      Golonka??


def fetch_CaoToyRodinia(load=True, model_case='NNR'):
    '''
    Load Toy Billion-year reconstructions from Cao et al (2020), doi:

    '''
    fnames = _retrieve(
        url="https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip?download=1",
        known_hash="md5:b7ea40c77826ef5d5e3b99affa3e9d66",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel()
    reconstruction_model.add_rotation_model('{:s}/Global_EB_250-0Ma_GK07_2017_ASM.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/triple_junction_superoceanic_plates.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/410-250_toy_introversion_simplified.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/1000-410_toy_introversion_simplified.rot'.format(dirname))

    reconstruction_model.add_continent_polygons('{:s}/COBfile_1000_0_Toy_introversion.gpml'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/coastline_file_1000_250_new_valid_time.gpml'.format(dirname))

    reconstruction_model.add_dynamic_polygons('{:s}/Global_EarthByte_Mesozoic-Cenozoic_plate_boundaries_2016_v5.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/TopologyBuildingBlocks_AREPS.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Toy_introversion_plate_boundaries_410_250_new_valid_time.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Toy_introversion_plate_boundaries_1000_410_new_valid_time.gpml'.format(dirname))

    if model_case == 'NNR':
        reconstruction_model.add_rotation_model('{:s}/NLR_SLOW_CONTINENT_0Ma_1000Ma_NNR.rot'.format(dirname))
    elif model_case == 'OV':
        reconstruction_model.add_rotation_model('{:s}/NLR_SLOW_CONTINENT_0Ma_1000Ma_OV.rot'.format(dirname))
    elif model_case == 'SSL':
        reconstruction_model.add_rotation_model('{:s}/NLR_SLOW_CONTINENT_0Ma_1000Ma_SSL.rot'.format(dirname))
    else:
        ValueError('Unrecognised model name {}'.format(model_case))

    return reconstruction_model


def fetch_Li2008(load=True):
    '''
    Load Rodinia reconstruction from Li et al (2008), doi:

    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp_data/Data_Collections/Li_etal_2008_RodiniaModel.zip",
        known_hash="sha256:e659371df79acfd7e599d0a358be0b154705b84d92388c042e6382ef78a3f4f6",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Li++2008')
    reconstruction_model.add_rotation_model('{:s}/Li_etal_2008_RodiniaModel/RodiniaModel_CompleteRotationFile.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Li_etal_2008_RodiniaModel/RodiniaBlocks_WithPlateIDColumnAndIDs.shp'.format(dirname))
    
    return reconstruction_model


def fetch_Matthews2016(load=True):
    '''
    Load 0-410 Ma reconstruction from Matthews et al, doi:

    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Matthews_etal_2016_Global_Plate_Model_GPC.zip",
        known_hash="sha256:c88acba32f7e5a00734d14d8c512a20392dc8e62d75fd1777d351eb7e6ada28f",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    for fname in fnames:
        if _os.path.split(fname)[1] == 'License.txt':
            dirname = _os.path.split(fname)[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Matthews++2016')
    reconstruction_model.add_rotation_model('{:s}/Matthews_etal_2016_Global_Plate_Model_GPC/Global_EB_250-0Ma_GK07_Matthews_etal.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Matthews_etal_2016_Global_Plate_Model_GPC/Global_EB_410-250Ma_GK07_Matthews_etal.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Matthews_etal_2016_Global_Plate_Model_GPC/StaticGeometries/StaticPolygons/PresentDay_StaticPlatePolygons_Matthews++.shp'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/Matthews_etal_2016_Global_Plate_Model_GPC/StaticGeometries/ContinentalPolygons/PresentDay_ContinentalPolygons_Matthews++.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/Matthews_etal_2016_Global_Plate_Model_GPC/StaticGeometries/Coastlines/Global_coastlines_2015_v1_low_res.shp'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Matthews_etal_2016_Global_Plate_Model_GPC/Global_EarthByte_Paleozoic_plate_boundaries_Matthews_etal.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Matthews_etal_2016_Global_Plate_Model_GPC/Global_EarthByte_Mesozoic-Cenozoic_plate_boundaries_Matthews_etal.gpml'.format(dirname))

    return reconstruction_model


def fetch_Merdith2021(load=True):
    '''
    Load Billion-year reconstruction from Merdith et al (2021), doi:

    '''
    fnames = _retrieve(
        url="https://zenodo.org/record/4320873/files/SM2-Merdith_et_al_1_Ga_reconstruction.zip?download=1",
        known_hash="md5:1786d68e949c4242de1801388c68cb8c",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    for fname in fnames:
        if '__MACOSX' not in _os.path.split(fname)[0]:
            dirname = _os.path.split(fname)[0]
            break

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Merdith++2021')
    reconstruction_model.add_rotation_model('{:s}/SM2/1000_0_rotfile_Merdith_et_al.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/SM2/shapes_static_polygons_Merdith_et_al.gpml'.format(dirname))
    #reconstruction_model.add_coastlines('{:s}/'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/SM2/shapes_continents_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/SM2/410-250_plate_boundaries_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/SM2/250-0_plate_boundaries_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/SM2/TopologyBuildingBlocks_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/SM2/1000-410-Transforms_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/SM2/1000-410-Convergence_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/SM2/1000-410-Divergence_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/SM2/1000-410-Topologies_Merdith_et_al.gpml'.format(dirname))

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

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Muller++2016')
    reconstruction_model.add_rotation_model('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Shapefiles/StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_2015_v1.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_Coastlines.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_PlateBoundaries.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Muller_etal_2016_AREPS_Supplement_v1.17/Global_EarthByte_230-0Ma_GK07_AREPS_Topology_BuildingBlocks.gpml'.format(dirname))

    return reconstruction_model


def fetch_Muller2019(load=True):
    '''
    Load Pangea breakup reconstruction from Muller et al, doi:

    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp_data/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip",
        known_hash="sha256:6e8d193f61ebeaa2f68cc55afe399a654bf31a8c564d5305245c91c32162814c",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = '{:s}/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics/'.format(_os.path.split(fnames[0])[0])

    # if downloading for first time, remove the unwanted MeshPoint files
    if _os.path.isdir('{:s}/DeformingMeshPoints'.format(dirname)):
        import shutil
        shutil.rmtree('{:s}/DeformingMeshPoints'.format(dirname))

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Muller++2019')
    reconstruction_model.add_rotation_model('{:s}/Alps_Mesh_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Andes_Flat_Slabs_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Andes_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Australia_Antarctica_Mesh_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Australia_North_Zealandia_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Eurasia_Arabia_Mesh_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Global_250-0Ma_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/North_America_Flat_Slabs_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/North_America_Mesh_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/North_China_Mesh_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/South_Atlantic_Rotations_2019_v2.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Southeast_Asia_Rotations_2019_v2.rot'.format(dirname))

    reconstruction_model.add_continent_polygons('{:s}/StaticGeometries/ContinentalPolygons/Global_EarthByte_GPlates_PresentDay_ContinentalPolygons_2019_v1.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/StaticGeometries/StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_2019_v1.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/StaticGeometries/Coastlines/Global_coastlines_2019_v1_low_res.shp'.format(dirname))
    
    reconstruction_model.add_dynamic_polygons('{:s}/Alps_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Alps_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/America_Anyui_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/America_Anyui_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Andes_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Andes_Flat_Slabs_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Andes_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Arctic_Eurasia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Australia_Antarctica_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Australia_Antarctica_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Australia_North_Zealandia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Australia_North_Zealandia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Baja_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Coral_Sea_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Coral_Sea_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/East_African_Rift_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/East-West_Gondwana_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Ellesmere__Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Eurasia_Arabia_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Global_Mesozoic-Cenozoic_PlateBoundaries_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Greater_India_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Greater_India_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    #reconstruction_model.add_dynamic_polygons('{:s}/Inactive_Meshes_and_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/North_America_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/North_Atlantic_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/North_Atlantic_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/North_China_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/North_China_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Northern_Andes_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Northern_Andes_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Papua_New_Guinea_Deforming_Meshes_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Papua_New_Guinea_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Scotia_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Siberia_Eurasia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Siberia_Eurasia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/South_Atlantic_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/South_Atlantic_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/South_China_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/South_China_Sea_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/South_Zealandia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/South_Zealandia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Southeast_Asia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Southeast_Asia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/West_Antarctic_Zealandia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/West_Antarctica_Zealandia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Western_North_America_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Western_Tethys_Deforming_Mesh_2019_v2.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Western_Tethys_Tectonic_Boundary_Topologies_2019_v2.gpml'.format(dirname))

    return reconstruction_model


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
    
    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Pehrsson++2015')
    reconstruction_model.add_rotation_model('{:s}/T_Rot_Model_Abs_25Ma_20131004_sort.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/PlatePolygons.shp'.format(dirname))

    return reconstruction_model


def fetch_Seton2012(load=True):
    '''
    Load Pangea breakup reconstruction from Seton et al (2012), doi:

    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp_data/Data_Collections/Seton_etal_2012_ESR.zip",
        known_hash="sha256:b117354f93296dc1035d6709c7d475bf9ad517dc3f882b1621ef68db712c603e",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Seton++2012')
    reconstruction_model.add_rotation_model('{:s}/Seton_etal_2012_ESR/Rotations/Seton_etal_ESR2012_2012.1.rot'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1.gpml'.format(dirname))

    return reconstruction_model


def fetch_TorsvikCocks2017(load=True):
    '''
    Load Phanerozoic reconstruction from Torvsik and Cocks (2017), doi:

    '''
    fnames = _retrieve(
        url="http://www.earthdynamics.org/earthhistory/bookdata/CEED6.zip",
        known_hash="sha256:9b6d6f8a9a6299a269fd16f07aeb48dc0b4d591743d6691b86fde7b550d1ce7b",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Torsvik+Cocks2017')
    reconstruction_model.add_rotation_model('{:s}/Torsvik_Cocks_HybridRotationFile.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/CEED6_TERRANES.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/CEED6_MICROCONTINENTS.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/CEED6_LAND.gpml'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/CEED6_LAND.gpml'.format(dirname))
    
    return reconstruction_model


def fetch_vanHinsbergen(load=True):
    '''
    Load global reconstructions compiled from Douwe van Hinsbergen's work

    '''
    fnames = _retrieve(
        url="http://www.geologist.nl/wp-content/uploads/2019/09/vanHinsbergen_GPlates_reconstructions.zip",
        known_hash="sha256:7ed6319f11b4f4626c8211359cfeb8b454cb4381a81ee368fa11effbf06c1eeb",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('vanHinsbergen++2017')
    reconstruction_model.add_rotation_model('{:s}/'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/.gpml'.format(dirname))
    
    return reconstruction_model


def fetch_Young2019(load=True):
    '''
    Load 0-410 Ma reconstruction from Young et al, doi:

    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Young_etal_2018_GeoscienceFrontiers/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel.zip",
        known_hash="sha256:3cffdd988b802ad8961aad65901a95890a7b0058a3de3c353cf46986cca9f1f1",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(),
    )

    for fname in fnames:
        if _os.path.split(fname)[1] == 'License.txt':
            dirname = _os.path.split(fname)[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Young++2019')
    reconstruction_model.add_rotation_model('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/Global_410-250Ma_Young_et_al.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/Global_250-0Ma_Young_et_al.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/StaticPolygons/Global_GPlates_PresentDay_StaticPlatePolygons_Young_et_al.shp'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/ContinentalPolygons/PresentDay_ContinentalPolygons_Young_et_al.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/Coastlines/Global_coastlines_Young_et_al_low_res.shp'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/Global_Mesozoic-Cenozoic_plate_boundaries_Young_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/Global_Paleozoic_plate_boundaries_Young_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel/TopologyBuildingBlocks_Young_et_al.gpml'.format(dirname))

    return reconstruction_model

