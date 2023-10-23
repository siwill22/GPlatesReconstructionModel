'''
MIT License

Copyright (c) 2017-2023 Simon Williams

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
'''

from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
import geopandas as _gpd
import os as _os
#from gprm import ReconstructionModel as _ReconstructionModel


# TODO
# Add: Domeier and Torsvik
#      Shephard?? (or a version of Seton with Static Polygons)



def fetch_CaoToyRodinia(load=True, model_case='NNR'):
    '''
    Load Toy Billion-year reconstructions from Cao et al (2020), Tectonics
    doi: 10.1029/2020GC009244

    model_case options: 'NNR' [default], 'OV', 'SSL'
    '''
    
    fnames = _retrieve(
        url="https://zenodo.org/record/3854549/files/1000Myr_synthetic_tectonic_reconstructions.zip?download=1",
        known_hash="md5:b7ea40c77826ef5d5e3b99affa3e9d66",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='CaoToyRodinia'),
    )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('CaoToyRodinia_'.format(model_case))
    reconstruction_model.add_rotation_model('{:s}/Global_EB_250-0Ma_GK07_2017_ASM.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/triple_junction_superoceanic_plates.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/410-250_toy_introversion_simplified.rot'.format(dirname))
    reconstruction_model.add_rotation_model('{:s}/1000-410_toy_introversion_simplified.rot'.format(dirname))

    reconstruction_model.add_continent_polygons('{:s}/COBfile_1000_0_Toy_introversion.gpml'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/coastline_file_1000_250_new_valid_time.gpml'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/coastline_file_250_0_new_valid_time.gpml'.format(dirname))

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
    Load Rodinia reconstruction from Li et al (2008), 
    doi: 10.1016/j.precamres.2007.04.021
    and updated in Li et al (2013) Sedimentary Geology
    doi: 10.1016/j.sedgeo.2013.05.016
    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Li_etal_2008_RodiniaModel.zip",
        known_hash="sha256:e659371df79acfd7e599d0a358be0b154705b84d92388c042e6382ef78a3f4f6",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Li2008'),
    )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Li++2008')
    reconstruction_model.add_rotation_model('{:s}/Li_etal_2008_RodiniaModel/RodiniaModel_CompleteRotationFile.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Li_etal_2008_RodiniaModel/RodiniaBlocks_WithPlateIDColumnAndIDs.shp'.format(dirname))
    
    return reconstruction_model


def fetch_Li2023(load=True, model_case='East'):
    '''
    Load the 2000-540 Ma reconstruction from Li et al (2023)
    doi:10.1016/j.earscirev.2023.104336

    NOTE 
    This rotation model does not contain zero-time rotations, and as such
    can give strange results when extracting rotations that specify the 'from_time' as 0
    '''    

    fnames = _retrieve(
            url="https://ars.els-cdn.com/content/image/1-s2.0-S0012825223000259-mmc10.zip",
            #known_hash="sha256:1e205c867feaed1b796b098bff4d65ecc18d0d6e1f37291ca3161a32957373a9", 
            known_hash="sha256:f49bcc6feb85684900035205e3d99be424eff80ef73caedec8f5d72c5dc620d8",
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='Li2023'),
        )

    dirname = _os.path.split(fnames[0])[0]

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Li++2023_{:s}'.format(model_case))
    reconstruction_model.add_continent_polygons('{:s}/Continental_outlines.shp'.format(dirname))
    #reconstruction_model.add_static_polygons('{:s}/Continental_outlines.shp'.format(dirname))
    if model_case=='East':
        reconstruction_model.add_rotation_model('{:s}/EDRG_90E_2000-540Ma.rot'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/EDRG_boundary_90E_2000-540Ma.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/EDRG_topology_90E_2000-540Ma.gpml'.format(dirname))
    elif model_case=='West':
        reconstruction_model.add_rotation_model('{:s}/EDRG_90W_2000-540Ma.rot'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/EDRG_boundary_90W_2000-540Ma.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/EDRG_topology_90W_2000-540Ma.gpml'.format(dirname))
    
    # The continent polygons contain many polygons that could be problematic for plate id assignment:
    # - the valid time does not span any of the time range covered by the reconstruction, 
    # - the valid time does not extend to 0 Ma
    # we remove these and create a new file for plate id assignment
    sp = _gpd.read_file('{:s}/Continental_outlines.shp'.format(dirname))
    sp = sp.query('FROMAGE>601')
    banned_list = 'Great India|Australia before Paterson-Petermann orogeny|Australia_outline|Pure Australia|South China'
    sp = sp[~sp.NAME.str.contains(banned_list, na=False)]
    #sp['TOAGE']=-999
    sp.to_file('{:s}/Static_polygon_outlines.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Static_polygon_outlines.shp'.format(dirname))

    return reconstruction_model


def fetch_DomeierTorsvik2014(load=True):
    '''
    Load 250-410 Ma reconstruction from Domeier and Torsvik (2014) 
    doi:

    '''
    fnames = _retrieve(
        url="http://www.earthdynamics.org/data/Domeier2014_data.zip",
        known_hash="sha256:e2bd29afc9bf9cfdbcaba286787526dc71f4e0b5e05415b01e11a953d8abfe61",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='DomeierTorsvik2014'),
    )

    from gprm import ReconstructionModel as _ReconstructionModel

    dirname = '{:s}/DomeierTorsvik2014/'.format(str(_os_cache('gprm')))

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('DomeierTorsvik2014')
    reconstruction_model.add_rotation_model('{:s}/Domeier2014_data/LP_TPW.rot'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Domeier2014_data/LP_Land.shp'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/Domeier2014_data/LP_Land.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/Domeier2014_data/LP_Land.shp'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Domeier2014_data/LP_topos.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Domeier2014_data/LP_transform.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Domeier2014_data/LP_subduction.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Domeier2014_data/LP_ridge.gpml'.format(dirname))
    
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
        processor=_Unzip(extract_dir='Matthews2016'),
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
    Load Billion-year reconstruction from Merdith et al (2021) Earth Science Reviews
    doi: https://doi.org/10.1016/j.earscirev.2020.103477
    '''
    fnames = _retrieve(
        url="https://zenodo.org/record/4320873/files/SM2-Merdith_et_al_1_Ga_reconstruction.zip?download=1",
        known_hash="md5:1786d68e949c4242de1801388c68cb8c",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Merdith2021'),
    )

    dirname = '{:s}/Merdith2021/'.format(fnames[0].split('Merdith2021')[0])
    #for fname in fnames:
    #    if '__MACOSX' not in _os.path.split(fname)[0]:
    #        dirname = _os.path.split(fname)[0]
    #        break

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


def fetch_Muller2022(NNR=False, load=True):
    '''
    Load Billion-year reconstruction from Muller et al (2022) Solid Earth
    with optimised reference frame generated for Merdith et al (2021) model
    doi:
    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2022_SE/Muller_etal_2022_SE_1Ga_Opt_PlateMotionModel.zip",
        known_hash="sha256:a1e37f0a201a827ffe1150a808984786c6a2982923362580e28718fe7bc716e7",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Muller2022'),
    )

    dirname = '{:s}/Muller2022/Muller_etal_2022_SE_1Ga_Opt_PlateMotionModel/'.format(str(_os_cache('gprm')))

    from gprm import ReconstructionModel as _ReconstructionModel
    if NNR:
        reconstruction_model = _ReconstructionModel('Muller++2022_NNR')
        reconstruction_model.add_rotation_model('{:s}/optimisation/no_net_rotation_model.rot'.format(dirname))
    else:
        reconstruction_model = _ReconstructionModel('Muller++2022_Opt')
        reconstruction_model.add_rotation_model('{:s}/optimisation/1000_0_rotfile_Merdith_et_al_optimised.rot'.format(dirname))
    #reconstruction_model.add_static_polygons('{:s}/shapes_static_polygons_Merdith_et_al.gpml'.format(dirname))
    #reconstruction_model.add_coastlines('{:s}/'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/shapes_continents_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/410-250_plate_boundaries_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/250-0_plate_boundaries_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/TopologyBuildingBlocks_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/1000-410-Transforms_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/1000-410-Convergence_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/1000-410-Divergence_Merdith_et_al.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/1000-410-Topologies_Merdith_et_al.gpml'.format(dirname))

    return reconstruction_model


def fetch_Muller2016(load=True):
    '''
    Load Pangea breakup reconstruction from Muller et al (2016) Ann Rev Earth & Plan Sci
    doi: 10.1146/annurev-earth-060115-012211
    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2016_AREPS/Muller_etal_2016_AREPS_Supplement/Muller_etal_2016_AREPS_Supplement_v1.17.zip",
        known_hash="sha256:a671d6f2318b329e6f633065771fe37d29b6932e805e619039c4405dcb0fb91a",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Muller2016'),
    )

    #dirname = _os.path.split(fnames[0])[0]
    dirname = '{:s}/Muller2016/'.format(fnames[0].split('Muller2016')[0])


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
    Load Pangea breakup reconstruction from Muller et al, (2019) Tectonics
    doi: https://doi.org/10.1029/2018TC005462
    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Muller_etal_2019_Tectonics/Muller_etal_2019_PlateMotionModel/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics.zip",
        known_hash="sha256:32c30c80cd165fe0d28b3fda44a8b7d42e660a2a95baf508bdf7d1666977be9d",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Muller2019'),
    )

    dirname = '{:s}/Muller2019/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics/'.format(fnames[0].split('Muller2019')[0])
    #dirname = '{:s}/Muller_etal_2019_PlateMotionModel_v2.0_Tectonics/'.format(_os.path.split(fnames[0])[0])

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
    Load Nuna reconstruction from Pehrsson et al, (2015) Geol Soc London Spec Pub 
    doi: http://doi.org/10.1144/SP424.5
    '''
    fnames = _retrieve(
        url="https://www.geolsoc.org.uk/~/media/Files/GSL/shared/Sup_pubs/2015/18822_7.zip",
        known_hash="sha256:12e7ed7f1f736b0421a60c60151fed7b46ce028b3348f8bf39ba6d7916651b6f",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Pehrsson2015'),
    )

    #dirname = _os.path.split(fnames[0])[0]
    dirname = '{:s}/Pehrsson2015/'.format(fnames[0].split('Pehrsson2015')[0])


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
    Load Pangea breakup reconstruction from Seton et al (2012)
    doi:10.1016/j.earscirev.2012.03.002
    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Seton_etal_2012_ESR.zip",
        known_hash="sha256:b117354f93296dc1035d6709c7d475bf9ad517dc3f882b1621ef68db712c603e",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Seton2012'),
    )

    #dirname = _os.path.split(fnames[0])[0]
    dirname = '{:s}/Seton2012/'.format(fnames[0].split('Seton2012')[0])

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Seton++2012')
    reconstruction_model.add_rotation_model('{:s}/Seton_etal_2012_ESR/Rotations/Seton_etal_ESR2012_2012.1.rot'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/Seton_etal_2012_ESR/Coastlines/Seton_etal_ESR2012_Coastline_2012.1.gpml'.format(dirname))
    reconstruction_model.add_dynamic_polygons('{:s}/Seton_etal_2012_ESR/Plate_polygons/Seton_etal_ESR2012_PP_2012.1.gpml'.format(dirname))

    return reconstruction_model


def fetch_TorsvikCocks2017(load=True, wrap_static_polygons=True):
    '''
    Load Phanerozoic reconstruction from the book 'Earth History and Paleogeography'
    by Torvsik and Cocks (2017)
    doi: https://doi.org/10.1017/9781316225523

    NOTE: Terranes which only exist in the Paleozoic are included in the continents layer, 
    which could cause a problem if this layer is used for plate_id assignment of present-day
    data. Use the static polygon layer instead.

    '''
    fnames = _retrieve(
        url="http://www.earthdynamics.org/earthhistory/bookdata/CEED6.zip",
        known_hash="sha256:c2964d3ed791ae4b05690b16ff917309ae9dde60256be7e892df2463c2587eb1",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='TorsvikCocks2017'),
    )

    #dirname = _os.path.split(fnames[0])[0]
    dirname = '{:s}/TorsvikCocks2017/'.format(fnames[0].split('TorsvikCocks2017')[0])


    

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Torsvik+Cocks2017')
    reconstruction_model.add_rotation_model('{:s}/Torsvik_Cocks_HybridRotationFile.rot'.format(dirname))
    
    reconstruction_model.add_continent_polygons(('{:s}/CEED6_LAND.gpml'.format(dirname)))
    reconstruction_model.add_continent_polygons('{:s}/CEED6_TERRANES.shp'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/CEED6_MICROCONTINENTS.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/CEED6_LAND.gpml'.format(dirname))
    if wrap_static_polygons:
        import pygplates as _pygplates
        _pygplates.reconstruct(['{:s}/CEED6_MICROCONTINENTS.shp'.format(dirname),
                                '{:s}/CEED6_LAND.gpml'.format(dirname)],
                                [],
                                '{:s}/CEED6_Static_polygons.shp'.format(dirname),
                                0, export_wrap_to_dateline=True)
        reconstruction_model.add_static_polygons('{:s}/CEED6_Static_polygons.shp'.format(dirname))                        
    else:
        reconstruction_model.add_static_polygons('{:s}/CEED6_MICROCONTINENTS.shp'.format(dirname))
        reconstruction_model.add_static_polygons('{:s}/CEED6_LAND.gpml'.format(dirname))
    
    return reconstruction_model


def fetch_vanHinsbergen(load=True):
    '''
    Load global reconstructions compiled from Douwe van Hinsbergen's work

    NB CURRENTLY INCOMPLETE - NEED TO ESTABLISH A SET OF POLYGONS TO USE

    '''
    fnames = _retrieve(
        url="http://www.geologist.nl/wp-content/uploads/2019/09/vanHinsbergen_GPlates_reconstructions.zip",
        known_hash="sha256:7ed6319f11b4f4626c8211359cfeb8b454cb4381a81ee368fa11effbf06c1eeb",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='vanHinsbergen'),
    )

    #dirname = _os.path.split(fnames[0])[0]
    dirname = '{:s}/vanHinsbergen/'.format(fnames[0].split('vanHinsbergen')[0])

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('vanHinsbergen++2017')
    reconstruction_model.add_rotation_model('{:s}/'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/.shp'.format(dirname))
    reconstruction_model.add_coastlines('{:s}/.gpml'.format(dirname))
    
    return reconstruction_model


def fetch_Young2019(load=True):
    '''
    Load 0-410 Ma reconstruction from Young et al (2019) Geoscience Frontiers 
    doi: https://doi.org/10.1016/j.gsf.2018.05.011
    '''
    fnames = _retrieve(
        url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Young_etal_2018_GeoscienceFrontiers/Young_etal_2018_GeoscienceFrontiers_GPlatesPlateMotionModel.zip",
        known_hash="sha256:3cffdd988b802ad8961aad65901a95890a7b0058a3de3c353cf46986cca9f1f1",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Young2019'),
    )

    #for fname in fnames:
    #    if _os.path.split(fname)[1] == 'License.txt':
    #        dirname = _os.path.split(fname)[0]
    dirname = '{:s}/Young2019/'.format(fnames[0].split('Young2019')[0])

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

def fetch_Scotese(load=True):
    '''
    Load 0-, doi:

    '''
    fnames = _retrieve(
        url="https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip",
        known_hash="sha256:e01b19cee7c65a011ca4c42f187aba0ec24c1a87b842e2061eab9d22dc52ca80",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Cao2018_SM'),
    )

    #print(fnames[0].split('.unzip'))
    dirname = '{:s}/'.format(_os.path.split(fnames[0])[0])
    #dirname = '{:s}.unzip/SupplementaryMaterial_Cao_etal/'.format(fnames[0].split('.unzip')[0])

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Scotese2008')
    reconstruction_model.add_rotation_model('{:s}/SupplementaryMaterial_Cao_etal/Rotation_models/Scotese_2008_Rotation.rot'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/SupplementaryMaterial_Cao_etal/Rotation_models/Scotese_2008_PresentDay_ContinentalPolygons.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/SupplementaryMaterial_Cao_etal/Rotation_models/Scotese_2008_PresentDay_ContinentalPolygons.shp'.format(dirname))
    
    return reconstruction_model


def fetch_Golonka(load=True):
    '''
    Load reconstruction of Golonka, 2007, spanning the time range 0-5XX Ma, 
    doi:

    '''
    fnames = _retrieve(
        url="https://static.cambridge.org/content/id/urn:cambridge.org:id:article:S0016756818000110/resource/name/S0016756818000110sup001.zip",
        known_hash="sha256:e01b19cee7c65a011ca4c42f187aba0ec24c1a87b842e2061eab9d22dc52ca80",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
        processor=_Unzip(extract_dir='Cao2018_SM'),
    )
    
    dirname = '{:s}/'.format(_os.path.split(fnames[0])[0])

    from gprm import ReconstructionModel as _ReconstructionModel
    reconstruction_model = _ReconstructionModel('Golonka2007')
    reconstruction_model.add_rotation_model('{:s}/Rotation_models/Golonka_2007_Rotation.rot'.format(dirname))
    reconstruction_model.add_continent_polygons('{:s}/Rotation_models/Golonka_2007_PresentDay_ContinentalPolygons.shp'.format(dirname))
    reconstruction_model.add_static_polygons('{:s}/Rotation_models/Golonka_2007_PresentDay_ContinentalPolygons.shp'.format(dirname))
    
    return reconstruction_model


def fetch_Clennett(load=True, model_case='M2019'):
    """
    Load reconstruction files associated with the study of Clennett et al (2020),
    Geochemistry, Geophysics, Geosystems
    doi: https://doi.org/10.1029/2020GC009117

    model case must be either:
     - 'M2019' (default, version based on deforming model of Muller et al 2019), or 
     - 'S2013' (rigid topological model based on Shephard et al, 2013)
    """
    if model_case=='M2019':
        fnames = _retrieve(
            url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_M2019.zip",
            known_hash="sha256:1ad6a29ceb396b581930734b1f6e8409e52dc4e2ae9658156ac2dd732cb82ab8",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='Clennett2020_M2019'),
        )

        #dirname = _os.path.split(fnames[0])[0]
        #dirname = '{:s}/Clennett_etal_2020_M2019/'.format(_os.path.split(fnames[0])[0])
        dirname = '{:s}/Clennett2020_M2019/Clennett_etal_2020_M2019/'.format(fnames[0].split('Clennett2020_M2019')[0])

        # if downloading for first time, remove the unwanted MeshPoint files
        if _os.path.isdir('{:s}/DeformingMeshPoints'.format(dirname)):
            import shutil
            shutil.rmtree('{:s}/DeformingMeshPoints'.format(dirname))

        from gprm import ReconstructionModel as _ReconstructionModel
        reconstruction_model = _ReconstructionModel('Clennett++M2019')
        
        reconstruction_model.add_rotation_model('{:s}/Global_250-0Ma_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/Clennett_etal_2020_Rotations.rot'.format(dirname))

        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/Alps_Mesh_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/Andes_Flat_Slabs_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/Andes_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/Australia_Antarctica_Mesh_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/Australia_North_Zealandia_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/Eurasia_Arabia_Mesh_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/North_America_Flat_Slabs_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/North_America_Mesh_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/North_China_Mesh_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/South_Atlantic_Rotations_2019_v2.rot'.format(dirname))
        reconstruction_model.add_rotation_model('{:s}/DeformingMeshes/Southeast_Asia_Rotations_2019_v2.rot'.format(dirname))

        reconstruction_model.add_continent_polygons('{:s}/Clennett_etal_2020_Terranes.gpml'.format(dirname))
        #reconstruction_model.add_static_polygons('{:s}/StaticGeometries/StaticPolygons/Global_EarthByte_GPlates_PresentDay_StaticPlatePolygons_2019_v1.shp'.format(dirname))
        reconstruction_model.add_coastlines('{:s}/Clennett_etal_2020_Coastlines.gpml'.format(dirname))
        
        reconstruction_model.add_dynamic_polygons('{:s}/Clennett__etal_2020_NAm_boundaries.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/Clennett_etal_2020_Plates.gpml'.format(dirname))

        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Alps_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Alps_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/America_Anyui_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/America_Anyui_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Andes_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Andes_Flat_Slabs_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Andes_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Arctic_Eurasia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Australia_Antarctica_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Australia_Antarctica_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Australia_North_Zealandia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Australia_North_Zealandia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Baja_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Coral_Sea_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Coral_Sea_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/East_African_Rift_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/East-West_Gondwana_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Ellesmere__Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Eurasia_Arabia_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Greater_India_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Greater_India_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        #reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Inactive_Meshes_and_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/North_America_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/North_Atlantic_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/North_Atlantic_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/North_China_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/North_China_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Northern_Andes_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Northern_Andes_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Papua_New_Guinea_Deforming_Meshes_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Papua_New_Guinea_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Scotia_Deforming_Mesh_and_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Siberia_Eurasia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Siberia_Eurasia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/South_Atlantic_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/South_Atlantic_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/South_China_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/South_China_Sea_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/South_Zealandia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/South_Zealandia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Southeast_Asia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Southeast_Asia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/West_Antarctic_Zealandia_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/West_Antarctica_Zealandia_Mesh_Topologies_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Western_North_America_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Western_Tethys_Deforming_Mesh_2019_v2.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/DeformingMeshes/Western_Tethys_Tectonic_Boundary_Topologies_2019_v2.gpml'.format(dirname))

        return reconstruction_model

    elif model_case=='S2013':
        fnames = _retrieve(
            url="https://www.earthbyte.org/webdav/ftp/Data_Collections/Clennett_etal_2020_G3/Clennett_etal_2020_S2013.zip",
            known_hash="sha256:7749aac19c2d07c80a2cd77ba6a9038c8911fa8804c1d4adb3c9da7cb635b691",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
            processor=_Unzip(extract_dir='Clennett2020_S2013'),
        )

        #dirname = '{:s}/Clennett_etal_2020_S2013/'.format(_os.path.split(fnames[0])[0])
        dirname = '{:s}/Clennett2020_S2013/Clennett_etal_2020_S2013/'.format(fnames[0].split('Clennett2020_S2013')[0])


        from gprm import ReconstructionModel as _ReconstructionModel
        reconstruction_model = _ReconstructionModel('Clennett++S2013')
        
        reconstruction_model.add_rotation_model('{:s}/Clennett_etal_2020_Rotations.rot'.format(dirname))

        reconstruction_model.add_continent_polygons('{:s}/Clennett_etal_2020_Terranes.gpml'.format(dirname))
        reconstruction_model.add_coastlines('{:s}/Clennett_etal_2020_Coastlines.gpml'.format(dirname))
        
        reconstruction_model.add_dynamic_polygons('{:s}/Clennett_etal_2020_NAm_boundaries.gpml'.format(dirname))
        reconstruction_model.add_dynamic_polygons('{:s}/Clennett_etal_2020_Plates.gpml'.format(dirname))

        return reconstruction_model 

    else:
        ValueError('Unrecognised model name {}'.format(model_case))

    
    