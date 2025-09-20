from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
import geopandas as _gpd
import os as _os
import numpy as _np


def loadDB(version=2021):
    '''
    Load the Puetz (2018) database, returning two dataframes:
     - the Sample Details
     - the Data table
    '''

    if version==2018:
        fname = _retrieve(
            url="https://ars.els-cdn.com/content/image/1-s2.0-S1674987117302141-mmc1.xlsx",
            known_hash="sha256:8ea19b08d5c8d3c6e7f3239471d27b6da440fcfa40994e712ac9ae95642fe3d9",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )

        # load the whole xls file once to avoid duplication
        xls = _pd.ExcelFile(fname)
            
        # Get dataframes for the sheets with information per site and per sample
        df_Data = _pd.read_excel(xls, sheet_name='Data',
                                          dtype={'Sample Key': _np.int64,
                                                  'Sample_ID': 'string',
                                                  'Sample & Grain': 'string',
                                                  u'206Pb /\n238U\nAge\n(Ma)': _np.float64,
                                                  u'206Pb /\n238U\n2σ\nPrecis': _np.float64,
                                                  u'207Pb /\n235U\nAge\n(Ma)': _np.float64,
                                                  u'207Pb /\n235U\n2σ\nPrecis': 'string',
                                                  u'207Pb /\n206Pb\nAge\n(Ma)': 'string',
                                                  u'207Pb /\n206Pb\n2σ\nPrecis': _np.float64})
        df_Data[u'207Pb /\n235U\n2σ\nPrecis']=df_Data[u'207Pb /\n235U\n2σ\nPrecis'].replace(',','.',regex=True).astype(float)
        df_Data[u'207Pb /\n206Pb\nAge\n(Ma)']=df_Data[u'207Pb /\n206Pb\nAge\n(Ma)'].replace(',','.',regex=True).astype(float)

        df_SampleDetails = _pd.read_excel(xls, sheet_name='Sample_Details', 
                                          dtype={'Sample Key': _np.int64,
                                                'Reference Key': _np.int64,
                                                'Sample_ID': 'string',
                                                'Continent': 'string',
                                                'Country': 'string',
                                                'Region': 'string',
                                                'Major Geographic-Geologic Description': 'string',
                                                'Minor Geologic-Geographic Unit': 'string',
                                                'Locality': 'string',
                                                'Latitude': _np.float64,
                                                'Longitude': _np.float64, 
                                                'Est. Depos. Age (Ma)': _np.float64,
                                                'Max. Depos. Age (Ma)': _np.float64,
                                                'Table': 'string',
                                                'Mineral Tested': 'string',
                                                'Mass Spectrometer': 'string',
                                                'Class-1 Rock Type': 'string',
                                                'Class-2 Rock Type': 'string',
                                                'Class-3 Rock Type': 'string',
                                                'Sample Count': _np.int64})
        
        # rename some fields for neatness
        df_Data.rename(columns = {u'206Pb /\n238U\nAge\n(Ma)': '206Pb_238U_Age_Ma',
                                u'206Pb /\n238U\n2σ\nPrecis': '206Pb_238U_Precis',
                                u'207Pb /\n235U\nAge\n(Ma)': '207Pb_235U_Age_Ma',
                                u'207Pb /\n235U\n2σ\nPrecis': '207Pb_235U_Precis',
                                u'207Pb /\n206Pb\nAge\n(Ma)': '207Pb_206Pb_Age_Ma',
                                u'207Pb /\n206Pb\n2σ\nPrecis': '207Pb_206Pb_Precis'},
                    inplace = True)
        
        df_SampleDetails = df_SampleDetails.rename(columns = {'Est. Depos. Age (Ma)': 'Est_Depos_Age_Ma',
                                                              'Max. Depos. Age (Ma)': 'Max_Depos_Age_Ma'})
        df_SampleDetails = df_SampleDetails.dropna(subset=['Sample Key'])

        return df_SampleDetails, df_Data


    elif version==2019:
        '''
        Files from Puetz and Condie (2019) Geoscience Frontiers
        https://www.sciencedirect.com/science/article/pii/S1674987119300751#appsec1
        '''

        fname = _retrieve(
            url="https://ars.els-cdn.com/content/image/1-s2.0-S1674987119300751-mmc4.xlsx",
            known_hash="sha256:17128e11dfc340f932ae322e88ea09977dbc9abed103a6f6effb7cd5541376c9",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )

        xls = _pd.ExcelFile(fname)
        df = _pd.read_excel(xls, sheet_name='U_Pb_Detrital_Zircon')

        df = df.rename(columns = {'GPS Longitude': 'Longitude',
                                  'GPS Latitude': 'Latitude'})

        gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

        return gdf


    elif version==2021:
        # Supplementary table 4 from Puetz et al, 2021, ESR
        # https://www.sciencedirect.com/science/article/pii/S0012825221002464#s0050

        fname = _retrieve(
            url="https://ars.els-cdn.com/content/image/1-s2.0-S0012825221002464-mmc4.xlsx",
            known_hash="sha256:afdcaa26698ba06f113f29fab311af0ec37f290db117fc85847966c28b78ff09",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )

        xls = _pd.ExcelFile(fname)
        df = _pd.read_excel(xls, sheet_name='UPb_Data',
                            dtype={'Database': 'string',
                                   'Ref-Sample Key': 'string',
                                   'Sample&Grain': 'string',
                                   'Concord Class        (1 to 5)': _np.int64,
                                   'Base Concord delineator': _np.float64,
                                   'Minimum Segmented Disc. (myr)': _np.float64,
                                   'Probability 207Pb/206Pb age is correct': _np.float64,
                                   'Non-iter          age                     (Ma)': _np.float64,
                                   '6x6    Seq. Index': _np.int64,
                                   '6x6 Grid Wt.': _np.float64,
                                   '12x12    Seq. Index': _np.int64,
                                   '12x12 Grid Wt.': _np.float64,
                                   'Country Wts.': _np.float64,
                                   'Region Wts.': _np.float64,
                                   'Country/Area': 'string',
                                   'Major Region': 'string',
                                   'Continent': 'string',
                                   'Major Geographic-Geologic Description': 'string',
                                   'Minor Geologic-Geographic Unit': 'string',
                                   'Locality': 'string',
                                   'Latitude': _np.float64,
                                   'Longitude': _np.float64,
                                   'Max. Depos. Age': _np.int64,
                                   'Est. Depos. Age': _np.float64,
                                   'Min. Depos. Age': _np.float64,
                                   'Rock Type':'string'})

        df = df.rename(columns = {'Est. Depos. Age': 'Est_Depos_Age_Ma',
                                  'Max. Depos. Age': 'Max_Depos_Age_Ma',
                                  'Min. Depos. Age': 'Min_Depos_Age_Ma',
                                  'Non-iter          age                     (Ma)': 'Non_Iter_Age_Ma'})

        gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

        return gdf
    
    elif version==2024:
        # Supplementary tables 1-3 from Puetz et al, 2024, Scientific Data
        # https://www.nature.com/articles/s41597-023-02902-9

        fname1 = _retrieve(
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41597-023-02902-9/MediaObjects/41597_2023_2902_MOESM1_ESM.xlsx",
            known_hash="sha256:37a882b39e2e17c58951f8d25e508ce5d99e39dbcc4e5e11c877eb6cac7501b6",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )
        xls1 = _pd.ExcelFile(fname1)

        fname2 = _retrieve(
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41597-023-02902-9/MediaObjects/41597_2023_2902_MOESM2_ESM.xlsx",
            known_hash="sha256:a30bb122cd761fd5c28d9811fd54d633dc722bc61b351b9caafad7395a55edd4",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )
        xls2 = _pd.ExcelFile(fname2)

        fname3 = _retrieve(
            url="https://static-content.springer.com/esm/art%3A10.1038%2Fs41597-023-02902-9/MediaObjects/41597_2023_2902_MOESM3_ESM.xlsx",
            known_hash="sha256:52ada19c30f809852cbafc53b5e48e9f637cd0ebd66f37b9b5472e602cf4dac6",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )
        xls3 = _pd.ExcelFile(fname3)

        print('Concatenating databases, please be patient....')
        df1_merge = _pd.merge(xls1.parse('Samples'),xls1.parse('UPb_Data'),on='Ref-Sample Key')
        df2_merge = _pd.merge(xls2.parse('Samples'),xls2.parse('UPb_Data'),on='Ref-Sample Key')
        df3_merge = _pd.merge(xls3.parse('Samples'),xls3.parse('UPb_Data'),on='Ref-Sample Key')
        df_full = _pd.concat([df1_merge, df2_merge, df3_merge])

        # TODO clean up columns and column names
        gdf = _gpd.GeoDataFrame(df_full, geometry=_gpd.points_from_xy(df_full.Longitude, df_full.Latitude), crs=4326)

        return gdf


def load_Hf():
    # Supplementary table 5 from Puetz et al, 2021, ESR
    # https://www.sciencedirect.com/science/article/pii/S0012825221002464#s0050

    fname = _retrieve(
        url="https://ars.els-cdn.com/content/image/1-s2.0-S0012825221002464-mmc5.xlsx",
        known_hash="sha256:7402470ca7e7319899d949207da3eae8c0ba20bdf21050a85ca70ff6c9be9b8c",  
        downloader=_HTTPDownloader(progressbar=True),
        path=_os_cache('gprm'),
    )

    xls = _pd.ExcelFile(fname)
    df = _pd.read_excel(xls, sheet_name='Hf_Data')

    gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return gdf


# These functions are specifically for the 2018 version, to match sample coordinates against age distributions
def get_igneous_samples(df_SampleDetails=None,df_Data=None,version=2018):
    
    if version==2018:
        if df_SampleDetails is None:
            df_SampleDetails, df_Data = loadDB(version=2018)

        # igneous samples
        df_IgneousSamples = df_SampleDetails[df_SampleDetails['Class-1 Rock Type'].str.contains('igneous')]
        df = _pd.merge(df_IgneousSamples,df_Data,on='Sample Key')
        IgneousZircons = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

        IgneousZircons['Sample_ID_x'] = IgneousZircons['Sample_ID_x'].astype("string")
        IgneousZircons['Sample_ID_y'] = IgneousZircons['Sample_ID_y'].astype("string")

    elif version==2019:
        fname = _retrieve(
            url="https://ars.els-cdn.com/content/image/1-s2.0-S1674987119300751-mmc5.xlsx",
            known_hash="sha256:10329fcb555d6d4cc7adbe8bd55f5306646e7f4edad167cb45285aee2a630219",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )

        xls = _pd.ExcelFile(fname)
        df = _pd.read_excel(xls, sheet_name='U_Pb_Igneous_Zircon',
                            dtype={'Sample-Ref ID': 'string',
                                   'Country': 'string',
                                   'Region': 'string',
                                   'Continent': 'string',
                                   'GPS Latitude': _np.float64,
                                   'GPS Longitude': _np.float64,
                                   'Minor Geologic-Geographic Unit': 'string',
                                   'Rock-type 3': 'string',
                                   'Mafic or Felsic ?': 'string',
                                   'Spectrometer': 'string',
                                   'Best Age (Ma)': _np.float64,
                                   'Uncertainty (2σ precision)': _np.float64,
                                   '235U Abs. Disc. (myr)': _np.float64,
                                   'Sample & Grain': 'string',
                                   'Core/Rim': 'string',
                                   '206Pb/238U Age (Ma)': _np.float64,
                                   'Uncert. (2σ)': _np.float64,
                                   '207Pb/235U Age (Ma)': _np.float64,
                                   'Uncert. (2σ).1': _np.float64,
                                   '207Pb/206Pb Age (Ma)': _np.float64,
                                   'Uncert. (2σ).2': _np.float64})

        df = df.rename(columns = {'GPS Longitude': 'Longitude',
                                  'GPS Latitude': 'Latitude'})

        # rename some fields for neatness
        df = df.rename(columns = {u'206Pb/238U Age (Ma)': '206Pb_238U_Age_Ma',
                                u'Uncert. (2σ)': '206Pb_238U_Precis',
                                u'207Pb/235U Age (Ma)': '207Pb_235U_Age_Ma',
                                u'Uncert. (2σ).1': '207Pb_235U_Precis',
                                u'207Pb/206Pb Age (Ma)': '207Pb_206Pb_Age_Ma',
                                u'Uncert. (2σ).2': '207Pb_206Pb_Precis'})

        IgneousZircons = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return IgneousZircons


def get_sedimentary_samples(df_SampleDetails=None,df_Data=None,version=2018):
    '''
    For the 2018 database, the loading could be faster if you supply df_SampleDetails and df_Data
    This is only really an advantage if the aim is to load the 2018 Sedimentary and Igneous Zircons together 
    '''
    
    if version==2018:
        if df_SampleDetails is None:
            df_SampleDetails, df_Data = loadDB(version=2018)

        df_SamplesWithDepositionalAge = df_SampleDetails.dropna(subset=['Est_Depos_Age_Ma'])
        df = _pd.merge(df_SamplesWithDepositionalAge,df_Data,on='Sample Key')
        SedimentaryZircons = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    elif version==2019:

        SedimentaryZircons = loadDB(version=2019)

    elif version==2021:

        SedimentaryZircons = loadDB(version=2021)

    elif version==2024:

        SedimentaryZircons = loadDB(version=2024)


    return SedimentaryZircons


###########################
# Some functions for data manipulation and classification

def get_record_spectrum(df,Sample_ID):
    return df[df.Sample_ID_x==Sample_ID].dropna(subset=['206Pb_238U_Age_Ma'])


# Tectonic classifications
def tectonic_category(SedimentaryZircons, 
                      sample_key='Ref-Sample Key',
                      grain_age_key='Non_Iter_Age_Ma',
                      depositional_age_key='Est_Depos_Age_Ma'):
    ''' 
    Given a set of zircon spectra, return categories according to the 
    method of Cawood et al (2012). See also Jian et al (2022)
    '''
    cdf_markers = _np.arange(0,1.0001,0.01)

    sample_groups = SedimentaryZircons.groupby(by=sample_key)

    category_list = []

    for sample_group in sample_groups:
    
        s = sample_group[1][grain_age_key]
    
        age_at_deposition = s - sample_group[1][depositional_age_key]

        # define the cdf
        dst = _np.sort(age_at_deposition)
        xtmp = _np.linspace(0,1,len(dst))

        # sample cdf at regular increment
        cdf_vals = _np.interp(cdf_markers,xtmp,dst)
        #cdf_vals_ma = _np.interp(cdf_markers_ma,dst,xtmp)

        # classify according to Cawood et al (2012)
        print(cdf_vals[1], cdf_vals[6])
        print('This looks like it wrongly assumes a 0.05 spacing in the cdf markers')
        category = 'C'
        if cdf_vals[1]<150.:
            category = 'B'
            if cdf_vals[6]<100.:
                category = 'A'

        category_list.append(category)

    return _gpd.GeoDataFrame(
        data={'Longitude': sample_groups.Longitude.median(),
              'Latitude': sample_groups.Latitude.median(),
              sample_key: list(sample_groups.groups.keys()),
              depositional_age_key: sample_groups[depositional_age_key].median(),
              'TectonicClass': category_list},
        geometry=_gpd.points_from_xy(sample_groups.Longitude.median(), 
                                     sample_groups.Latitude.median()), crs=4326)


def tectonic_fingerprint(SedimentaryZircons, 
                         sample_key='Ref-Sample Key',
                         grain_age_key='Non_Iter_Age_Ma',
                         depositional_age_key='Est_Depos_Age_Ma',
                         chi_square_age_bins=None):
    """
    Given a set of zircon spectra, return values useful to fingerprint the tectonic
    setting following the work of Barham et al (2022, EPSL)
    """

    sample_groups = SedimentaryZircons.groupby(by=sample_key)

    chi_square_list = []
    percentile_list = []

    for sample_group in sample_groups:
        grain_ages = sample_group[1][grain_age_key]
        
        # define the cdf
        if not _np.any(_np.isfinite(grain_ages)):
            chi_square_list.append(_np.nan)
            percentile_list.append(_np.nan)
        if len(grain_ages)<2:
            chi_square_list.append(_np.nan)
            percentile_list.append(_np.nan)
        else:
            dst = _np.sort(grain_ages)
            xtmp = _np.linspace(0,1,len(dst))
            
            # sample cdf at 10% and 50%
            cdf_vals = _np.interp([0.1,0.5],xtmp,dst)

            if chi_square_age_bins is None:
                chi_square_age_bins = len(dst)

            chi_square_list.append(chi_square(dst, chi_square_age_bins))
            percentile_list.append(cdf_vals[1]-cdf_vals[0])

    return _gpd.GeoDataFrame(
        data={'Longitude': sample_groups.Longitude.median(),
              'Latitude': sample_groups.Latitude.median(),
              sample_key: list(sample_groups.groups.keys()),
              depositional_age_key: sample_groups[depositional_age_key].median(),
              'chi_square': chi_square_list,
              'percentile': percentile_list},
        geometry=_gpd.points_from_xy(sample_groups.Longitude.median(), 
                                     sample_groups.Latitude.median()), crs=4326)


def chi_square(dst, num_age_bins=100):

    actual_bin_counts = _np.histogram(dst, bins=_np.linspace(_np.min(dst),4501,num_age_bins))

    expected_bin_counts = len(dst)/len(actual_bin_counts[0])
    
    # chi square test, scaled by division by number of bins
    return _np.sum(((actual_bin_counts[0] - expected_bin_counts)**2) / expected_bin_counts) / len(actual_bin_counts[0])

