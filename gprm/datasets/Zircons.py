from pooch import os_cache as _os_cache
from pooch import retrieve as _retrieve
from pooch import HTTPDownloader as _HTTPDownloader
from pooch import Unzip as _Unzip
import pandas as _pd
import geopandas as _gpd
import os as _os


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
        df_Data = _pd.read_excel(xls, sheet_name='Data')
        df_SampleDetails = _pd.read_excel(xls, sheet_name='Sample_Details')
        
        # rename some fields for neatness
        df_Data.rename(columns = {u'206Pb /\n238U\nAge\n(Ma)': '206Pb_238U_Age_Ma',
                                u'206Pb /\n238U\n2σ\nPrecis': '206Pb_238U_Precis',
                                u'207Pb /\n235U\nAge\n(Ma)': '207Pb_235U_Age_Ma',
                                u'207Pb /\n235U\n2σ\nPrecis': '207Pb_235U_Precis',
                                u'207Pb /\n206Pb\nAge\n(Ma)': '207Pb_206Pb_Age_Ma',
                                u'207Pb /\n206Pb\n2σ\nPrecis': '207Pb_206Pb_Precis'},
                    inplace = True)
        
        df_SampleDetails.rename(columns = {'Est. Depos. Age (Ma)': 'Est_Depos_Age_Ma'}, inplace = True)
        df_SampleDetails = df_SampleDetails.dropna(subset=['Sample Key'])

        return df_SampleDetails, df_Data

    elif version==2021:

        fname = _retrieve(
            url="https://ars.els-cdn.com/content/image/1-s2.0-S0012825221002464-mmc4.xlsx",
            known_hash="sha256:afdcaa26698ba06f113f29fab311af0ec37f290db117fc85847966c28b78ff09",  
            downloader=_HTTPDownloader(progressbar=True),
            path=_os_cache('gprm'),
        )

        xls = _pd.ExcelFile(fname)
        df = _pd.read_excel(xls, sheet_name='UPb_Data')

        gdf = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

        return gdf


def load_Hf():

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
def get_igneous_samples(df_SampleDetails,df_Data):       
    # igneous samples
    df_IgneousSamples = df_SampleDetails[df_SampleDetails['Class-1 Rock Type'].str.contains('igneous')]
    df = _pd.merge(df_IgneousSamples,df_Data,on='Sample Key')
    IgneousZircons = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return IgneousZircons


def get_sedimentary_samples(df_SampleDetails,df_Data):
    # sedimentary samples
    df_SamplesWithDepositionalAge = df_SampleDetails.dropna(subset=['Est_Depos_Age_Ma'])
    df = _pd.merge(df_SamplesWithDepositionalAge,df_Data,on='Sample Key')
    SedimentaryZircons = _gpd.GeoDataFrame(df, geometry=_gpd.points_from_xy(df.Longitude, df.Latitude), crs=4326)

    return SedimentaryZircons


def get_record_spectrum(df,Sample_ID):
    return df[df.Sample_ID_x==Sample_ID].dropna(subset=['206Pb_238U_Age_Ma'])

