"""The OneDZ loader, run against a small zip laid out like the real Zenodo archive.

The real archive is 161 MB, so these tests stand in a few rows that carry the quirks the
loader has to deal with: a UTF-8 BOM, data split across numbered parts, misspelled Lu-Hf
headers, text among the numbers, and rows whose values have shifted into the coordinate columns.
"""
import zipfile

import numpy as np
import pytest


UPB_HEADER = ('Lead_Author,Year,Published_Sample_ID,Latitude,Longitude,'
              'Max_Depos_Age_Ma,Est_Depos_Age_Ma,Min_Depos_Age_Ma,'
              'Grain,Spot_Location,Spot_diam,Pb206U238_age,Pb207U235_age,Pb207Pb206_age,Best_age\n')
LUHF_HEADER = ('Lead_Author,Year,Web Link,Published Sample_ID,Latitud,Logitud,'
               'Max_Depos_Age_Ma,Est_Depos_Age_Ma,Min_Depos_Age_Ma,'
               'Grain,Spot_Location,Spot_diam,Upb_Age,Upb_Age_two_sigma,TDM1_Ma,TDM2_Ma\n')


@pytest.fixture
def onedz_zip(tmp_path, monkeypatch):
    path = tmp_path / 'onedz_datasets_csv.zip'
    with zipfile.ZipFile(path, 'w') as zf:
        # part 02 is written first, to check the loader sorts parts rather than trusting zip order
        zf.writestr('Experts_checked_UPb_split_parts/expert_upb_part_02.csv',
                    '﻿' + UPB_HEADER
                    + '"Smith, A. et al.",2019.0,S2,-30.5,150.2,260,Early Miocene,240,'
                      'S2.1,core,20/25,300,305,−241,nan\n')
        zf.writestr('Experts_checked_UPb_split_parts/expert_upb_part_01.csv',
                    '﻿' + UPB_HEADER
                    + '"Rosenbaum, G. et al.",2020,GI22,-25.7,152.3,277,263,257,'
                      'GI22.001,core,25.0,272.0,280,310,280\n')
        zf.writestr('Total_UPb_split_parts/zircon_upb_part_01.csv',
                    '﻿' + UPB_HEADER + 'Dou Yu Xin,2025.0,22XZ,,,,nan,,nan,nan,,1.4,,,\n')
        zf.writestr('Total_LuHf_split_parts/zircon_luhf_part_01.csv',
                    '﻿' + LUHF_HEADER
                    + '"Jiang, X et al.",2015,10.1016/x,12HN-87,19.4,109.6,145,105.5,66,'
                      '12HN-87@1,,20,239.6,5,900,1474\n'
                    + 'Yang DB et al,2010,,YD1,117.531,Jurassic,,200,,YD1@1,,20,210,4,800,1200\n'
                    + 'Wu et al,2012,,SH1,Shuanghe,30.16,,150,,SH1@1,,20,180,6,700,1100\n')
        zf.writestr('README.md', '# not a csv\n')

    import gprm.datasets.Zircons as Zircons
    monkeypatch.setattr(Zircons, '_retrieve', lambda **kwargs: str(path))
    return Zircons


def test_expert_parts_are_concatenated_in_part_order(onedz_zip):
    gdf = onedz_zip.load_OneDZ()

    assert list(gdf.Published_Sample_ID) == ['GI22', 'S2']
    assert gdf.crs.to_epsg() == 4326
    assert gdf.geometry.iloc[0].x == pytest.approx(152.3)
    assert gdf.geometry.iloc[0].y == pytest.approx(-25.7)


def test_sample_key_identifies_samples_and_year_has_no_trailing_decimal(onedz_zip):
    gdf = onedz_zip.load_OneDZ('expert')

    assert list(gdf.Sample_Key) == ['Rosenbaum, G. et al.|2020|GI22', 'Smith, A. et al.|2019|S2']
    assert list(gdf.Year) == [2020, 2019]                    # the source column, as read


def test_source_columns_are_kept_and_cleaned_numbers_go_in_aliases(onedz_zip):
    """The source keeps its text ('Early Miocene', a unicode minus); the aliases are numbers,
    and a number sitting beside text in the same column survives the cleaning."""
    gdf = onedz_zip.load_OneDZ('expert')

    assert gdf.Est_Depos_Age_Ma.iloc[1] == 'Early Miocene'
    assert gdf.Est_Depos_Age_Ma_numeric.iloc[0] == 263
    assert np.isnan(gdf.Est_Depos_Age_Ma_numeric.iloc[1])

    assert gdf.Pb207Pb206_age.iloc[1] == '−241'
    assert list(gdf['207Pb_206Pb_Age_Ma']) == [310, -241]
    assert list(gdf['206Pb_238U_Age_Ma']) == [272, 300]
    assert gdf.Best_Age_Ma.iloc[0] == 280 and np.isnan(gdf.Best_Age_Ma.iloc[1])

    assert list(gdf.Spot_diam.astype(str)) == ['25.0', '20/25']   # ranges are not numbers


def test_luhf_gets_the_upb_names_alongside_its_own(onedz_zip):
    gdf = onedz_zip.load_OneDZ('LuHf')

    for source, alias in [('Latitud', 'Latitude'), ('Logitud', 'Longitude'),
                          ('Published Sample_ID', 'Published_Sample_ID'), ('Web Link', 'Web_Link')]:
        assert source in gdf.columns and alias in gdf.columns
    assert list(gdf.Published_Sample_ID) == list(gdf['Published Sample_ID'])
    assert gdf.Latitude.iloc[0] == pytest.approx(19.4)


def test_shifted_rows_are_kept_with_missing_geometry(onedz_zip):
    gdf = onedz_zip.load_OneDZ('LuHf')

    assert len(gdf) == 3
    assert gdf.Latitud.iloc[2] == 'Shuanghe'                           # source untouched
    assert np.isnan(gdf.Latitude.iloc[1]) and np.isnan(gdf.Longitude.iloc[1])  # lat 117
    assert np.isnan(gdf.Latitude.iloc[2])                                      # a place name
    assert list(gdf.geometry.isna()) == [False, True, True]


def test_total_includes_rows_without_coordinates(onedz_zip):
    gdf = onedz_zip.load_OneDZ('total')

    assert len(gdf) == 1
    assert gdf.geometry.isna().all()


@pytest.mark.parametrize('dataset', ['expert', 'total', 'LuHf'])
def test_age_description_is_stamped_and_names_only_real_columns(onedz_zip, dataset):
    gdf = onedz_zip.load_OneDZ(dataset)

    description = gdf.attrs['gprm_age']
    assert description['dataset'] == 'Zircons.load_OneDZ:' + dataset
    named = [description['field'], description['source_field'], description['uncertainty']]
    named += list(description['other_ages'])
    assert [c for c in named if c is not None and c not in gdf.columns] == []


def test_unknown_dataset_names_the_valid_options(onedz_zip):
    with pytest.raises(ValueError, match="'expert', 'total' and 'LuHf'"):
        onedz_zip.load_OneDZ('Expert')
