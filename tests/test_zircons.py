"""The OneDZ loader, run against a small zip laid out like the real Zenodo archive.

The real archive is 161 MB, so these tests stand in a few rows that carry the quirks the
loader has to deal with: a UTF-8 BOM, data split across numbered parts, misspelled Lu-Hf
headers, and rows whose values have shifted into the coordinate columns.
"""
import zipfile

import numpy as np
import pytest


UPB_HEADER = ('Lead_Author,Year,Published_Sample_ID,Latitude,Longitude,'
              'Max_Depos_Age_Ma,Est_Depos_Age_Ma,Min_Depos_Age_Ma,'
              'Grain,Spot_Location,Spot_diam,Pb206U238_age,Pb207Pb206_age,Best_age\n')
LUHF_HEADER = ('Lead_Author,Year,Published Sample_ID,Latitud,Logitud,'
               'Max_Depos_Age_Ma,Est_Depos_Age_Ma,Min_Depos_Age_Ma,'
               'Grain,Spot_Location,Spot_diam,Upb_Age,epsilon_Hf_t\n')


@pytest.fixture
def onedz_zip(tmp_path, monkeypatch):
    path = tmp_path / 'onedz_datasets_csv.zip'
    with zipfile.ZipFile(path, 'w') as zf:
        # part 02 is written first, to check the loader sorts parts rather than trusting zip order
        zf.writestr('Experts_checked_UPb_split_parts/expert_upb_part_02.csv',
                    '﻿' + UPB_HEADER + '"Smith, A. et al.",2019.0,S2,-30.5,150.2,260,250,240,S2.1,core,20/25,300,−241,nan\n')
        zf.writestr('Experts_checked_UPb_split_parts/expert_upb_part_01.csv',
                    '﻿' + UPB_HEADER + '"Rosenbaum, G. et al.",2020,GI22,-25.7,152.3,277,263,257,GI22.001,core,25.0,272.0,310,280\n')
        zf.writestr('Total_UPb_split_parts/zircon_upb_part_01.csv',
                    '﻿' + UPB_HEADER + 'Dou Yu Xin,2025.0,22XZ,,,,nan,,nan,nan,,1.4,,\n')
        zf.writestr('Total_LuHf_split_parts/zircon_luhf_part_01.csv',
                    '﻿' + LUHF_HEADER
                    + '"Jiang, X et al.",2015,12HN-87,19.4,109.6,145,105.5,66,12HN-87@1,,20,239.6,-7.98\n'
                    + 'Yang DB et al,2010,YD1,117.531,Jurassic,,200,,YD1@1,,20,210,1.5\n'
                    + 'Wu et al,2012,SH1,Shuanghe,30.16,,150,,SH1@1,,20,abc,2.0\n')
        zf.writestr('README.md', '# not a csv\n')

    import gprm.datasets.Zircons as Zircons
    monkeypatch.setattr(Zircons, '_retrieve', lambda **kwargs: str(path))
    return Zircons


def test_expert_parts_are_concatenated_in_part_order(onedz_zip):
    gdf = onedz_zip.load_OneDZ()

    assert list(gdf.Published_Sample_ID) == ['GI22', 'S2']
    assert gdf.crs.to_epsg() == 4326
    assert gdf.geometry.iloc[0].x == pytest.approx(152.3)


def test_sample_key_identifies_samples_and_year_has_no_trailing_decimal(onedz_zip):
    gdf = onedz_zip.load_OneDZ('expert')

    assert list(gdf.Sample_Key) == ['Rosenbaum, G. et al.|2020|GI22', 'Smith, A. et al.|2019|S2']


def test_luhf_headers_are_renamed_to_match_upb(onedz_zip):
    gdf = onedz_zip.load_OneDZ('LuHf')

    for column in ('Latitude', 'Longitude', 'Published_Sample_ID'):
        assert column in gdf.columns
    for column in ('Latitud', 'Logitud', 'Published Sample_ID'):
        assert column not in gdf.columns


def test_shifted_rows_are_kept_with_missing_coordinates(onedz_zip):
    gdf = onedz_zip.load_OneDZ('LuHf')

    assert len(gdf) == 3
    assert gdf.Latitude.iloc[0] == pytest.approx(19.4)
    assert np.isnan(gdf.Latitude.iloc[1]) and np.isnan(gdf.Longitude.iloc[1])  # lat 117
    assert np.isnan(gdf.Latitude.iloc[2])                                      # a place name
    assert list(gdf.geometry.isna()) == [False, True, True]
    assert np.isnan(gdf.Upb_Age.iloc[2])                    # stray text in a measurement column
    assert gdf.epsilon_Hf_t.dtype.kind == 'f'


def test_unicode_minus_is_read_and_spot_size_ranges_are_kept(onedz_zip):
    gdf = onedz_zip.load_OneDZ('expert')

    assert list(gdf.Pb207Pb206_age) == [310, -241]   # the number beside the text must survive
    assert list(gdf.Est_Depos_Age_Ma) == [263, 250]
    assert list(gdf.Spot_diam.astype(str)) == ['25.0', '20/25']


def test_total_includes_rows_without_coordinates(onedz_zip):
    gdf = onedz_zip.load_OneDZ('total')

    assert len(gdf) == 1
    assert gdf.geometry.isna().all()


def test_unknown_dataset_names_the_valid_options(onedz_zip):
    with pytest.raises(ValueError, match="'expert', 'total' and 'LuHf'"):
        onedz_zip.load_OneDZ('Expert')
