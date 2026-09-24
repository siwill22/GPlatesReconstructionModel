"""Unit tests for the small GPML point reader behind the Tappe 2018 and Johansson 2018 loaders.

No network: the fixture below is a two-feature GPML file written out to tmp_path, matching the
structure of the published files (one feature with a shapefile attribute table, one without).

The test that matters is the coordinate order. GPML writes gml:pos as latitude then longitude,
and a reader that takes them as (x, y) still produces a full, plausible-looking point cloud --
just in the wrong places. Nothing downstream would notice.
"""
import numpy as np
import pandas as pd
import pytest

from gprm.datasets._gpml import read_gpml_points

# Igwisi Hills, Tanzania (lat -4.8667, lon 31.9167) and Mendeleev Ridge, Arctic Ocean
# (lat 79.56, lon 177.48): two real points from the published files, both unambiguous about
# which way round the coordinates go.
SAMPLE_GPML = '''<?xml version="1.0" encoding="UTF-8"?>
<gpml:FeatureCollection xmlns:gpml="http://www.gplates.org/gplates"
                        xmlns:gml="http://www.opengis.net/gml"
                        gpml:version="1.6.0318">
    <gml:featureMember>
        <gpml:UnclassifiedFeature>
            <gpml:identity>GPlates-923df0d6-225f-4d12-8efd-8dbbd4d587ee</gpml:identity>
            <gpml:reconstructionPlateId>
                <gpml:ConstantValue>
                    <gpml:value>701</gpml:value>
                    <gpml:valueType>gpml:plateId</gpml:valueType>
                </gpml:ConstantValue>
            </gpml:reconstructionPlateId>
            <gpml:shapefileAttributes>
                <gpml:KeyValueDictionary>
                    <gpml:element>
                        <gpml:KeyValueDictionaryElement>
                            <gpml:key>Kimberlite</gpml:key>
                            <gpml:valueType>xsi:string</gpml:valueType>
                            <gpml:value>Igwisi Hills</gpml:value>
                        </gpml:KeyValueDictionaryElement>
                    </gpml:element>
                    <gpml:element>
                        <gpml:KeyValueDictionaryElement>
                            <gpml:key>Age (recom</gpml:key>
                            <gpml:valueType>xsi:double</gpml:valueType>
                            <gpml:value>0</gpml:value>
                        </gpml:KeyValueDictionaryElement>
                    </gpml:element>
                </gpml:KeyValueDictionary>
            </gpml:shapefileAttributes>
            <gpml:unclassifiedGeometry>
                <gpml:ConstantValue>
                    <gpml:value>
                        <gml:Point><gml:pos>-4.8667 31.9167</gml:pos></gml:Point>
                    </gpml:value>
                    <gpml:valueType>gml:Point</gpml:valueType>
                </gpml:ConstantValue>
            </gpml:unclassifiedGeometry>
            <gml:name>Igwisi Hills</gml:name>
            <gml:description>U-Th-Pb</gml:description>
            <gml:validTime>
                <gml:TimePeriod>
                    <gml:begin><gml:TimeInstant><gml:timePosition>10</gml:timePosition></gml:TimeInstant></gml:begin>
                    <gml:end><gml:TimeInstant><gml:timePosition>-10</gml:timePosition></gml:TimeInstant></gml:end>
                </gml:TimePeriod>
            </gml:validTime>
        </gpml:UnclassifiedFeature>
    </gml:featureMember>
    <gml:featureMember>
        <gpml:UnclassifiedFeature>
            <gpml:identity>GPlates-6d03a603-6a7e-49f6-80b2-1b76162cd6c7</gpml:identity>
            <gpml:unclassifiedGeometry>
                <gpml:ConstantValue>
                    <gpml:value>
                        <gml:Point><gml:pos>79.559953910412972 177.48025476520445</gml:pos></gml:Point>
                    </gpml:value>
                    <gpml:valueType>gml:Point</gpml:valueType>
                </gpml:ConstantValue>
            </gpml:unclassifiedGeometry>
            <gml:validTime>
                <gml:TimePeriod>
                    <gml:begin><gml:TimeInstant><gml:timePosition>120</gml:timePosition></gml:TimeInstant></gml:begin>
                    <gml:end><gml:TimeInstant><gml:timePosition>0</gml:timePosition></gml:TimeInstant></gml:end>
                </gml:TimePeriod>
            </gml:validTime>
            <gml:name>Mendeleyev Ridge</gml:name>
            <gpml:reconstructionPlateId>
                <gpml:ConstantValue>
                    <gpml:value>42100</gpml:value>
                    <gpml:valueType>gpml:plateId</gpml:valueType>
                </gpml:ConstantValue>
            </gpml:reconstructionPlateId>
        </gpml:UnclassifiedFeature>
    </gml:featureMember>
</gpml:FeatureCollection>
'''


def write_sample(tmp_path):
    path = tmp_path / 'sample.gpml'
    path.write_text(SAMPLE_GPML)
    return path


def test_gml_pos_is_read_as_lat_lon_not_lon_lat(tmp_path):
    """The whole point of this file. Igwisi Hills is in Tanzania, not in the Indian Ocean."""
    gdf = read_gpml_points(write_sample(tmp_path))

    igwisi = gdf[gdf.NAME == 'Igwisi Hills'].iloc[0]
    assert np.isclose(igwisi.geometry.x, 31.9167)
    assert np.isclose(igwisi.geometry.y, -4.8667)

    mendeleyev = gdf[gdf.NAME == 'Mendeleyev Ridge'].iloc[0]
    assert np.isclose(mendeleyev.geometry.x, 177.48025476520445)
    assert np.isclose(mendeleyev.geometry.y, 79.559953910412972)

    assert gdf.crs == 'EPSG:4326'
    assert gdf.geometry.y.between(-90, 90).all()


def test_gpml_properties_get_the_column_names_gplates_exports_them_under(tmp_path):
    gdf = read_gpml_points(write_sample(tmp_path))

    assert gdf.loc[0, 'NAME'] == 'Igwisi Hills'
    assert gdf.loc[0, 'DESCR'] == 'U-Th-Pb'
    assert gdf.loc[0, 'FEATURE_ID'] == 'GPlates-923df0d6-225f-4d12-8efd-8dbbd4d587ee'
    assert (gdf.GPGIM_TYPE == 'gpml:UnclassifiedFeature').all()
    for column in ('PLATEID2', 'L_PLATE', 'R_PLATE', 'RECON_METH', 'SPREAD_ASY', 'IMPORT_AGE'):
        assert gdf[column].isna().all()
    # Nothing is invented: no Longitude/Latitude columns unless the file carries them
    assert 'Longitude' not in gdf.columns


def test_unbounded_validity_is_999_and_minus_999_as_gplates_exports_it(tmp_path):
    path = tmp_path / 'unbounded.gpml'
    path.write_text(SAMPLE_GPML
                    .replace('<gml:timePosition>120</gml:timePosition>',
                             '<gml:timePosition>http://gplates.org/times/distantPast</gml:timePosition>')
                    .replace('<gml:timePosition>-10</gml:timePosition>',
                             '<gml:timePosition>http://gplates.org/times/distantFuture</gml:timePosition>'))

    gdf = read_gpml_points(path)
    assert gdf.loc[0, 'TOAGE'] == -999.
    assert gdf.loc[1, 'FROMAGE'] == 999.


def test_plate_ids_and_valid_times_are_read(tmp_path):
    gdf = read_gpml_points(write_sample(tmp_path))

    assert list(gdf.PLATEID1) == [701, 42100]
    assert list(gdf.FROMAGE) == [10.0, 120.0]
    assert list(gdf.TOAGE) == [-10.0, 0.0]


def test_shapefile_attributes_keep_their_source_names_and_types(tmp_path):
    """Including the 10-character truncation the shapefile format imposed on 'Age (recommended)'
    -- expanding it here would invent a field name the source does not use."""
    gdf = read_gpml_points(write_sample(tmp_path))

    assert gdf.loc[0, 'Kimberlite'] == 'Igwisi Hills'
    assert gdf.loc[0, 'Age (recom'] == 0.0
    assert isinstance(gdf.loc[0, 'Age (recom'], float)
    # The second feature has no attribute table at all, so its columns fill with NaN
    assert pd.isna(gdf.loc[1, 'Kimberlite'])
    assert pd.isna(gdf.loc[1, 'Age (recom'])


def test_features_without_a_point_geometry_are_skipped(tmp_path):
    path = tmp_path / 'nogeom.gpml'
    path.write_text(SAMPLE_GPML.replace(
        '<gml:Point><gml:pos>79.559953910412972 177.48025476520445</gml:pos></gml:Point>', ''))

    gdf = read_gpml_points(path)
    assert len(gdf) == 1
    assert gdf.iloc[0].NAME == 'Igwisi Hills'


def test_matches_what_pygplates_exports_for_the_same_file(tmp_path):
    """The contract: this reader exists only so the loaders need not import pygplates, and must
    give the same table pygplates would. Compared on every column pygplates fills in.
    (pygplates' writer also emits an empty 'TYPE' column for some files; that is not in the
    GPML at all, so it is not reproduced.)"""
    pygplates = pytest.importorskip('pygplates')
    from gprm.utils.create_gpml import gpml2gdf

    path = write_sample(tmp_path)
    ours = read_gpml_points(path)
    theirs = gpml2gdf(pygplates.FeatureCollection(str(path)))

    _assert_same_table(ours, theirs)


@pytest.mark.parametrize('fname', ['T18_centroids_M21_plateIDs.gpml',
                                   'J18_centroids_M21_plateIDs.gpml'])
def test_matches_pygplates_on_the_published_files(fname):
    """The same check on the real Tappe 2018 / Johansson 2018 files, when they are in the cache
    (run the network tests once to fetch them)."""
    pygplates = pytest.importorskip('pygplates')
    from gprm.datasets import cache_path
    from gprm.utils.create_gpml import gpml2gdf

    path = cache_path(fname)
    if not path.exists():
        pytest.skip('{} is not in the gprm cache'.format(fname))

    _assert_same_table(read_gpml_points(path), gpml2gdf(pygplates.FeatureCollection(str(path))))


def _assert_same_table(ours, theirs):
    filled = [c for c in theirs.columns
              if c != 'geometry' and theirs[c].replace('', None).notna().any()]
    missing = set(filled) - set(ours.columns)
    assert not missing, 'columns pygplates writes but this reader does not: {}'.format(missing)

    ours = ours.sort_values('FEATURE_ID').reset_index(drop=True)
    theirs = theirs.sort_values('FEATURE_ID').reset_index(drop=True)
    assert len(ours) == len(theirs)
    for column in filled:
        a, b = ours[column], theirs[column]
        if pd.api.types.is_numeric_dtype(b):
            assert np.allclose(pd.to_numeric(a), b, equal_nan=True), column
        else:
            assert (a.astype(str) == b.astype(str)).all(), column
    assert np.allclose(ours.geometry.x, theirs.geometry.x)
    assert np.allclose(ours.geometry.y, theirs.geometry.y)
