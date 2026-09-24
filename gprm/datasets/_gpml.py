"""Read point features out of a GPML file into a GeoDataFrame, without pygplates.

pygplates is the canonical GPML reader, but the dataset loaders in this subpackage are
contractually importable without it (see ``tests/test_import_contract.py``), and the files
this is used for are frozen published supplements whose structure cannot drift. So for those,
a small targeted parser is enough, and it keeps Rocks and Seafloor dependency-free.

The output is meant to be indistinguishable from what GPlates itself exports, i.e. from
``gprm.utils.create_gpml.gpml2gdf`` on the same file: the same column names (``NAME``,
``DESCR``, ``PLATEID1``, ``FROMAGE``, ``TOAGE``, ``FEATURE_ID``, ``GPGIM_TYPE`` and the empty
GPlates standard columns), and GPlates' shapefile convention of 999 / -999 for a validity that
runs to the distant past / future. ``tests/test_gpml_reader.py`` checks this against pygplates.

The one thing worth being careful about: GPML writes ``gml:pos`` as **latitude then
longitude**, the opposite order to the (x, y) that geopandas wants. Getting that backwards
produces a plausible-looking global point cloud in entirely the wrong places, so the swap is
done in exactly one place below and is covered by the same tests.

MIT License

Copyright (c) 2017-2021 Simon Williams
"""

import xml.etree.ElementTree as _ET

import geopandas as _gpd
import numpy as _np
import pandas as _pd

_GML = '{http://www.opengis.net/gml}'
_GPML = '{http://www.gplates.org/gplates}'

# GPML writes these in place of a number when a feature's validity is unbounded. GPlates
# exports them as 999 and -999, and so does this reader.
_UNBOUNDED = {
    'http://gplates.org/times/distantPast': 999.,
    'http://gplates.org/times/distantFuture': -999.,
}

# Attribute columns GPlates always writes on export, empty unless the feature sets them.
_EMPTY_STANDARD_COLUMNS = ('PLATEID2', 'L_PLATE', 'R_PLATE', 'RECON_METH', 'SPREAD_ASY',
                           'IMPORT_AGE')


def _time(element):
    if element is None or element.text is None:
        return _np.nan
    text = element.text.strip()
    if text in _UNBOUNDED:
        return _UNBOUNDED[text]
    try:
        return float(text)
    except ValueError:
        return _np.nan


def _shapefile_attributes(feature):
    """The key/value dictionary GPlates carries over from a shapefile's attribute table.

    Keys are whatever the original shapefile had, including the 10-character truncation that
    the shapefile format imposes; they are returned verbatim rather than guessed at.
    """
    attributes = {}
    dictionary = feature.find('{}shapefileAttributes/{}KeyValueDictionary'.format(_GPML, _GPML))
    if dictionary is None:
        return attributes
    for element in dictionary.iterfind('{}element/{}KeyValueDictionaryElement'.format(_GPML, _GPML)):
        key = element.findtext('{}key'.format(_GPML))
        value = element.findtext('{}value'.format(_GPML))
        value_type = (element.findtext('{}valueType'.format(_GPML)) or '').lower()
        if key is None:
            continue
        if value_type.endswith(('double', 'float', 'int', 'integer')) and value is not None:
            try:
                value = float(value)
            except ValueError:
                pass
        attributes[key] = value
    return attributes


def read_gpml_points(fname):
    """Load the point features of a GPML file as a GeoDataFrame in lon/lat (EPSG:4326).

    Columns are those GPlates writes when it exports the same file: each feature's shapefile
    attributes (if any, under their original names), then ``PLATEID1``, ``FROMAGE``, ``TOAGE``
    (the validity window, in Ma; 999 / -999 where it is unbounded), ``NAME``, ``DESCR``,
    ``FEATURE_ID``, ``GPGIM_TYPE``, and the GPlates standard columns that are empty here.

    Features with no point geometry -- lines, polygons, or unclassified features holding
    neither -- are skipped.

    :param fname: path to a .gpml file.
    :returns: geopandas.GeoDataFrame
    """
    root = _ET.parse(str(fname)).getroot()

    records = []
    for member in root.iterfind('{}featureMember'.format(_GML)):
        for feature in member:
            position = feature.find('.//{}Point/{}pos'.format(_GML, _GML))
            if position is None or position.text is None:
                continue
            # GPML orders gml:pos as latitude then longitude. Do not swap these.
            latitude, longitude = (float(value) for value in position.text.split())

            record = _shapefile_attributes(feature)

            plate_id = feature.findtext(
                '{}reconstructionPlateId/{}ConstantValue/{}value'.format(_GPML, _GPML, _GPML))
            record['PLATEID1'] = int(plate_id) if plate_id is not None else 0

            period = feature.find('{}validTime/{}TimePeriod'.format(_GML, _GML))
            if period is None:
                # GPlates treats a feature with no validTime as valid for all time
                record['FROMAGE'], record['TOAGE'] = 999., -999.
            else:
                record['FROMAGE'] = _time(period.find(
                    '{}begin/{}TimeInstant/{}timePosition'.format(_GML, _GML, _GML)))
                record['TOAGE'] = _time(period.find(
                    '{}end/{}TimeInstant/{}timePosition'.format(_GML, _GML, _GML)))

            record['NAME'] = feature.findtext('{}name'.format(_GML))
            record['DESCR'] = feature.findtext('{}description'.format(_GML))
            record['FEATURE_ID'] = feature.findtext('{}identity'.format(_GPML))
            record['GPGIM_TYPE'] = 'gpml:' + feature.tag.split('}', 1)[-1]
            for column in _EMPTY_STANDARD_COLUMNS:
                record.setdefault(column, None)

            record['_lon'], record['_lat'] = longitude, latitude
            records.append(record)

    df = _pd.DataFrame.from_records(records)
    if df.empty:
        return _gpd.GeoDataFrame(df, geometry=_gpd.GeoSeries([], dtype='geometry'), crs=4326)

    geometry = _gpd.points_from_xy(df.pop('_lon'), df.pop('_lat'))
    return _gpd.GeoDataFrame(df, geometry=geometry, crs=4326)
