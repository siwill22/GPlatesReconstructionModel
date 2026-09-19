"""
Rebuild gprm/Data/boucot_paleolithology.gpkg from the published source tables.

The bundled shapefile (gprm/Data/boucot_paleolithology_combined.shp) had five text
columns destroyed by a numeric type coercion at some point in its history: Lithology,
Formation, LithComm, Period and AgeComm are all 0.0 or null. This script recovers them
by joining the 28 original CSV tables of Boucot, Chen & Scotese (2013) back onto the
shapefile's GPlates attributes.

Usage:
    python tools/build_paleolithology.py /path/to/"Boucot et al. 2013 Lithology Data Tables"

See docs/adr/0001-regenerate-boucot-from-source-tables.md for why this join is shaped
the way it is.
"""
import sys
import glob
import os
import math

import numpy as np
import pandas as pd
import geopandas as gpd

SHAPEFILE = os.path.join(os.path.dirname(__file__), '..', 'gprm', 'Data',
                          'boucot_paleolithology_combined.shp')
OUTPUT = os.path.join(os.path.dirname(__file__), '..', 'gprm', 'Data',
                       'boucot_paleolithology.gpkg')

# Header spellings vary across the 28 source files (e.g. 'LithNumber' vs
# 'LithologyIDNumber', 'AgeComment' vs 'AgeComments', a leading space in one file).
# Anything not listed here is dropped (stray unnamed/blank columns from trailing commas).
CANONICAL_COLUMNS = {
    'lithologycode': 'LithCode',
    'lithologyidnumber': 'LithID',
    'lithnumber': 'LithID',
    'oldidnumber': 'OldID',
    'lat': 'LAT',
    'long': 'LONG',
    'continent': 'Continent',
    'country': 'Country',
    'geogcomments': 'GeogComm',
    'geogcomment': 'GeogComm',
    'lmu': 'LMU',
    'period': 'Period',
    'stage': 'Stage',
    'agecomments': 'AgeComm',
    'agecomment': 'AgeComm',
    'lithology': 'Lithology',
    'formation': 'Formation',
    'lithcomments': 'LithComm',
    'lithcomment': 'LithComm',
    'primaryreference': 'PrimRef',
    'seealso': 'SeeAlso',
}

RECOVERED_TEXT_COLUMNS = ['Lithology', 'Formation', 'LithComm', 'Period', 'AgeComm', 'OldID']

# GPlates reconstruction-specific columns from the shapefile: PLATEID1/PLATEID2/L_PLATE/
# R_PLATE/SPREAD_ASY are dropped rather than kept, because they were assigned by
# partitioning against an unrecorded static polygon set (see ADR 0002) and would silently
# mislead anyone reconstructing against a different plate model. NAME/DESCR/RECON_METH/TYPE
# are empty in every row; GPGIM_TYPE is constant.
SHAPEFILE_COLUMNS = ['LithCode', 'LithID', 'GeogComm', 'Continent', 'Country', 'LMU',
                      'Stage', 'PrimRef', 'SeeAlso', 'FEATURE_ID',
                      'FROMAGE', 'TOAGE', 'geometry']


def load_source_tables(source_dir):
    files = sorted(glob.glob(os.path.join(source_dir, '*.csv')))
    if len(files) != 28:
        raise ValueError(f'expected 28 source CSVs in {source_dir!r}, found {len(files)}')

    frames = []
    for map_no, f in enumerate(files, start=1):
        df = pd.read_csv(f, dtype=str, keep_default_na=False, encoding='cp1252')
        df.columns = [CANONICAL_COLUMNS.get(c.strip().lower(), '__drop_' + c) for c in df.columns]
        df = df[[c for c in df.columns if not c.startswith('__drop_')]]
        df['MapNo'] = map_no
        frames.append(df)

    source = pd.concat(frames, ignore_index=True)
    source['LithCode'] = source['LithCode'].str.strip().str.upper()
    # A single editorial marker row ('DO NOT USE') precedes a stray appended block in the
    # Miocene table; it carries no coordinates and is not a real observation.
    source = source[source['LithCode'] != 'DO NOT USE'].reset_index(drop=True)
    source['LithID'] = pd.to_numeric(source['LithID'], errors='coerce')
    source['LAT_n'] = pd.to_numeric(source['LAT'], errors='coerce')
    source['LONG_n'] = pd.to_numeric(source['LONG'], errors='coerce')

    for col in RECOVERED_TEXT_COLUMNS:
        source[col] = source[col].str.strip().replace('', np.nan)

    source['_srcidx'] = source.index
    return source


def round_half_up(x):
    return math.floor(x + 0.5)


def build(source_dir):
    source = load_source_tables(source_dir)

    gdf = gpd.read_file(SHAPEFILE)
    gdf = gdf[SHAPEFILE_COLUMNS].copy()

    # MapNo isn't in the shapefile; recover it from FROMAGE, whose 28 distinct values
    # correspond 1:1 with the 28 source files (oldest = Map01, youngest = Map28).
    ages_desc = sorted(gdf['FROMAGE'].unique(), reverse=True)
    if len(ages_desc) != 28:
        raise ValueError(f'expected 28 distinct FROMAGE values, found {len(ages_desc)}')
    gdf['MapNo'] = gdf['FROMAGE'].map({age: i + 1 for i, age in enumerate(ages_desc)})

    gdf['_glon'] = gdf.geometry.x
    gdf['_glat'] = gdf.geometry.y
    gdf['_gidx'] = gdf.index

    source_cols = ['MapNo', 'LithCode', 'LithID', '_srcidx', 'LAT_n', 'LONG_n'] + RECOVERED_TEXT_COLUMNS
    merged = gdf[['_gidx', 'MapNo', 'LithCode', 'LithID', '_glon', '_glat']].merge(
        source[source_cols], on=['MapNo', 'LithCode', 'LithID'], how='left')

    # (MapNo, LithCode, LithID) is not always unique in the source (a handful of stray
    # duplicate/appended entries) — break ties by nearest coordinate to the shapefile point.
    merged['_dist'] = np.hypot(merged['_glon'] - merged['LONG_n'], merged['_glat'] - merged['LAT_n'])
    merged['_dist'] = merged['_dist'].fillna(1e9)
    best = merged.sort_values('_dist').drop_duplicates('_gidx', keep='first').sort_values('_gidx')

    n_unmatched = best['_srcidx'].isna().sum()
    if n_unmatched:
        raise ValueError(f'{n_unmatched} shapefile rows had no matching source row')

    gdf = gdf.merge(best[['_gidx'] + RECOVERED_TEXT_COLUMNS], on='_gidx', how='left')
    gdf = gdf.drop(columns=['_glon', '_glat', '_gidx'])

    gdf['ReconstructionAge'] = gdf.apply(
        lambda r: round_half_up((r['FROMAGE'] + r['TOAGE']) / 2), axis=1)

    assert len(gdf) == 8698, f'expected 8698 rows, got {len(gdf)}'
    assert gdf['Lithology'].notna().sum() == 7586, \
        f"expected 7586 non-null Lithology, got {gdf['Lithology'].notna().sum()}"

    gdf = gpd.GeoDataFrame(gdf, geometry='geometry', crs='EPSG:4326')
    gdf.to_file(OUTPUT, driver='GPKG')
    print(f'wrote {OUTPUT} ({len(gdf)} rows)')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    build(sys.argv[1])
