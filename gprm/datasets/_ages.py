"""What "age" means in each dataset gprm loads.

Every loader returns its source columns under their source names (see ``_columns.py``). Those
names say little about what an age actually dates: ``FROMAGE`` is an emplacement age in one
dataset, age + error in another, a stratigraphic interval in a third and all zeros in a fourth,
and some tables carry two different kinds of age on every row. This module records, per dataset,
which column holds the sample's age and what event it dates, so that gprm's own functions (and the
Geode viewer) do not have to guess. The evidence for each entry is in ``docs/dataset_ages.md``.

Each loader stamps its description onto the returned table as ``gdf.attrs['gprm_age']``. pandas
keeps ``attrs`` through filtering, column selection, copy, sort and similar, but drops it on
``merge``, on ``concat`` of tables with different attrs, and on any file round trip. So the
registry below is the source of truth, and can always be looked up by name with
:func:`age_description`.

A description has these keys:

``dataset``
    the registry key, e.g. ``'Rocks.Carbonatites'``.
``field``
    the column holding the sample's own age, **in Ma**, or ``None`` where the dataset has no single
    such age (only a range, only the age of something else, or two kinds of age with no reason to
    prefer one). A source column name wherever the source gives one in Ma.
``source_field``
    where ``field`` is a column gprm derived (e.g. Ga converted to Ma), the source column it came
    from; otherwise ``None``.
``event``
    what the age dates -- one of :data:`EVENTS` -- or ``None`` where the source does not say.
``uncertainty``
    column holding the age's uncertainty, or ``None``.
``range``
    ``(min_column, max_column)`` bounding the age, or ``None``.
``validity_window``
    what the GPlates ``FROMAGE``/``TOAGE`` columns hold, if present: ``'age to present'`` (FROMAGE
    is the age), ``'age +/- error'``, ``'interval'`` (a stratigraphic range), or ``'none'`` (absent,
    or present but carrying no age information).
``other_ages``
    ``{column: what it dates}`` for further age-like columns that are *not* the sample's age in
    ``field`` -- the depositional age on a zircon grain row, the seafloor age beneath a seamount,
    a model age.
``note``
    anything a user should know before trusting ``field``.

MIT License

Copyright (c) 2017-2021 Simon Williams
"""

import copy as _copy
import warnings as _warnings

#: What an age can date. Taken from the sources' own wording; 'emplacement' was chosen for all
#: igneous point catalogues (kimberlites, LIPs, carbonatites, seamounts) whether intrusive or not.
EVENTS = ('emplacement', 'crystallisation', 'magmatism', 'metamorphism', 'deposition',
          'mineralisation', 'seafloor formation', 'Hf model age')

VALIDITY_WINDOWS = ('age to present', 'age +/- error', 'interval', 'none')

ATTRS_KEY = 'gprm_age'


def _entry(field=None, event=None, source_field=None, uncertainty=None, range=None,
           validity_window='none', other_ages=None, note=None):
    return dict(field=field, source_field=source_field, event=event, uncertainty=uncertainty,
                range=range, validity_window=validity_window, other_ages=other_ages or {},
                note=note)


# The 2018 source's 207Pb/206Pb column holds decimal-comma strings in places, so these point at
# gprm's cleaned numeric aliases rather than the multi-line source headers.
_ZIRCON_2018_ISOTOPE_AGES = {'206Pb_238U_Age_Ma': 'crystallisation (206Pb/238U)',
                             '207Pb_235U_Age_Ma': 'crystallisation (207Pb/235U)',
                             '207Pb_206Pb_Age_Ma': 'crystallisation (207Pb/206Pb)'}
_ZIRCON_2019_ISOTOPE_AGES = {'206Pb/238U Age (Ma)': 'crystallisation (206Pb/238U)',
                             '207Pb/235U Age (Ma)': 'crystallisation (207Pb/235U)',
                             '207Pb/206Pb Age (Ma)': 'crystallisation (207Pb/206Pb)'}
_NO_BEST_AGE_2018 = ('The 2018 compilation reports each grain\'s age by three isotope systems and '
                     'does not choose one, so neither does gprm: pass the column you want.')

_JOHANSSON_ZERO_AGES = ('Age is FROMAGE, except that an emplacement age of 0 is taken to be an '
                        'error (144 of 2526 polygons, 2 of 185 centroids): Age is NaN there, and '
                        'the loader drops those rows unless keep_unknown_age_samples=True.')

_REGISTRY = {
    # --- Rocks --------------------------------------------------------------------------------
    'Rocks.Geochem': _entry(
        'Age', source_field='age', uncertainty='age_sd', range=('age_min', 'age_max'),
        note='What the age dates is not stated per sample; age_method is mostly "estimated". Age '
             'is the source age except: the 2 values above 4567 Ma and the 3 negatives larger than '
             '10 kyr are NaN; the 119 negatives within 10 kyr of 0 (historical eruptions entered as '
             'negative years) are 0. age_min/age_max are the source columns, uncleaned. The 794 '
             'zeros are kept (478 have age_max > 0, so are likely young rather than missing).'),
    'Rocks.BaseMetalDeposits': _entry(
        'Age', 'mineralisation', source_field='Age (Ga)',
        note='The source table does not itself say what the age dates; it is taken to be the '
             'mineralisation age (a working assumption, not checked against Hoggard et al. 2020). '
             '"ND" in Age (Ga) becomes NaN in Age.'),
    'Rocks.Kimberlites:Faure2010': _entry(
        'Age1_Ma', 'emplacement', uncertainty='Age1_Error_Plus_Minus',
        other_ages={'Age2_Ma': 'emplacement, as free text (e.g. "115-135", "Cretaceous (65-144)")'},
        note='Only 641 of 4287 rows have a numeric Age1_Ma.'),
    'Rocks.Kimberlites:Tappe2018': _entry(
        'Age (recom', 'emplacement', validity_window='age to present',
        note='"Age (recom" is the source\'s recommended age, its name truncated by shapefile export.'),
    'Rocks.Carbonatites': _entry(
        'Age', 'emplacement', source_field='Age_ma', uncertainty='Error_ma',
        validity_window='age +/- error',
        note='FROMAGE is Age_ma + Error_ma and TOAGE is Age_ma - Error_ma, so FROMAGE is not the age. '
             'Age is Age_ma except where Age_ma = 0 with Error_ma = 0 (217 of 593 rows, all but one '
             'without a reference), which means "no age determined": Age is NaN there, and the '
             'loader drops those rows unless keep_unknown_age_samples=True.'),
    'Rocks.Metamorphism': _entry(
        'Age', 'metamorphism', source_field='AGE(Ga)',
        note='"Peak P/P-T" says which point on the P-T path the age and conditions refer to.'),

    # --- Seafloor -----------------------------------------------------------------------------
    'Seafloor.MagneticPicks': _entry(
        'GeeK2007', 'seafloor formation',
        note='Age of the picked chron boundary on the Gee & Kent (2007) timescale.'),
    'Seafloor.PacificSeamountAges:2013': _entry(
        'Average_age(Ma)', 'emplacement', uncertainty='Average_error(Ma)'),
    'Seafloor.PacificSeamountAges:2021': _entry(
        'Age', 'emplacement', uncertainty='Error'),
    'Seafloor.Seamounts:KimWessel': _entry(
        other_ages={'CrustAge': 'seafloor formation (of the crust beneath the seamount)'},
        note='No seamount ages. CrustAge is the age of the underlying seafloor from the AGE 3.2 grid.'),
    'Seafloor.Seamounts:SIO': _entry(note='No ages.'),
    'Seafloor.Seamounts:HillierWatts': _entry(note='No ages.'),
    'Seafloor.LargeIgneousProvinces:Whittaker': _entry(
        'FROMAGE', 'emplacement', validity_window='age to present'),
    'Seafloor.LargeIgneousProvinces:Johansson': _entry(
        'Age', 'emplacement', source_field='FROMAGE', validity_window='age to present',
        note=_JOHANSSON_ZERO_AGES),
    'Seafloor.LargeIgneousProvinces:Johansson_centroids': _entry(
        'Age', 'emplacement', source_field='FROMAGE', validity_window='age to present',
        note=_JOHANSSON_ZERO_AGES),
    'Seafloor.LargeIgneousProvinces:UTIG': _entry(
        note='No ages: FROMAGE and TOAGE are 0 on every feature.'),

    # --- Strat --------------------------------------------------------------------------------
    'Strat.pbdb': _entry(
        event='deposition', range=('min_ma', 'max_ma'),
        note='Only the stratigraphic interval of each collection is known.'),
    'Strat.PaleoReefs': _entry(
        event='deposition',
        note='No numeric age: only chronostratigraphic names in system, series and intervall '
             '(free text, stage- to series-level). Converting them to Ma needs a timescale '
             'chosen by the user.'),
    'Strat.PaleoLithology': _entry(
        event='deposition', range=('TOAGE', 'FROMAGE'), validity_window='interval',
        other_ages={'ReconstructionAge': 'mid-age of the time slice, not a measured age'}),

    # --- Geology ------------------------------------------------------------------------------
    'Geology.GlobalTectonicMap': _entry(
        other_ages={'mag_min_age': 'magmatism (minimum)', 'mag_max_age': 'magmatism (maximum)',
                    'met_min_age': 'metamorphism (minimum)', 'met_max_age': 'metamorphism (maximum)'},
        note='Each province carries separate magmatic and metamorphic age ranges.'),
    'Geology.SurfaceGeology': _entry(
        range=('TOAGE', 'FROMAGE'), validity_window='interval',
        note='Age range of each map unit, derived from era names; units are of mixed rock types.'),

    # --- Zircons ------------------------------------------------------------------------------
    'Zircons.loadDB:2018:samples': _entry(
        'Est. Depos. Age (Ma)', 'deposition',
        other_ages={'Max. Depos. Age (Ma)': 'deposition (maximum)'}),
    'Zircons.loadDB:2018:data': _entry(
        event='crystallisation', other_ages=_ZIRCON_2018_ISOTOPE_AGES, note=_NO_BEST_AGE_2018),
    'Zircons.loadDB:2019': _entry(
        'Best Age (Ma)', 'crystallisation',
        other_ages=dict({'Estimated Dep_Age': 'deposition'}, **_ZIRCON_2019_ISOTOPE_AGES)),
    'Zircons.loadDB:2021': _entry(
        'Non-iter          age                     (Ma)', 'crystallisation',
        other_ages={'Est. Depos. Age': 'deposition', 'Min. Depos. Age': 'deposition (minimum)',
                    'Max. Depos. Age': 'deposition (maximum)'}),
    'Zircons.loadDB:2024': _entry(
        'Non-Iter. Probability age (Ma)', 'crystallisation',
        other_ages={'Est. Stratigraphic Age (Ma)': 'deposition',
                    'Min. Stratigraphic Age (Ma)': 'deposition (minimum)',
                    'Max. Stratigraphic Age (Ma)': 'deposition (maximum)'},
        note='This release has no schema sheet; the grain age is identified by analogy with the '
             'same authors\' 2026 schema, which defines its "U-Pb Non-Iter. Prob. age (Ma)" as '
             '"The U-Pb age of the detrital zircon".'),
    'Zircons.loadDB:2026': _entry(
        'U-Pb Non-Iter. Prob. age (Ma)', 'crystallisation',
        other_ages={'Est. Stratigraphic Age (Ma)': 'deposition',
                    'Min. Stratigraphic Age (Ma)': 'deposition (minimum)',
                    'Max. Stratigraphic Age (Ma)': 'deposition (maximum)',
                    'Strat/Dep  Age (Ma)': 'deposition (of the sample the grain came from)'}),
    'Zircons.load_Hf': _entry(
        'U-Pb    Age  (Ma)', 'crystallisation',
        other_ages={'Est. Depos. Age (Ma)': 'deposition', 'Min. Depos. Age (Ma)': 'deposition (minimum)',
                    'Max. Depos. Age (Ma)': 'deposition (maximum)',
                    'TDM1 (Ma)': 'Hf model age', 'TDM2 (Ma)': 'Hf model age'}),
    'Zircons.get_igneous_samples:2018': _entry(
        event='crystallisation', other_ages=_ZIRCON_2018_ISOTOPE_AGES, note=_NO_BEST_AGE_2018),
    'Zircons.get_igneous_samples:2019': _entry(
        'Best Age (Ma)', 'crystallisation', uncertainty='Uncertainty (2σ precision)',
        other_ages=_ZIRCON_2019_ISOTOPE_AGES),
    'Zircons.get_igneous_samples:2026': _entry(
        'Magm. / crystal age (Ma)', 'crystallisation',
        other_ages={'Map age (Ma)': '50-Myr bin the sample is mapped in, not a measured age'},
        note='"Type of age" gives the U-Pb method (WMA, UI, LI, CA, TuffZirc).'),
    'Zircons.get_mafic_felsic_samples': _entry(
        'Magm. / crystal age (Ma)', 'crystallisation',
        other_ages={'Map age (Ma)': '50-Myr bin the sample is mapped in, not a measured age'},
        note='"Type of age" gives the U-Pb method (WMA, UI, LI, CA, TuffZirc).'),
    'Zircons.get_sedimentary_samples:2018': _entry(
        event='crystallisation',
        other_ages=dict({'Est. Depos. Age (Ma)': 'deposition',
                         'Max. Depos. Age (Ma)': 'deposition (maximum)'}, **_ZIRCON_2018_ISOTOPE_AGES),
        note=_NO_BEST_AGE_2018),
}


def datasets():
    """Names of every dataset with an age description."""
    return sorted(_REGISTRY)


def age_description(dataset):
    """Return what "age" means in a dataset, as a dict (see the module docstring for the keys).

    :param dataset: registry name, e.g. ``'Rocks.Carbonatites'`` or
        ``'Seafloor.LargeIgneousProvinces:Johansson'``. :func:`datasets` lists them.
    """
    if dataset not in _REGISTRY:
        raise KeyError('No age description for {!r}. Known datasets: {}'.format(
            dataset, ', '.join(datasets())))
    description = _copy.deepcopy(_REGISTRY[dataset])
    description['dataset'] = dataset
    return description


def stamp(gdf, dataset):
    """Attach ``dataset``'s age description to ``gdf.attrs['gprm_age']`` and return ``gdf``."""
    gdf.attrs[ATTRS_KEY] = age_description(dataset)
    return gdf


def resolve_age_field(df, age_field=None, default='age', argument='age_field'):
    """Decide which column holds each sample's own age, for a function about to use it.

    In order of precedence:

    1. ``age_field``, if the caller passed one -- an explicit choice always wins.
    2. The ``field`` of the age description a gprm loader stamped on ``df.attrs['gprm_age']``.
       If that description says the dataset has no single sample age (``field`` is None), or
       names a column that is no longer in ``df``, this raises rather than guess.
    3. ``default``, with a warning that it is being assumed. This keeps gprm's long-standing
       convention (e.g. the sample age held in ``FROMAGE``) working for tables that did not come
       from a gprm loader, or whose ``attrs`` were dropped by a merge or a file round trip.

    :param df: the table of samples.
    :param age_field: the column the caller named, or None.
    :param default: the column assumed when there is neither an argument nor a description.
    :param argument: the calling function's name for ``age_field``, used in messages.
    :returns: a column name.
    """
    if age_field is not None:
        return age_field

    description = getattr(df, 'attrs', {}).get(ATTRS_KEY)
    if description is not None:
        field = description.get('field')
        dataset = description.get('dataset', 'this dataset')
        if field is None:
            others = list(description.get('other_ages') or {}) + list(description.get('range') or ())
            raise ValueError(
                "{} has no single sample age to use{}. Pass {}='<column>' to choose one{}.".format(
                    dataset,
                    ' ({})'.format(description['note']) if description.get('note') else '',
                    argument,
                    ' -- its age-like columns are {}'.format(others) if others else ''))
        if field not in df.columns:
            raise ValueError(
                "The age description for {} says the sample age is in column {!r}, which is not in "
                "this table (was it left out with usecols, or dropped?). Pass {}='<column>' to use "
                "another.".format(dataset, field, argument))
        return field

    _warnings.warn(
        "No age description on this table and no {0} given, so {1!r} is assumed to hold each "
        "sample's own age. That is true only if you put it there: in several published datasets "
        "FROMAGE holds something else (age + error, a stratigraphic interval, nothing). Pass "
        "{0}='<column>' to say which column holds the age, or load the data with a "
        "gprm.datasets loader, which records it.".format(argument, default),
        UserWarning, stacklevel=3)
    return default
