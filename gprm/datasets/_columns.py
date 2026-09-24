"""Add gprm's conventional column names alongside a dataset's own, never in place of them.

A loader returns every column under the name the source gives it, so that the table looks the
same in Python as it does when the file is opened in any other program. Where gprm wants its own
name for something -- ``Longitude``/``Latitude`` for building geometry, or the short isotope-age
names the Zircons analysis functions expect -- the value is copied into an additional column
(an *alias*). Any cleaning a loader does (wrapping longitudes, fixing decimal commas) is applied
to the alias only, leaving the source column exactly as it was read.

MIT License

Copyright (c) 2017-2021 Simon Williams
"""

import pandas as _pd


def add_aliases(df, aliases):
    """Return ``df`` with extra columns copying existing ones under new names.

    :param df: DataFrame or GeoDataFrame.
    :param aliases: mapping ``{source column: alias}``. Source columns that are not in ``df``
        are skipped, so one mapping can cover several versions of a file (e.g. ``'Lon.'`` in
        some sheets, ``'Lon'`` in others).
    :raises ValueError: if an alias would overwrite a column the source already has -- the
        source's own column always wins, so that is a mistake in the loader, not the data.
    """
    present = {source: alias for source, alias in aliases.items() if source in df.columns}
    targets = list(present.values())
    doubled = sorted(set(a for a in targets if targets.count(a) > 1))
    if doubled:
        raise ValueError('More than one source column maps to alias(es) {}'.format(doubled))
    clash = [alias for alias in targets if alias in df.columns]
    if clash:
        raise ValueError('Alias(es) {} would overwrite columns the source already has'.format(clash))
    if not present:
        return df
    copies = df[list(present)].copy()
    copies.columns = list(present.values())
    # One concat rather than a column at a time: these tables run to a million rows
    return _pd.concat([df, copies], axis=1)
