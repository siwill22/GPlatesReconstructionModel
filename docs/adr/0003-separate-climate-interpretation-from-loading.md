# Keep the palaeoclimate interpretation out of the dataset loader

`PaleoLithology()` decodes `LithCode` into a human-readable `Indicator` name, which is a
lookup, but does not assign a climate meaning. The Humid/Arid/Glacial grouping lives in a
separate `paleolithology_climate_mapping()` function, mirroring the existing
`pbdb()` / `pbdb_elevation_mapping()` split in the same module. Climate attribution is a
scientific judgement and should be something a user opts into and can cite, not something
that silently arrives with the data.

## Considered Options

We adopt the three-group scheme of Cao et al. (2018, *Geol. Mag.*) — Humid = coals,
Arid = evaporites, Glacial = tillites + dropstones + glendonites — and leave the remaining
codes unclassified. A fuller scheme covering all 16 codes was rejected: that paper argues
explicitly that palms, mangroves and crocodilians suffer sampling bias and that laterites
and oolitic ironstones have samples too small to be reliable latitude indicators, so
classifying them would contradict the citation we are leaning on.

## Consequences

About 24% of points (2093 of 8698) carry no climate group. This is intended — an unlabelled
point is an honest statement that the compilation does not support a confident attribution,
and users who disagree can apply their own mapping to the `Indicator` column.

Note that the compilation's own coding is imperfect and this is not corrected: `M` conflates
Mangroves (22) with Lateritic Manganese (9), and the free-text `Lithology` sometimes
disagrees with `LithCode` (54 records coded `K` describe bauxites; `T` and `D` between them
hold 15 calcretes). The codes are preserved as published.
