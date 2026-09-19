# GPlatesReconstructionModel

Tools for building plate tectonic reconstructions with pygplates and for loading the
published palaeo-datasets that get reconstructed against them. This glossary covers
terms whose meaning here is narrower or different from their ordinary geological sense.

## Language

### Reconstruction

**Reconstruction Model**:
A rotation model plus its associated static polygons, coastlines and continent polygons,
treated as one addressable unit.
_Avoid_: rotation file, plate model, rotation model (when the polygons are also meant)

**Static Polygon**:
A present-day polygon carrying a plate ID, used to partition point data onto plates.
Partitioning against static polygons is what assigns a plate ID; it is not itself a
reconstruction.

**Reconstruction Age**:
The single age a whole time interval is reconstructed to — the rounded mid-age of that
interval, not an age measured on any individual sample.
_Avoid_: reconstruction time (which means the age passed to a reconstruct call)

**Validity Window**:
The `FROMAGE`–`TOAGE` span over which a feature exists. A feature is selected at age *t*
when `TOAGE <= t <= FROMAGE`. Ages are in Ma, increasing into the past.
_Avoid_: age range, time range, valid time

### Palaeoclimate indicators

**Palaeoclimate Indicator**:
An observed occurrence whose presence constrains past climate at that site. Includes both
rock types and organisms, so it is broader than "lithology".
_Avoid_: climate proxy, palaeoclimate data

**Lithologic Indicator**:
A palaeoclimate indicator that is a rock, sediment or soil — coals, evaporites, bauxites,
calcretes, tillites, dropstones, glendonites, kaolinites, laterites, oolites.

**Biotic Indicator**:
A palaeoclimate indicator that is an organism or trace fossil — crocodilians, palms,
mangroves, lungfish burrows. Present in the Boucot compilation despite its "lithology"
name, and excluded from the Cao et al. (2018) climate grouping because of sampling bias.

**Indicator Code**:
The one- or two-letter `LithCode` identifying an indicator type in the Boucot compilation
(`C`, `E`, `B`, `K`, `CA`, `T`, `CR`, `L`, `D`, `PA`, `M`, `G`, `O`, `I`, `LF`, `H`).
Distinct from the free-text `Lithology` description, which is uncontrolled vocabulary and
occasionally disagrees with the code.
_Avoid_: lithology code, lithology type

**Climate Group**:
The three-way grouping of indicator codes used by Cao et al. (2018) — Humid (coals), Arid
(evaporites), Glacial (tillites, dropstones, glendonites). An interpretation layered onto
the data, not a property of it, and deliberately leaves the remaining codes unclassified.
_Avoid_: climate class, climate zone, climate type

**Time Slice**:
One of the 28 intervals the Boucot atlas is drawn on, from Cambrian to Miocene. Each point
belongs to exactly one; the intervals do not overlap.
_Avoid_: time bin, stage, map interval
