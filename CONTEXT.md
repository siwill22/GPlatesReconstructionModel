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

**Anchor Plate**:
The plate id held fixed when reconstructing (`anchor_plate_id`, default 0). Which Reference
Frame it gives depends on the model; the plates worth anchoring on, and what each means, are
recorded per model in `gprm/datasets/_frames.py` and carried as `model.reference_frames`.
Intermediate plates in a rotation chain (a true polar wander correction, a longitude shift, a
zero-rotation placeholder) are anchorable but are not Reference Frames and are not listed.
_Avoid_: fixed plate (which means the parent plate of one rotation pole)

**Reference Frame**:
What an Anchor Plate is fixed to — one of two things. A **mantle** frame is fixed to the deep
mantle, however it was built: hotspots, palaeomagnetism corrected for true polar wander,
optimisation, no-net-rotation, or Pacific hotspots. A **spin axis** frame is palaeomagnetic,
without a true polar wander correction. Classification is by target, not method: a frame that
uses moving hotspots when young and TPW-corrected palaeomagnetism when old is one mantle frame.
_Avoid_: hybrid frame (as a category), absolute frame (without saying which), palaeomagnetic
frame (for a TPW-corrected one)

**Reconstruction Age**:
The single age a whole time interval is reconstructed to — the rounded mid-age of that
interval, not an age measured on any individual sample.
_Avoid_: reconstruction time (which means the age passed to a reconstruct call)

**Validity Window**:
The `FROMAGE`–`TOAGE` span over which a feature exists. A feature is selected at age *t*
when `TOAGE <= t <= FROMAGE`. Ages are in Ma, increasing into the past. It is **not** the
sample's age, even though published datasets often fill it from one: in different datasets
`FROMAGE` holds the age, the age plus its error, the start of a stratigraphic interval, or
nothing at all. Open-ended windows are written 999 (distant past) and -999 (distant future),
as GPlates exports them.
_Avoid_: age range, time range, valid time

### Sample ages

**Sample Age**:
The age measured or estimated for one sample (a pipe, a province, a zircon grain), in Ma.
Which column holds it varies by dataset and is recorded in its Age Description; it is never
assumed to be `FROMAGE`.
_Avoid_: age (on its own, when the column or dataset is ambiguous), FROMAGE

**Age Description**:
The per-dataset record, kept in `gprm/datasets/_ages.py` and looked up with
`gprm.datasets.age_description()`, of which column holds the Sample Age, what event it dates,
its uncertainty or range, what the Validity Window holds, and any other age-like columns.
Every loader attaches it to the table it returns as `attrs['gprm_age']`; functions that need a
Sample Age read it from there unless told a column explicitly.
_Avoid_: age metadata, age stamp (the attached copy is "the stamp" only in code comments)

**Emplacement Age**:
The event label for the igneous point catalogues — kimberlites, LIPs, carbonatites, seamounts —
whether the body is intrusive or extrusive. Chosen over "eruption age", which is wrong for
carbonatite intrusions and LIP sills and dykes.
_Avoid_: eruption age

**Source Column / Alias**:
A loader returns every column under the name its source gives it (the Source Column), so the
table matches what another program shows for the same file. Where gprm wants its own name —
`Longitude`/`Latitude`, or a short isotope-age name — it adds a copy under that name (an Alias).
Cleaning (wrapping longitudes, fixing decimal commas, blanking a "no age" code) is applied to
the Alias or to a derived column, never to the Source Column.
_Avoid_: renamed column

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
