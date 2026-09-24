# What "age" means in each gprm dataset

"Age" means different things in different datasets. Most of these datasets are GPlates-flavoured,
so `FROMAGE` is always there to tempt you, but it holds the age in some datasets, age + error in
another, a stratigraphic interval in others, and all zeros in one. This page explains how gprm
deals with that, and gives the evidence behind each dataset's entry. Terms follow
[CONTEXT.md](../CONTEXT.md) (*Sample Age*, *Age Description*, *Validity Window*,
*Source Column / Alias*).

## How gprm handles it

1. **Loaders keep source column names.** Every column comes back under the name its source gives
   it, so the table matches what you see opening the file in another program. Where gprm wants its
   own name (`Longitude`, `Latitude`, `206Pb_238U_Age_Ma`, …) it adds a copy alongside. Cleaning
   is done on the copy or on a derived column, never on the source column.
2. **Every loader attaches an Age Description** to its result as `gdf.attrs['gprm_age']`. It says
   which column holds the sample's age in Ma (`field`), what event that age dates (`event`), its
   uncertainty or range, what `FROMAGE`/`TOAGE` hold (`validity_window`), and which other columns
   are ages of something else (`other_ages`). The registry is in
   [gprm/datasets/_ages.py](../gprm/datasets/_ages.py), and any entry can be read without loading
   the data:

   ```python
   from gprm.datasets import age_description
   age_description('Rocks.Carbonatites')
   # {'field': 'Age', 'source_field': 'Age_ma', 'event': 'emplacement',
   #  'uncertainty': 'Error_ma', 'validity_window': 'age +/- error', ...}
   ```

3. **gprm functions that need a sample age read it from the description.** This applies to
   `ReconstructionModel.reconstruct_to_time_of_appearance` (for GeoDataFrames),
   `molchan.sample_distance_analysis` and `molchan.space_time_distances`. The order of precedence:
   - a column you pass explicitly always wins;
   - otherwise the description's `field` is used. If the description says the dataset has *no*
     single sample age (`field` is `None`: pbdb, Zircons 2018, …), or names a column that is no
     longer in the table, the function raises instead of guessing;
   - otherwise (a table you built yourself, or one whose `attrs` were dropped) the old default
     (`FROMAGE`, or `age` in `molchan`) is used, **with a warning**.

   pandas keeps `attrs` through filtering, column selection, copy, sort, `reset_index`, `to_crs`
   and `groupby().first()`. It drops them on `merge`, on `concat` of tables with different `attrs`,
   on row-wise `apply`, and on any file round trip. After any of those, pass the column explicitly.

4. **"No age" codes are dataset-specific and handled per loader.** Carbonatites codes "no age
   determined" as `Age_ma = 0` with `Error_ma = 0`. The Johansson LIP catalogues have
   emplacement ages of 0 that are errors (the Cenozoic Tuamotu seamounts among them). In both
   cases those rows are dropped by default; `keep_unknown_age_samples=True` keeps them, with the
   source column untouched and the derived `Age` NaN. In the 2021 Pacific seamount file, by
   contrast, the zero ages belong to active volcanoes and are real.

## Registry entries at a glance

`field` is the column holding the sample age in Ma. `None` means there is no single sample age.

| Dataset (`age_description` key) | `field` | `event` | uncertainty / range | `validity_window` | other ages |
|---|---|---|---|---|---|
| `Rocks.Geochem` | `Age` (from `age`) | not stated | `age_sd`; `age_min`–`age_max` | none | — |
| `Rocks.BaseMetalDeposits` | `Age` (from `Age (Ga)`) | mineralisation (assumed) | — | none | — |
| `Rocks.Kimberlites:Faure2010` | `Age1_Ma` | emplacement | `Age1_Error_Plus_Minus` | none | `Age2_Ma` (free text) |
| `Rocks.Kimberlites:Tappe2018` | `Age (recom` | emplacement | — | age to present | — |
| `Rocks.Carbonatites` | `Age` (from `Age_ma`) | emplacement | `Error_ma` | age +/- error | — |
| `Rocks.Metamorphism` | `Age` (from `AGE(Ga)`) | metamorphism | — | none | — |
| `Seafloor.MagneticPicks` | `GeeK2007` | seafloor formation | — | none | — |
| `Seafloor.PacificSeamountAges:2013` | `Average_age(Ma)` | emplacement | `Average_error(Ma)` | none | — |
| `Seafloor.PacificSeamountAges:2021` | `Age` | emplacement | `Error` | none | — |
| `Seafloor.Seamounts:KimWessel` | None | — | — | none | `CrustAge` (seafloor beneath) |
| `Seafloor.Seamounts:SIO`, `:HillierWatts` | None | — | — | none | — |
| `Seafloor.LargeIgneousProvinces:Whittaker` | `FROMAGE` | emplacement | — | age to present | — |
| `Seafloor.LargeIgneousProvinces:Johansson`, `:Johansson_centroids` | `Age` (from `FROMAGE`) | emplacement | — | age to present | — |
| `Seafloor.LargeIgneousProvinces:UTIG` | None | — | — | none (all zeros) | — |
| `Strat.pbdb` | None | deposition | `min_ma`–`max_ma` | none | — |
| `Strat.PaleoLithology` | None | deposition | `TOAGE`–`FROMAGE` | interval | `ReconstructionAge` |
| `Geology.GlobalTectonicMap` | None | — | — | none | `mag_*_age`, `met_*_age` |
| `Geology.SurfaceGeology` | None | not stated | `TOAGE`–`FROMAGE` | interval | — |
| `Zircons.loadDB:2018:samples` | `Est. Depos. Age (Ma)` | deposition | — | none | `Max. Depos. Age (Ma)` |
| `Zircons.loadDB:2018:data`, `get_igneous_samples:2018`, `get_sedimentary_samples:2018` | None | crystallisation | — | none | three isotope-system ages (+ deposition) |
| `Zircons.loadDB:2019` | `Best Age (Ma)` | crystallisation | — | none | `Estimated Dep_Age`, isotope ages |
| `Zircons.loadDB:2021` | `Non-iter … age … (Ma)` | crystallisation | — | none | `Est./Min./Max. Depos. Age` |
| `Zircons.loadDB:2024` | `Non-Iter. Probability age (Ma)` | crystallisation | — | none | `Est./Min./Max. Stratigraphic Age (Ma)` |
| `Zircons.loadDB:2026` | `U-Pb Non-Iter. Prob. age (Ma)` | crystallisation | — | none | stratigraphic ages, `Strat/Dep  Age (Ma)` |
| `Zircons.load_Hf` | `U-Pb    Age  (Ma)` | crystallisation | — | none | depositional ages, `TDM1`/`TDM2` (Hf model) |
| `Zircons.get_igneous_samples:2019` | `Best Age (Ma)` | crystallisation | `Uncertainty (2σ precision)` | none | isotope ages |
| `Zircons.get_igneous_samples:2026`, `get_mafic_felsic_samples` | `Magm. / crystal age (Ma)` | crystallisation | — | none | `Map age (Ma)` (50-Myr bin) |

The `note` on each entry, which `age_description` returns, carries the caveats: for example,
how Geochem's `Age` was cleaned, that Hoggard's ages are assumed to be mineralisation ages, and
that Zircons 2024 has no schema sheet.

## The evidence

Each loader was run from the local cache (2026-09-23/24), and every column whose name suggests an
age, time or interval was measured (non-null count, min, max). What each column dates is quoted
from the source's own file header, schema sheet or README wherever there is one. Where the source
says nothing, the entry says **not stated** and makes no guess. `Strat.PaleoCurrents` could not
be surveyed (Dryad refuses the automated download, HTTP 403) and is not in the registry.

### Point-event catalogues (one event, one age per row)

| Loader | Rows | Age column(s) | What it dates, per the source | Validity window |
|---|---|---|---|---|
| Kimberlites Faure2010 | 4287 | `Age1_Ma` (641 non-null, 0.8–2675) | not stated; `Age1_Methode` (Rb-Sr, U-Pb, Ar-Ar, K-Ar) and `Age1_Mineral` imply a radiometric age | none |
| Kimberlites Tappe2018 | 1129 (M21) / 1174 (Y19, M16) | `Age (recom` (0.01–2848) | "recommended" age; `Geochronol` = mineral, `Geochron_1` = method | `FROMAGE` = age, `TOAGE` = −999 |
| Carbonatites | 593 (376 dated) | `Age_ma` (0–3007) | not stated | `FROMAGE` = `Age_ma + Error_ma`, `TOAGE` = `Age_ma − Error_ma` on 99.3% of rows |
| Metamorphism | 564 | `AGE(Ga)` (0.005–3.669) | the metamorphic event; `Peak P/P-T` says which point on the P-T path | none |
| BaseMetalDeposits | e.g. VMS 947 | `Age (Ga)` | not stated (columns are just deposit, country, lon, lat, age, tonnages); taken to be mineralisation | none |
| PacificSeamountAges 2013 / 2021 | 409 / 419 | `Average_age(Ma)` / `Age` (0–81.2) | average radiometric ages (PHT2021 header: "average ages (w = 1/s^2)") | none |
| Seamounts KimWessel | 24643 | `CrustAge` (0–179) | "Age of underlying seafloor from the AGE 3.2 grid" — not the seamount's age | none |
| Seamounts SIO | 39399 (good) | none | — | none |
| MagneticPicks | 101806 | `GeeK2007` (0–166.8), `Chron` | chron-boundary age on the Gee & Kent (2007) timescale | none |
| LIPs Whittaker | 80 polygons | `FROMAGE` (5–200) | not stated as such; `AgeMethod` (Ar/Ar, K/Ar, magnetic lineations, …) | `FROMAGE` = age, `TOAGE` = −999 |
| LIPs Johansson | 2526 polygons | `FROMAGE` (0–650; 144 are 0) | not stated | `FROMAGE` = age, `TOAGE` = −999 … 0 |
| LIPs Johansson_centroids | 185 | `FROMAGE` (0–650; 2 are 0) | not stated | `FROMAGE` = age, `TOAGE` = 0 or −999 |
| LIPs UTIG | 2501 polygons | none | — | `FROMAGE` = `TOAGE` = 0 everywhere |

### Samples carrying two kinds of age

Every zircon dataset has one row per grain. Each row carries the grain's U-Pb
**crystallisation** age and the **depositional** age of the sediment sample it came from. The
Puetz 2026 schema defines its `Est./Max./Min. Stratigraphic Age` as the sample's
"stratigraphic/depositional age", `Strat/Dep Age` as "the depositional age of the sample from
which the zircons were extracted", and `U-Pb Non-Iter. Prob. age` as "the U-Pb age of the detrital
zircon". The 2018 release gives three isotope-system ages and no preferred one, so gprm sets no
`field` for it. `load_Hf` adds Hf depleted-mantle model ages (`TDM1`, `TDM2`). The Global
Tectonic Map gives each province separate magmatic and metamorphic age ranges.

### Interval-only data

pbdb (`max_ma`/`min_ma`), PaleoLithology (`FROMAGE`/`TOAGE` = the Boucot time slice, with
`ReconstructionAge` its mid-age), and SurfaceGeology (`FROMAGE`/`TOAGE` derived from era names)
have only a range, no single age.

## Problems found along the way

Fixed:
- `Seamounts('SIO_*')` skipped 17 lines of files that have no header, silently dropping the first
  17 seamounts (39382 returned against 39399 in the file and README).
- `PacificSeamountAges('2013')` had its name and two-letter code columns swapped, because it
  followed the file's commented header, which lists them in the wrong order relative to the data.
- Carbonatites' "no age" rows (`Age_ma = 0`, `Error_ma = 0`, no reference) were indistinguishable
  from 0 Ma. They are now dropped by default, and `Age` is NaN for them.
- Johansson's zero emplacement ages (144 polygons, 2 centroids) are treated as errors, handled the
  same way.
- Geochem `age`: 380,655 of its 381,573 numeric ages lie in 0–4567 Ma. The derived `Age` sets
  the 2 values above 4567 Ma (20076, 127000) and the 3 negatives larger than 10 kyr (−7.1, −5.7,
  −0.2) to NaN. It sets the 119 negatives within 10 kyr of 0 (young volcanics, apparently
  historical eruption dates entered as negative years) to 0. The 794 zeros are kept. No rows are
  dropped.

Recorded in the registry notes, not changed:
- Zircons 2018/2019 grain ages include negative values and values above 4.57 Ga (raw, unfiltered).
- Hoggard et al.'s table does not say what its age dates. It is taken to be the mineralisation
  age, a working assumption not checked against the paper.
