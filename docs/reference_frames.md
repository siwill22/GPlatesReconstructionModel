# What the anchor plate ids mean in each gprm reconstruction model

Anchoring a reconstruction on plate 0 puts it in *an* absolute reference frame, but which one
depends on the model, and some models carry more than one. This page gives the rules gprm uses
to classify them and the evidence behind each model's entry. The registry itself is in
[gprm/datasets/_frames.py](../gprm/datasets/_frames.py). Terms follow
[CONTEXT.md](../CONTEXT.md) (*Reference Frame*, *Anchor Plate*).

## Using it

Every `fetch_*` model carries its frames, and `info()` prints them:

```python
from gprm.datasets import Reconstructions
model = Reconstructions.fetch_TorsvikCocks2017()
model.info()
# ...
# Reference Frames (anchor_plate_id):
#    - 0: mantle -- global moving hotspot frame to 120 Ma, then palaeomagnetic frame corrected for true polar wander
#    - 1: spin axis -- PM-HYBRID palaeomagnetic frame
#    - 2: mantle -- Pacific hotspots, Duncan & Clague (1985), 0-149.5 Ma only
#    ...
model.reference_frames        # the same, as a list of dicts with source and note
```

An entry can be read without downloading the model:

```python
from gprm.datasets import reference_frames
reference_frames('Muller2025:Opt')
```

For a model you assemble yourself, declare its frames with
`ReconstructionModel.add_reference_frame(plate_id, reference, description, valid=None,
source=None, note=None)`. A model with none prints "not documented".

## Rules

1. **A frame is classified by what it is fixed to, not by how it was built.** There are two
   references:
   - **mantle**: fixed to the deep mantle. Frames built from moving or fixed hotspots, from
     palaeomagnetic data corrected for true polar wander (TPW), by optimisation, as
     no-net-rotation, and Pacific hotspot frames are all mantle frames. A frame that changes
     method with age (moving hotspots to ~100 Ma, TPW-corrected palaeomagnetism before) is
     still one mantle frame, not a mantle frame at one time and a palaeomagnetic one at another.
   - **spin axis**: palaeomagnetic, without a TPW correction.
2. **Only frames a user would deliberately anchor on are listed**: each must relate to a
   tangible mantle or core concept. Intermediate plates in a rotation chain are left out, as
   are zero-rotation placeholders that merely duplicate plate 0: in the EarthByte-lineage
   models, 004 (TPW correction), 070 (longitude shift) and 001; in Torsvik & Cocks (2017), 11
   (PM-ZeroAfr), 3 (Atlantic hotspots) and 4 (Kerguelen–Réunion hotspots relative to India).
3. **A frame's valid span is recorded where it is shorter than the model's.** Outside it,
   anchoring still works but the anchor has no rotation, so the plate it hangs from (the
   Pacific, for a Pacific hotspot frame) is held at its present-day position. Every span below
   ends at the last pole the rotation file gives that plate.
4. **Evidence is a quotation from the rotation file** wherever the file says; where it does
   not, the classification is the gprm maintainer's and is marked so.

## Evidence per model

Quotations are from the pole comments of the rotation files each fetcher loads.

| Model (registry name) | Plate | Reference | Evidence |
|---|---|---|---|
| `TorsvikCocks2017` | 0 | mantle | "Super duper HYBRID (PID=0) is the GMHRF for ages <= 120 Ma, and PM-HYBRID corrected for TPW for 120-550 Ma" |
| | 1 | spin axis | "PM-HYBRID (PID=1) is a PM frame with longitudes adjusted to match 1) TPW observed in GMHFR for the last 120 Ma, 2) LIPs and kimberlites for 250-550 Ma, and 3) a smooth transition imposed between 120 and 250 Ma"; all plates hang from 001 |
| | 2 | mantle, 0–149.5 Ma | 002→901 "HOTSPOTS-PAC DUNCAN & CLAGUE (1985)"; poles repeat unchanged from 149.5 to 600 Ma |
| | 5 | mantle, 0–140 Ma | 005→901 "PAC HS-PAC Koppers01" |
| | 6 | mantle, 0–145 Ma | 006→901 "PAC HS-PAC Wessel and Kroenke 97" |
| `DomeierTorsvik2014` | 0 | mantle | 001→000 "TPW corrections FROM PM" |
| | 1 | spin axis | plates attach to 001 as "PM60"; 001→000 is the TPW correction |
| `Seton2012` | 0 | mantle | 701→001 "AFR-AHS O'Neill et.al 2005", then "AFR-AHS PM corrected TPW derived from Steinberger & Torsvik 2008"; 001→000 zero, "AHS-HOT Present day Atlantic-Indian hotspots fixed to 000" |
| | 2 | mantle, 0–200 Ma | 002→901 "PHS-PAC WK08-A Wessel & Kroenke 2008" |
| `Clennett:S2013` | 0 | mantle | 701→001 "AFR-AHS O'Neill et.al 2005", then "AFR-AHS Stage rot from Muller et.al. 1993" |
| | 2 | mantle, 0–144 Ma | 002→901 "PHS-PAC WK08-A Wessel & Kroenke 2008" |
| `Muller2016` | 0 | mantle | 701→001 "AFR-GLH … Global Moving Hotspot", then "TPW-corrected PMAG modified to include a 10 deg. longitudinal shift"; 001→000 zero |
| | 2 | mantle, 0–144 Ma | 002→901 "PHS-PAC Pacific hotspot-Pacific … Model WK08-A" |
| `Matthews2016` | 0 | mantle | 701→070 "AFR-LGS … Global Moving Hotspot", then "PMAG reference frame"; 070→004 longitude shift; 004→001 "TPW correction FROM PM" (non-zero through 410 Ma); 001→000 zero |
| | 2 | mantle, 0–144 Ma | as Muller2016 |
| `Young2019` | 0 | mantle | as Matthews2016 to 250 Ma. **Note:** the 410–250 Ma file sets 004→001 to zero ("No true polar wander correction") and attaches Africa via "Torsvik and Van der Voo … GAD" poles, so before 250 Ma the frame rests on palaeomagnetism alone. Classified mantle by the maintainer. |
| | 2 | mantle, 0–144 Ma | as Muller2016 |
| `CaoToyRodinia:NNR` | 0 | mantle | 005→001 "Removing net rotation - calculated outside GPlates" |
| `CaoToyRodinia:OV` | 0 | mantle | pre-250 Ma frame node labelled "orthoversion"; maintainer |
| `CaoToyRodinia:SSL` | 0 | mantle | pre-250 Ma frame node labelled "slow continent"; maintainer |
| (all three) | 2 | mantle, 0–144 Ma | as Muller2016 |
| `Muller2019`, `Clennett:M2019` | 0 | mantle | 701→000 "AFR-000 … Optimised Absolute Reference Frame, model v2.0" |
| `Muller2022:Opt` | 0 | mantle | file name `1000_0_rotfile_MantleOpt.rot`; the frame plates are merged away ("Removed fixed plate 5", "Removed fixed plate 70") |
| `Muller2022:NNR` | 0 | mantle | identified by the file name `no_net_rotation_model.rot`. **Note:** its pole comments read "Optimised mantle reference frame" |
| `Muller2025:Opt` | 0 | mantle | 005→000 "optAPM"; archive README |
| | 5 | spin axis | archive README: "setting 005 as anchor plate ID removes the optimised rotations", i.e. the palaeomagnetic frame of Cao et al (2024) |
| `Muller2025:NNR` | 0 | mantle | `no_net_rotation_model_20240725.rot`, 005→000 |
| | 5 | spin axis | as `Muller2025:Opt` |
| `Merdith2021`, `Cao2024` | 0 | spin axis | 701→070 "AFR-070 Africa-longitude correction (GAPWAP)" with 070→000 zero; older cratons →000 "… Spin Axis" |
| `Li2008` | 0 | spin axis | not stated in the file; maintainer |
| `Li2023:East`, `Li2023:West` | 0 | spin axis | the file has no comments; maintainer |
| `Pehrsson2015` | 0 | spin axis | not stated in the file (plate ids 1–99 there are real plates, not frames); maintainer |
| `Scotese2008` | 0 | spin axis | maintainer. **Note:** loosely defined — an all-zero 001→000 pole is labelled "Hot Spot to PMAG", and the model is sometimes used as a mantle frame. Use with caution in either role. |

`fetch_Golonka` has no entry: its file path is currently broken, so its rotation file was never
examined.
