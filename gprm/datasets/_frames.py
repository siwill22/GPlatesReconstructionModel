"""What the anchor plate ids mean in each reconstruction model gprm loads.

Anchoring a reconstruction on plate 0 puts it in *some* absolute reference frame, but which one
depends on the model, and a few models carry more than one: Torsvik & Cocks (2017) put a mantle
frame on plate 0 and a palaeomagnetic one on plate 1, Muller et al (2025) keep Cao et al's
palaeomagnetic frame on plate 5. Nothing in the rotation model says so short of reading its
comments. This module records, per model, the plate ids a user would deliberately anchor on and
what each is fixed to, so that the meaning is discoverable from the loaded model. The evidence
for each entry is in ``docs/reference_frames.md``.

A frame is classified by what it is fixed to, not by how it was built. There are two:

``'mantle'``
    fixed to the deep mantle. Includes frames built from moving or fixed hotspots, from
    palaeomagnetic data corrected for true polar wander, from optimisation (net rotation,
    trench migration, hotspot fit), no-net-rotation frames, and Pacific hotspot frames. A frame
    that switches method with age ("hybrid") is still one mantle frame.
``'spin axis'``
    palaeomagnetic, fixed to the Earth's spin axis, without a true polar wander correction.

Intermediate plates in a rotation chain (EarthByte's 004 true polar wander correction and 070
longitude shift, zero-rotation placeholders such as 001) are deliberately not listed: anchoring
on them gives a frame that is not fixed to anything tangible.

Each fetch_ function sets ``model.reference_frames`` from here. An entry has these keys:

``plate_id``
    the anchor plate id.
``reference``
    ``'mantle'`` or ``'spin axis'`` -- one of :data:`REFERENCES`.
``description``
    which frame it is, in a few words.
``valid``
    ``(young, old)`` ages in Ma over which the model defines the frame, or ``None`` if it spans
    the whole model. Anchoring outside it does not fail, but the anchor then has no rotation,
    so the plate it hangs from (the Pacific, for a Pacific hotspot frame) is simply held at its
    present-day position: the result is not in the frame described.
``source``
    where the classification comes from: a quotation from the rotation file, or the gprm
    maintainer where the files do not say.
``note``
    anything a user should know before relying on the frame, or ``None``.

MIT License

Copyright (c) 2017-2021 Simon Williams
"""

import copy as _copy

#: What a reference frame can be fixed to.
REFERENCES = ('mantle', 'spin axis')

_MAINTAINER = 'gprm maintainer (not stated in the model files)'


def _frame(plate_id, reference, description, source, valid=None, note=None):
    if reference not in REFERENCES:
        raise ValueError('reference must be one of {}, not {!r}'.format(REFERENCES, reference))
    return dict(plate_id=plate_id, reference=reference, description=description,
                valid=valid, source=source, note=note)


# Pacific hotspot frame (plate 002 relative to Pacific 901) shared by the EarthByte-lineage models
def _pacific_wk08a(young_old):
    return _frame(2, 'mantle', 'Pacific hotspots, Wessel & Kroenke (2008) WK08-A',
                  source='rotation file: "PHS-PAC Pacific hotspot-Pacific ... Model WK08-A"',
                  valid=young_old)


_EARTHBYTE_MANTLE_SOURCE = ('rotation file: Africa-to-frame poles "Global Moving Hotspot" '
                            '(Torsvik et al 2008) at young ages, then "TPW correction FROM PM"')

_MULLER2019_MANTLE = _frame(
    0, 'mantle', 'optimised mantle reference frame (Muller et al 2019)',
    source='rotation file: "AFR-000 ... Optimised Absolute Reference Frame, model v2.0"')

_MULLER2022_MANTLE_SOURCE = ('file name 1000_0_rotfile_MantleOpt.rot; the frame plates are merged away '
                             '("Removed fixed plate 5", "Removed fixed plate 70")')

_CAO_SPIN_AXIS = _frame(
    0, 'spin axis', 'palaeomagnetic reference frame',
    source='rotation file: "AFR-070 Africa-longitude correction (GAPWAP)", then '
           'cratons to "000 ... Spin Axis", with 070-000 a zero rotation')

_MAINTAINER_SPIN_AXIS = _frame(0, 'spin axis', 'palaeomagnetic reference frame', source=_MAINTAINER)


def _cao_toy(description, source):
    return [_frame(0, 'mantle', description, source=source), _pacific_wk08a((0., 144.))]


_REGISTRY = {
    'TorsvikCocks2017': [
        _frame(0, 'mantle', 'global moving hotspot frame to 120 Ma, then palaeomagnetic frame '
                            'corrected for true polar wander',
               source='rotation file: "Super duper HYBRID (PID=0) is the GMHRF for ages <= 120 '
                      'Ma, and PM-HYBRID corrected for TPW for 120-550 Ma"'),
        _frame(1, 'spin axis', 'PM-HYBRID palaeomagnetic frame',
               source='rotation file: "PM-HYBRID (PID=1) is a PM frame with longitudes adjusted '
                      'to match 1) TPW observed in GMHFR for the last 120 Ma, 2) LIPs and '
                      'kimberlites for 250-550 Ma, and 3) a smooth transition imposed between '
                      '120 and 250 Ma"'),
        _frame(2, 'mantle', 'Pacific hotspots, Duncan & Clague (1985)',
               source='rotation file: "HOTSPOTS-PAC DUNCAN & CLAGUE (1985)"', valid=(0., 149.5)),
        _frame(5, 'mantle', 'Pacific hotspots, Koppers et al (2001)',
               source='rotation file: "PAC HS-PAC Koppers01"', valid=(0., 140.)),
        _frame(6, 'mantle', 'Pacific hotspots, Wessel & Kroenke (1997)',
               source='rotation file: "PAC HS-PAC Wessel and Kroenke 97"', valid=(0., 145.)),
    ],
    'DomeierTorsvik2014': [
        _frame(0, 'mantle', 'palaeomagnetic frame corrected for true polar wander',
               source='rotation file: 001-000 "TPW corrections FROM PM"'),
        _frame(1, 'spin axis', 'palaeomagnetic frame (PM60)',
               source='rotation file: plates attached to 001 as "PM60"; 001-000 "no TPW" / '
                      '"TPW corrections FROM PM"'),
    ],
    'Seton2012': [
        _frame(0, 'mantle', 'Indo-Atlantic moving hotspots (O\'Neill et al 2005), then '
                            'palaeomagnetic frame corrected for true polar wander',
               source='rotation file: "AFR-AHS O\'Neill et.al 2005", then "AFR-AHS PM corrected '
                      'TPW derived from Steinberger & Torsvik 2008"'),
        _frame(2, 'mantle', 'Pacific hotspots, Wessel & Kroenke (2008) WK08-A',
               source='rotation file: "PHS-PAC WK08-A Wessel & Kroenke 2008"', valid=(0., 200.)),
    ],
    'Muller2016': [
        _frame(0, 'mantle', 'global moving hotspots (Torsvik et al 2008), then palaeomagnetic '
                            'frame corrected for true polar wander',
               source='rotation file: "AFR-GLH ... Global Moving Hotspot", then "TPW-corrected '
                      'PMAG modified to include a 10 deg. longitudinal shift"'),
        _pacific_wk08a((0., 144.)),
    ],
    'Matthews2016': [
        _frame(0, 'mantle', 'global moving hotspots (Torsvik et al 2008), then palaeomagnetic '
                            'frame corrected for true polar wander',
               source=_EARTHBYTE_MANTLE_SOURCE),
        _pacific_wk08a((0., 144.)),
    ],
    'Young2019': [
        _frame(0, 'mantle', 'global moving hotspots (Torsvik et al 2008), then palaeomagnetic '
                            'frame corrected for true polar wander',
               source=_EARTHBYTE_MANTLE_SOURCE,
               note='The 410-250 Ma rotation file applies no true polar wander correction '
                    '("No true polar wander correction"), so before 250 Ma the frame rests on '
                    'palaeomagnetic poles alone. Classified as mantle by the gprm maintainer.'),
        _pacific_wk08a((0., 144.)),
    ],
    'CaoToyRodinia:NNR': _cao_toy('no-net-rotation frame', 'rotation file: 005-001 "Removing '
                                  'net rotation - calculated outside GPlates"'),
    'CaoToyRodinia:OV': _cao_toy('orthoversion frame', _MAINTAINER + '; pre-250 Ma frame node '
                                 'labelled "orthoversion"'),
    'CaoToyRodinia:SSL': _cao_toy('slow-continent frame', _MAINTAINER + '; pre-250 Ma frame node '
                                  'labelled "slow continent"'),
    'Muller2019': [_MULLER2019_MANTLE],
    'Clennett:M2019': [_MULLER2019_MANTLE],
    'Clennett:S2013': [
        _frame(0, 'mantle', 'Indo-Atlantic hotspots: O\'Neill et al (2005) moving hotspots, '
                            'then Muller et al (1993)',
               source='rotation file: "AFR-AHS O\'Neill et.al 2005", then "AFR-AHS Stage rot '
                      'from Muller et.al. 1993"'),
        _frame(2, 'mantle', 'Pacific hotspots, Wessel & Kroenke (2008) WK08-A',
               source='rotation file: "PHS-PAC WK08-A Wessel & Kroenke 2008"', valid=(0., 144.)),
    ],
    'Muller2022:Opt': [
        _frame(0, 'mantle', 'optimised mantle reference frame (Muller et al 2022)',
               source=_MULLER2022_MANTLE_SOURCE)],
    'Muller2022:NNR': [
        _frame(0, 'mantle', 'no-net-rotation frame (Muller et al 2022)',
               source='file name no_net_rotation_model.rot',
               note='The pole comments in this file read "Optimised mantle reference frame"; '
                    'the file name is what identifies it as no-net-rotation.')],
    'Muller2025:Opt': [
        _frame(0, 'mantle', 'optimised mantle reference frame (Muller et al 2025)',
               source='rotation file: 005-000 "optAPM"; archive README'),
        _frame(5, 'spin axis', 'palaeomagnetic reference frame of Cao et al (2024)',
               source='archive README: "setting 005 as anchor plate ID removes the optimised '
                      'rotations"'),
    ],
    'Muller2025:NNR': [
        _frame(0, 'mantle', 'no-net-rotation frame (Muller et al 2025)',
               source='rotation file no_net_rotation_model_20240725.rot, 005-000'),
        _frame(5, 'spin axis', 'palaeomagnetic reference frame of Cao et al (2024)',
               source='archive README: "setting 005 as anchor plate ID removes the optimised '
                      'rotations"'),
    ],
    'Merdith2021': [_CAO_SPIN_AXIS],
    'Cao2024': [_CAO_SPIN_AXIS],
    'Li2008': [_MAINTAINER_SPIN_AXIS],
    'Li2023:East': [_MAINTAINER_SPIN_AXIS],
    'Li2023:West': [_MAINTAINER_SPIN_AXIS],
    'Pehrsson2015': [_MAINTAINER_SPIN_AXIS],
    'Golonka2007': [
        _frame(0, 'spin axis', 'palaeomagnetic reference frame',
               source='rotation file: 701-000 "Africa - Spin Axis", "Torsvik and van der Voo '
                      '(2002), spherical spline option GAD (dipole)"')],
    'Scotese2008': [
        _frame(0, 'spin axis', 'PALEOMAP palaeomagnetic reference frame', source=_MAINTAINER,
               note='Loosely defined: the rotation file carries an all-zero 001-000 pole '
                    'labelled "Hot Spot to PMAG", and the model is sometimes used as a mantle '
                    'frame. Use with caution in either role.')],
}


def models():
    """Names of every reconstruction model with recorded reference frames."""
    return sorted(_REGISTRY)


def reference_frames(model):
    """Return the anchor plates of a reconstruction model and what each is fixed to.

    :param model: registry name, e.g. ``'TorsvikCocks2017'`` or ``'Muller2025:NNR'``.
        :func:`models` lists them.
    :returns: list of dicts (see the module docstring for the keys), ordered by plate id.
    """
    if model not in _REGISTRY:
        raise KeyError('No reference frames recorded for {!r}. Known models: {}'.format(
            model, ', '.join(models())))
    return _copy.deepcopy(_REGISTRY[model])
