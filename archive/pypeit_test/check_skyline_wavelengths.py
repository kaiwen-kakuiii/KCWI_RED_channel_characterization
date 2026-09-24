#!/usr/bin/env python
"""Test a wavelength solution against the night sky, which owes it nothing.

The concern this answers: the RH1 template was seeded by holy-grail, which is
unreliable, and every slice was then solved against that one template.  A low
RMS cannot rule out a *common-mode* error -- if the template names each feature
as the neighbouring catalogue line, every slice inherits the same offset and
each fit is still tight.  The FeAr arc cannot detect that, because the arc is
what defined the scale.

Night-sky emission lines can.  They are in the science exposure, not the
calibration, they arise in the atmosphere at wavelengths fixed by atomic
physics, and nothing in the reduction ever fit to them.  If the wavelength scale
is right they land at their catalogue values; a misidentification by one arc
line would put them ~4 A out, far more than any plausible measurement error.

Catalogue sky lines are quoted in air; PypeIt works in vacuum.  They are
converted here, because a stray air value fakes a ~1.7 A offset -- the same
size as the error being looked for.

Usage:
    python check_skyline_wavelengths.py <spec2d_*.fits> [more ...]
"""
import sys

import numpy as np
from astropy import units
from astropy.io import fits
from pypeit.core.wave import airtovac

# Strong, isolated airglow lines.  Catalogues quote these in AIR, so they are
# converted once here with pypeit's own airtovac -- mixing an air value into a
# vacuum comparison fakes a ~1.7 A offset at these wavelengths, which is the
# same size as the error this test exists to look for.
# [OI] are forbidden atomic lines -- true singlets, and the only ones whose
# catalogue wavelength stays meaningful at any resolution.  The OH features are
# molecular doublets: at the Large slicer (1.42 A) they merge and a single
# catalogue value describes them, but at the Small slicer (0.69 A) they resolve
# and the centroid lands on one component, scattering the answer by +-1.5 A.
# Measured on 2023-11-08: OH gave 1.01 A scatter where [OI] alone gave 0.02 A.
# So the verdict below rests on singlets whenever any is in range.
SINGLETS = {"[OI] 6300": 6300.304, "[OI] 6364": 6363.776}
BLENDS = {"OH 6499": 6498.729, "OH 6554": 6553.617,
          "OH 6828": 6828.478, "OH 6922": 6922.089}
_AIR = {**SINGLETS, **BLENDS}
SKY_LINES = {k: float(airtovac(np.array([v]) * units.AA)[0].value)
             for k, v in _AIR.items()}


def measure(path, search=3.0):
    hdu = fits.open(path)
    det = [e.name for e in hdu if e.name.endswith("SKYMODEL")][0].split("-")[0]
    sky = hdu[f"{det}-SKYMODEL"].data
    wave = hdu[f"{det}-WAVEIMG"].data
    good = np.isfinite(wave) & (wave > 100)

    print(f"\n=== {path.split('/')[-1]}")
    print(f"  wavelength image spans {wave[good].min():.1f}-{wave[good].max():.1f} A")

    # Collapse to a 1D sky spectrum on a regular grid: every pixel carries its
    # own wavelength, so the sky is binned by wavelength rather than by column.
    lo, hi = wave[good].min(), wave[good].max()
    step = 0.3
    edges = np.arange(lo, hi, step)
    idx = np.digitize(wave[good], edges) - 1
    vals = sky[good]
    ok = (idx >= 0) & (idx < len(edges) - 1)
    prof = np.zeros(len(edges) - 1)
    cnt = np.zeros(len(edges) - 1)
    np.add.at(prof, idx[ok], vals[ok])
    np.add.at(cnt, idx[ok], 1)
    prof = np.divide(prof, np.maximum(cnt, 1))
    centres = 0.5 * (edges[:-1] + edges[1:])

    offsets = []
    for name, lam in SKY_LINES.items():
        if not (lo + search < lam < hi - search):
            continue
        sel = np.abs(centres - lam) < search
        if sel.sum() < 5 or cnt[sel].min() == 0:
            continue
        w, f = centres[sel], prof[sel]
        base = np.percentile(f, 20)
        f = f - base
        if f.max() <= 0:
            continue
        # flux-weighted centroid over the pixels above half maximum
        m = f > 0.5 * f.max()
        if m.sum() < 2:
            continue
        meas = float(np.sum(w[m] * f[m]) / np.sum(f[m]))
        offsets.append((name, meas - lam))
        kind = "singlet" if name in SINGLETS else "blend  "
        print(f"  {name:12s} {kind} expected {lam:8.2f}  measured {meas:8.2f}  "
              f"offset {meas - lam:+6.2f} A")

    if not offsets:
        print("  no catalogued sky line inside this wavelength range")
        return None
    sing = np.array([d for n, d in offsets if n in SINGLETS])
    allo = np.array([d for _, d in offsets])
    if sing.size:
        use, label = sing, f"{sing.size} singlet"
    else:
        use, label = allo, f"{allo.size} blended (no singlet in range)"
    print(f"  --> all {allo.size} lines: mean {allo.mean():+.2f} A, "
          f"scatter {allo.std():.2f} A")
    print(f"  --> verdict from {label}: mean {use.mean():+.2f} A, "
          f"scatter {use.std():.2f} A")
    verdict = ("scale is correct; a one-line misidentification would show ~4 A"
               if np.abs(use.mean()) < 1.5 else
               "OFFSET IS LARGE -- suspect a misidentified line")
    if not sing.size:
        verdict += " (blends only -- treat as indicative, not conclusive)"
    print(f"  --> {verdict}")
    return use.mean()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for p in sys.argv[1:]:
        measure(p)
