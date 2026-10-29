# Bug: local beam passed downstream after the mask

In all four scripts, the element after the 10 m mask receives the beam in the mask's
(or screen's) **local** coordinates instead of the global beam. The crystals then
effectively sit 10 m farther from the source than intended. The Fcurv optimization hides
this by weakening the sagittal curvature, which is why Fcurv came out above 1 in every
run.

The bug is **fixed** in the five scripts (`XRT_SSU_BLM.py`, `XRT_SSU_BLM_PP.py`,
`XRT_APS_BLM.py`, `XRT_SSU_SBDL.py`, `XRT_BNL_SBDL.py`). The results in `run_XRT_*` and
`RUN_XRT.md` were computed before the fix. The scans have been rerun with the fix
(`run_XRT_*_fixed/`); see `RERUN_RESULTS.md` and `RUN_XRT.md`.

## The bug

`XRT_SSU_BLM.py`, `run_process`:

```python
beamSource = (beamLine.source.shine())
beamMask1 = (beamLine.mask1.propagate(beamSource))
beamSCR1 = (beamLine.screen1.expose(beamMask1))
beamBraggGlobal, beamBraggLocal = (beamLine.Bragg.reflect(beam=rs.Beam(copyFrom=beamMask1)))
```

In xrt (`xrt/backends/raycing/apertures.py`), `RectangularAperture.propagate(beam)`
works as follows:

- It marks the blocked rays as lost in the **global** beam it receives (`beam.state`).
- It **returns the local beam** in the aperture frame: "Returned values: beamLocal".

`Screen.expose()` also returns a local beam.

**What the local beam looks like:** its rays sit at y = 0 (checked: the local y range at
the mask is exactly 0), with the transverse positions they had at the mask plane.

**What the next element makes of it:** it treats this beam as global, i.e. as the real
beam shifted 10 m back along the beam axis.

- Every downstream element is effectively 10 m (`d_slit1`) farther from the source.
- For SSU_BLM: Bragg at 41 m instead of 31 m, Laue at 41.53 m instead of 31.53 m.
- The screen distance from the Laue crystal is unchanged.

**Why the Fcurv optimization hid it:**

- The Laue sagittal radius RsL is computed for p = 31.53 m. With the effective
  p = 41.53 m, the focus moves.
- The optimizer then finds the Fcurv that puts the focus back on the screen, 1.1175 at
  80 keV.
- With the global beam, the analytic RsL already focuses on the screen with Fcurv = 1.0.

## Where it occurs

| Script | Line | Element receiving the local beam |
|---|---|---|
| `XRT_SSU_BLM.py` | 498 | Bragg crystal gets `beamMask1` (`mask1.propagate`) |
| `XRT_APS_BLM.py` | 511 | Bragg crystal gets `beamMask1` (`mask1.propagate`) |
| `XRT_SSU_SBDL.py` | 332 | Filter1 gets `beamMask1` (`Mask1.propagate`) |
| `XRT_BNL_SBDL.py` | 559 | Filter1 gets `beamSCR1` (`Screen1.expose`) |

## Test: SSU_BLM at 80 keV

- **Setup:** the script's own beamline (`build_beamline` plus `energy_scan` at 80 keV,
  200k rays), at the script's 80 keV values: RsL = 1217.64 mm, RmL = 76039.41 mm,
  RmB = 647668 mm.
- **Two cases:**
  - Script: the Bragg crystal gets the local mask beam.
  - Correct: the mask is applied, and the Bragg crystal gets the global `beamSource`.
- **Measurement:** FWHM from 200-bin histograms over x = ±0.2 mm and z = ±10 mm, as in
  the script's plots.

| Beam into Bragg | Fcurv | x FWHM at 56 m | Horizontal focus | z FWHM at 56 m | Σ intensity |
|---|---|---|---|---|---|
| script (local) | 1.1175 | 28 µm | 56.00 m | 4.0 mm | 4857 |
| script (local) | 1.0 | 374 µm | 52.25 m | 5.1 mm | 4874 |
| **correct (global)** | **1.0** | **40 µm** | **56.00 m** | **2.9 mm** | 4776 |
| correct (global) | 1.1175 | 380 µm | 61.75 m | 2.6 mm | 4777 |

- **The script's case reproduces its reported values.** The full run gave a 34 µm and
  3.9 mm FWHM, with the horizontal focus at 56.008 m.
- **With the correct geometry the analytic radii work.**
  - Fcurv = 1.0 focuses exactly on the screen.
  - The horizontal spot is about 40 µm.
  - The vertical size at the screen is smaller: 2.9 mm instead of 4.0 mm.
- **Flux changes by about 2%,** because the ray directions are the same and only the
  positions shift.
- **OASYS gives the correct result.** The OASYS2-XRT Runner keeps the global beam after
  an aperture, which is why the OASYS workspace
  (`oasys_SSU_BLM/SSU_BLM_80keV.ows`) gives the correct result with Fcurv = 1.0. Run
  through the Runner with 100k rays, it gives flux 3.38×10¹², horizontal FWHM 52 µm
  (1024-bin plot) and vertical FWHM 3.4 mm at the screen.

## Consequences for the runs

- **Affected results:** focus sizes, focus positions and the Fcurv values reported in
  `RUN_XRT.md`, for all four beamlines. They come from the distorted geometry.
  - The note "Fcurv: in every case the optimum is above 1" is an artefact of this bug.
  - The vertical focusing of the Bragg crystal (RmB, computed for 31 m) and the M1
    mirror of BNL_SBDL are also evaluated at the wrong source distance.
- **Flux changes depend on clipping downstream.** Ray directions are unaffected, but
  the beam is wider at every element.
  - SSU_BLM at 80 keV: about 2%.
  - APS_BLM at 44 keV: about 2%. With the global beam, Fcurv = 1.0 focuses horizontally
    at 52.1 m, on the 52 m screen; the script used 1.1505. See
    `oasys_APS_BLM/INFO.md`.
  - SSU_SBDL at 50 keV: about +30% intensity with the global beam. Fcurv ≈ 0.985 focuses
    on the 56 m screen; the script used 1.106. See `oasys_SSU_SBDL/INFO.md`.
  - **BNL_SBDL at 40 keV: the correct geometry gives about 2× the flux** (1.04×10¹³ vs
    5.24×10¹² ph/s). In the script's geometry the horizontal beam is wider at the 50 mm
    M1 mirror, so more is clipped. See `oasys_BNL_SBDL/INFO.md`.
  - BNL_SBDL also needs a much smaller correction in the correct geometry: Fcurv ≈ 1.02
    instead of 1.199 at 40 keV.
  - Bandwidths should change little; this was not checked separately.
- **This bug is independent of the ML reflectivity bug** (`BUG_ML.md`).

## Fix

Pass the global beam downstream. The aperture has already marked the blocked rays in it:

```python
beamSource = beamLine.source.shine()
beamMask1 = beamLine.mask1.propagate(beamSource)    # local beam, for plots only
beamSCR1 = beamLine.screen1.expose(beamSource)       # local beam, for plots only
beamBraggGlobal, beamBraggLocal = beamLine.Bragg.reflect(beam=rs.Beam(copyFrom=beamSource))
```

- Use the same pattern in `XRT_APS_BLM.py`.
- In `XRT_SSU_SBDL.py`, pass `beamSource` to `Filter1.double_refract`.
- In `XRT_BNL_SBDL.py`, pass `beamSource` to `Filter1.double_refract`.
- Alternatively, use `propagate(beamSource, needNewGlobal=True)` and use its global
  output.

## Check of the fix

The scripts' own `run_process` was run at one energy each (100k rays), including the broad
Fcurv search:

| Script | Fcurv found (before the fix) | Horizontal focus (script caustic) | Screen |
|---|---|---|---|
| SSU_BLM, 80 keV | 1.0001 (1.1175) | 54.80 m | 56.0 m |
| SSU_SBDL, 50 keV | 0.986 (1.106) | 55.95 m | 56.0 m |

- The optimizer now finds Fcurv ≈ 1, as expected with the analytic radii.
- For SSU_BLM the script's FWHM caustic puts the narrowest point at 54.8 m, with 15 µm
  FWHM. A direct trace with Fcurv = 1.0 gave the minimum rms width at 56.0 m. This
  suggests a narrow core focusing slightly upstream of broader wings; it was not
  investigated further.
- APS_BLM and BNL_SBDL were not rerun after the fix.

## To do

Done: full energy scans rerun and `RUN_XRT.md` updated. The remaining open points are in
the TODO list of `RUN_XRT.md`.
