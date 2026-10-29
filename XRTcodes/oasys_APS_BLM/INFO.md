# APS_BLM in OASYS2-XRT (44 keV)

This folder holds an OASYS2 workspace for the APS_BLM beamline of `XRT_APS_BLM.py` at one
energy, 44 keV, in the middle of the script's flat flux region (38–54 keV). The workspace
uses OASYS2-XRT widgets only. It has the same structure as `../oasys_SSU_BLM`.

> **Status (2026-10-29): with the fixed ML reflectivity this workspace, like the script,
> transmits almost nothing.** The Bragg and Laue reflection windows no longer overlap: the
> Laue crystal sits at its nominal angle, and the model uses the crystal mid-plane as the
> depth reference. See `../RUN_XRT.md` ("Bragg+Laue beamlines") and its TODO list. The
> workspace needs the same changes as the script (entrance-surface depth reference, Laue
> pitch offset; +15 µrad at 44 keV in the test run).

## Files

| File | Content |
|---|---|
| `APS_BLM_44keV.ows` | OASYS2 workspace |
| `APS_BLM_44keV_runner_script.py` | Python script that the XRT Runner widget generates from this workspace (500k rays) |

## Workspace

Six widgets in a chain, each sending XRTData to the next:

| # | Widget | Name | Settings |
|---|---|---|---|
| 0 | Bending Magnet Light Source | `DW` | B0 = 0.6 T, 7 GeV, 0.1 A; σx = 91 µm, σz = 30 µm; εx = 6.37 nm·rad, εz = 1.29 nm·rad; acceptance 1.0 × 0.5 mrad; E = 43.78–44.22 keV (±0.5%); 500 000 rays |
| 1 | Rectangular Aperture | `Mask1` | at [0, 10000, 0]; gap 10 × 5 mm |
| 2 | Johann Cylinder Crystal | `Bragg` | Si(311) Bragg, symmetric, t = 10 mm; at [0, 26500, 0]; pitch `np.pi - alpha - thetaB`; Rm = 217527.5769 mm; ML; tensor oriented with nrm [5,−1,−1], sag [0,1,−1]; limPhysY ±500 mm |
| 3 | Bent Crystal 2D | `Laue` | Si(311) Laue, nrm [5,−1,−1], sag [0,1,−1], t = 0.6 mm; alpha `-alpha_hkl_nrm` (−48.97°); pitch `np.pi/2 + alpha - thetaB`; at [0, 26787.303927, −50]; Rs = 1686.0604 mm, Rm = −28546.8644 mm, **Fcurv = 1.0**; ML, `ML_sign` −1 |
| 4 | Screen | `SCRexp` | at [0, 52000, −50]; plotted with limits x ±1000 µm, z ±5000 µm |
| 5 | XRT Runner Python Script | | beamline name `APS_BLM_44keV`, 1 repetition, scores written to `APS_BLM_44keV/scores/` |

- **Where the values come from:** the script's own `energy_scan` at 44 keV, with θB from
  xrt (refraction-corrected) and Bragg–Laue distance 287.30 mm. The formulas are:
  - RsL = 2 sin θB sin α · fs, with p = 26.79 m and q = 25.21 m
  - RmL = RsL / (c_fac · poL), with c_fac = 0.7
  - RmB from the script's vertical-focusing formula
- **Running it:** the crystal, aperture and screen widgets have automatic execution
  switched on. Clicking "Send Data" on the Bending Magnet widget should run the whole
  chain.
- **Changing the energy:** update the following, using the formulas in
  `XRT_APS_BLM.py` (`energy_scan`):
  - `eMin` and `eMax`, and both `alignment_energy` fields
  - RmB, Rs and Rm of the Laue crystal
  - the Laue crystal y position, 26500 + 50/tan(2θB)

## Results (44 keV)

The generated Runner script was run outside the GUI with 200 000 rays. The 500k-ray
script in this folder has not been run.

| Quantity | Workspace, current code | Script after the fixes | Script test run (entrance ref. + Laue rocking) | Workspace before the ML fix | Script before the fixes (Fcurv = 1.1505) |
|---|---|---|---|---|---|
| Flux at screen (ph/s) | **1.20×10⁸** | 1.20×10⁸ | 1.90×10¹¹ | 1.09×10¹² | 1.05×10¹² |
| Horizontal FWHM at screen | 0.209 mm | 0.203 mm | 0.206 mm | 0.198 mm | 0.148 mm |
| Vertical FWHM at screen | 1.06 mm | 0.85 mm | 0.70 mm | 0.90 mm | 1.95 mm |

- **Current code:** the workspace agrees with the script after the fixes: about 10⁻⁴ of
  the old flux, so the sizes come from few rays.
- **The test run** shows what this beamline gives once the crystals are aligned. It needs
  changes that are not yet in the code or the widgets.
- The notes below compare the workspace and the script before the ML fix.

- **Measurement:**
  - The workspace values come from the Runner's plot: 1024 bins, `plot.dx` and `plot.dy`.
  - The script's values come from `fwhm_from_hist` on 200-bin plots over x = ±0.2 mm.
  - At 0.15–0.20 mm the horizontal profile fills most of the script's ±0.2 mm window, so
    that window probably narrows the script's value.
- **Flux** agrees within 3%.
- **Vertical size:** the vertical beam is about half as large as in the script (see
  difference 1 below).
- **Reflectivity:** these pre-fix values used the old internal ML reflectivity: Laue too
  high, Bragg without absorption (see `../BUG_ML.md`).

## Differences from `XRT_APS_BLM.py`

1. **Geometry: global beam after the mask (the main difference).**
   - The script passes the mask's local beam to the Bragg crystal. That effectively moves
     both crystals 10 m farther from the source (`../BUG_GEOMETRY.md`).
   - The OASYS Runner keeps the global beam, which is correct.
   - The table compares the two at 44 keV, with the same 200k source rays and FWHM from
     200-bin histograms over x = ±1 mm:

     | | Fcurv | Horizontal focus | x FWHM at 52 m | z FWHM at 52 m | Vertical focus |
     |---|---|---|---|---|---|
     | script (local beam) | 1.1505 | 51.9 m | 0.16 mm | 1.8–2.0 mm | beyond the screen (broad) |
     | correct (global beam) | 0.98 | 51.0 m | 1.01 mm | 1.1 mm | |
     | **correct (global beam)** | **1.0** | **52.1 m** | **0.20 mm** | **1.2 mm** | **about 46.5 m** |
     | correct (global beam) | 1.02 | 53.0 m | 0.92 mm | 1.1 mm | |
     | correct (global beam) | 1.1505 | 60.0 m | 1.78 mm | 1.1 mm | |

   - **The analytic Laue radius focuses on the screen with Fcurv = 1.0.** The script
     needed 1.15 to compensate for the bug.
   - **Vertical focus:**
     - With the correct geometry, the Bragg radius RmB puts the vertical focus upstream
       of the screen, around 46.5 m, with 0.95 mm FWHM.
     - The vertical profile is broad, so the size at the screen is only slightly larger
       (1.2 mm) and still smaller than in the script's geometry.
     - RmB was not adjusted.
   - **Flux changes by about 2%** between the two geometries.
2. **One energy, no optimization.**
   - There is no energy scan, no Fcurv optimization and no caustic search.
   - The Runner outputs flux, sizes and energy centroid and peak for the screen, but not
     the script's `res_*_energy.txt` table.
3. **Compliance tensors:**
   - Neither crystal in the script passes `ST=`, so both read `S_tensor.txt`, which holds
     the Laue tensor (nrm [5,−1,−1], sag [0,1,−1]).
   - The workspace reproduces this by giving both widgets that nrm and sag. For the actual
     Bragg crystal orientation, change these fields in the Johann widget.
4. **Laue crystal `aLrad`:**
   - The Bent Crystal 2D widget passes −48.97° to the crystal material. The script passes
     π/2 − aLrad (41.03°).
   - In `bentSi_Crys` this angle only selects an absorption threshold, so it has no
     effect.
5. **Elements not included:**
   - The screen at the mask (`SCR1`) is omitted.
   - The diamond and SiC filter materials are defined in the script but never placed in
     the beamline, so they are omitted too.
6. **Bragg angle:**
   - The widgets and the updated script both use xrt's refraction-corrected Bragg angle.
   - The script's results above were computed with optlnls, which differs by about
     1 µrad for Si(311) at 44 keV.

## Verification

- **Settings.** Every workspace property matches an existing widget setting. The Runner
  script generated from the `.ows` builds the chain shine → propagate → reflect (Bragg) →
  reflect (Laue) → expose, passing the global beam at each step.
- **Geometry.** Positions, pitches, radii and asymmetries were taken from the script's
  own `build_beamline` and `energy_scan` at 44 keV.
- **Element-by-element check.** The same element code was verified element by element
  for SSU_BLM (`../oasys_SSU_BLM/INFO.md`).
- **Not tested in the OASYS GUI.** The workspace was checked by loading each widget
  offscreen with the `.ows` properties and generating and running the Runner script.
