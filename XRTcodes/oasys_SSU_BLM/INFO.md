# SSU_BLM in OASYS2-XRT (80 keV)

This folder holds an OASYS2 workspace for the SSU_BLM beamline of `XRT_SSU_BLM.py` at one
energy, 80 keV. The workspace uses OASYS2-XRT widgets only, including two new ones:
Bending Magnet Light Source and Johann Cylinder Crystal.

## Files

| File | Content |
|---|---|
| `SSU_BLM_80keV.ows` | OASYS2 workspace |
| `SSU_BLM_80keV_runner_script.py` | Python script that the XRT Runner widget generates from this workspace (500k rays) |

## Workspace

Six widgets in a chain, each sending XRTData to the next:

| # | Widget | Name | Settings |
|---|---|---|---|
| 0 | Bending Magnet Light Source | `DW` | B0 = 6.6 T, 3 GeV, 0.1 A; σx = 19.3 µm, σz = 1.9 µm; εx = 0.24897 nm·rad, εz = 0.00247 nm·rad; acceptance 1.6 × 0.66 mrad; E = 79.2–80.8 keV (±1%); 500 000 rays |
| 1 | Rectangular Aperture | `Mask1` | at [0, 10000, 0] mm; gap 16 × 6.6 mm |
| 2 | Johann Cylinder Crystal | `Bragg` | Si(311) Bragg, symmetric, t = 10 mm; at [0, 31000, 0]; pitch `np.pi - alpha - thetaB`; Rm = 647668.0188 mm; ML; tensor oriented with nrm [19,1,1], sag [0,1,−1]; limPhysY ±500 mm |
| 3 | Bent Crystal 2D | `Laue` | Si(311) Laue, nrm [19,1,1], sag [0,1,−1], t = 1 mm; alpha `-alpha_hkl_nrm` (−69.02°); pitch `np.pi/2 + alpha - thetaB`; at [0, 31526.4724688, −50]; Rs = 1217.6367 mm, Rm = −76039.4075 mm, **Fcurv = 1.0**; ML, `ML_sign` −1 |
| 4 | Screen | `SCRexp` | at [0, 56000, −50]; plotted with limits x ±200 µm, z ±5000 µm |
| 5 | XRT Runner Python Script | | beamline name `SSU_BLM_80keV`, 1 repetition, scores written to `SSU_BLM_80keV/scores/` |

- **Where the values come from:**
  - The radii and positions are the script's formulas evaluated at 80 keV:
    θB = 47.344 mrad, Bragg–Laue distance 526.47 mm, RsL = 1217.64 mm,
    RmL = 76039.41 mm, RmB = 647668.02 mm.
  - Both crystals use 80 keV as the alignment energy for thetaB.
- **Running it:** the crystal, aperture and screen widgets have automatic execution
  switched on. Clicking "Send Data" on the Bending Magnet widget should run the whole
  chain.
- **Changing the energy:** update the following, using the formulas in
  `XRT_SSU_BLM.py` (`energy_scan`):
  - `eMin` and `eMax`, and both `alignment_energy` fields
  - RmB, Rs and Rm of the Laue crystal
  - the Laue crystal y position, 31000 + 50/tan(2θB)

## Results (80 keV)

The generated Runner script was run outside the GUI with 100 000 rays. The 500k-ray
script in this folder has not been run.

| Quantity | OASYS workspace (Fcurv = 1.0) | Original script (Fcurv = 1.1175) |
|---|---|---|
| Flux at screen (ph/s) | 3.38×10¹² | 3.38×10¹² |
| Horizontal FWHM at screen | 52 µm | 34 µm |
| Vertical FWHM at screen | 3.4 mm | 3.9 mm |
| Energy centroid / peak | 80.022 / 80.029 keV | — |
| Energy FWHM | not output by the Runner | 0.163 keV |
| Horizontal focus | at the screen (56.00 m) | 56.008 m |

- **How the sizes are measured differs.**
  - The workspace values come from the Runner's plot: 1024 bins, `plot.dx` and `plot.dy`.
  - The script's values come from `fwhm_from_hist` on 200-bin plots over x = ±0.2 mm and
    z = ±10 mm.
- **Same measurement, same geometry.** Measured both ways on the same correct
  (global-beam) geometry at 80 keV with 200k rays, the horizontal FWHM is 40 µm and the
  vertical 2.9 mm (see `../BUG_GEOMETRY.md`).
- **Flux** agrees with the script within statistics.
- **Reflectivity:** both use the internal ML Laue reflectivity, which is too high (see
  `../BUG_ML.md`).

## Differences from `XRT_SSU_BLM.py`

1. **Geometry: global beam after the mask (the main difference).**
   - The script passes the mask's local beam to the Bragg crystal. That effectively moves
     the crystals 10 m farther from the source (`../BUG_GEOMETRY.md`).
   - The OASYS Runner keeps the global beam after an aperture, which is correct.
   - As a result, the analytic Laue radius focuses on the screen with **Fcurv = 1.0**.
     The script needed Fcurv = 1.1175 at 80 keV to compensate for the bug.
   - With the script's Fcurv = 1.1175, the workspace focuses at 61.75 m and the
     horizontal spot at the screen is about 0.4 mm.
2. **One energy, no optimization.**
   - There is no energy scan, no Fcurv optimization and no caustic search: Fcurv is
     entered by hand.
   - The Runner outputs flux, sizes and energy centroid and peak for the screen, but not
     the script's `res_*_energy.txt` table: no flux0, no divergences, no energy FWHM.
3. **Bragg angle:**
   - The crystal widgets compute thetaB with xrt
     (`get_Bragg_angle - get_dtheta_symmetric_Bragg`).
   - The results above come from the script run with optlnls
     `calc_Bragg_angle(corrected=True)`. That differs by 0.6 µrad at 80 keV, which shifts
     the energy centroid by about 1 eV. Reflected intensity changes by 0.004%.
   - The scripts now use the xrt formula too, so new runs use the same angle.
4. **Laue crystal `aLrad`:**
   - The Bent Crystal 2D widget passes the same alpha to the optical element and to the
     crystal material (−69.02°). The script gives the material π/2 − aLrad (20.98°).
   - In `bentSi_Crys` this angle only selects an absorption threshold (`aLrad < π`), so
     it has no effect.
5. **Bragg crystal compliance tensor:**
   - The script's Bragg material has no `ST=`, so it reads `S_tensor.txt`, which holds
     the Laue tensor (nrm [19,1,1], sag [0,1,−1]).
   - The workspace reproduces this with the same nrm and sag in the Johann widget. For
     the actual Bragg crystal orientation, change these two fields.
6. **Elements not included:**
   - The screen at the mask (`SCR1`) is omitted.
   - The diamond and SiC filter materials are defined in the script but never placed in
     the beamline, so they are omitted too.
7. **BeamLine alignment:** the script creates `BeamLine(alignE='auto')` and sets
   `alignE` per energy, while the Runner uses `BeamLine()`. All pitches are numeric, so
   this has no effect.

## Verification

- **Element-by-element check.** The workspace's elements were traced with the same
  source rays as the script's elements (mask, Bragg, Laue, screen) at Fcurv = 1.1175.
  Ray counts and intensities agree within 0.1%, with identical distributions. The Johann
  Cylinder widget alone reproduces the script's Bragg crystal to 0.004% in reflected
  intensity.
- **Origin of the different Fcurv.** The script's own `build_beamline` and
  `energy_scan` at 80 keV, traced with the global beam, focus on the screen with
  Fcurv = 1.0 and at 61.75 m with Fcurv = 1.1175. This is the same as the workspace,
  so the remaining difference is entirely the local-beam bug.
- **Not tested in the OASYS GUI.** The workspace was checked by loading each widget
  offscreen with the `.ows` properties and generating and running the Runner script.
