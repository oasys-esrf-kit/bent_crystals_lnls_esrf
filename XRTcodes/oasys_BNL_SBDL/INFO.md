# BNL_SBDL in OASYS2-XRT (40 keV)

This folder holds an OASYS2 workspace for the BNL_SBDL beamline of `XRT_BNL_SBDL.py` at one
energy, 40 keV, where the script's flux peaks. The workspace uses OASYS2-XRT widgets only.

## Files

| File | Content |
|---|---|
| `BNL_SBDL_40keV.ows` | OASYS2 workspace |
| `BNL_SBDL_40keV_runner_script.py` | Python script that the XRT Runner widget generates from this workspace (500k rays) |

## Workspace

Nine widgets in a chain, each sending XRTData to the next:

| # | Widget | Name | Settings |
|---|---|---|---|
| 0 | Wiggler Light Source | `DW100` | K = 16.5, period 100 mm, 68 periods, 3 GeV, 0.1 A; σx = 137 µm, σz = 4.9 µm; εx = 0.9042 nm·rad, εz = 0.00784 nm·rad; acceptance 1.1 × 0.1 mrad; E = 39.92–40.08 keV (±0.2%); 500 000 rays |
| 1 | Rectangular Aperture | `Mask1` | at [0, 10000, 0]; gap 11000 × 1000 mm (as in the script, see below) |
| 2 | Plate | `Filter1` | diamond, `Material('C', rho=3.52, kind='plate')`, t = 2 mm, at [0, 10000, 0], pitch π/2 |
| 3 | Plate | `Filter2` | SiC, `Material(elements=['Si', 'C'], quantities=[1, 1], kind='plate', rho=3.21)`, t = 4.8 mm, at [0, 10004, 0], pitch π/2 |
| 4 | Bent Crystal 2D | `LaueC1` | Si(111) Laue, nrm [1,0,0], sag [0,−1,1], t = 0.7 mm; alpha `alpha_hkl_nrm` (35.26°); pitch `np.pi/2 + alpha + thetaB`; at [0, 32300, 0]; Rs = −1575.6831 mm, Rm = 50670.4309 mm, **Fcurv = 1.02**; ML, `ML_sign` +1 |
| 5 | Bent Crystal 2D | `LaueC2` | same crystal; alpha `-alpha_hkl_nrm`; pitch `np.pi/2 + alpha + thetaB`; at [0, 32803.8652734, 50]; same radii and Fcurv; ML, `ML_sign` −1 |
| 6 | Toroid Mirror Distorted | `M1` | Pt, `Material('Pt', rho=19.305)`; at [0, 40100, 50]; pitch 1.74886 mrad; R = 13253444.23 mm, r = 1e30 (flat sagittally); 50 × 1300 mm; no distortion file |
| 7 | Screen | `SCRexp` | at [0, 56400, 107.0127]; plotted with limits x ±1500 µm, z ±500 µm |
| 8 | XRT Runner Python Script | | beamline name `BNL_SBDL_40keV`, 1 repetition, scores written to `BNL_SBDL_40keV/scores/` |

- **Where the values come from:**
  - Taken from the script's own `energy_scan` at 40 keV: θB = 49.45 mrad (xrt,
    refraction-corrected), C1–C2 distance 503.87 mm.
  - Rs = 4 sin θB sin α · fs, and Rm = Rs / (c_fac · Poisson ratio), with c_fac = 0.5.
  - M1 pitch from the script's polynomial in energy. Analytic focusing radius
    R = 2 / (1/40.1 m + 1/16.3 m) / sin(pitch).
  - The screen height follows the M1 deflection: 50 + 16.3 m · 2 · pitch = 107.01 mm.
- **Running it:** all optical elements have automatic execution switched on. Clicking
  "Send Data" on the Wiggler widget should run the whole chain.
- **Changing the energy:** update the following, using the formulas in
  `XRT_BNL_SBDL.py` (`energy_scan`):
  - `eMin` and `eMax`, and the two `alignment_energy` fields
  - Rs and Rm
  - the C2 position, 32300 + 50/tan(2θB)
  - the M1 pitch and R
  - the screen height

## Results (40 keV)

The generated Runner script was run outside the GUI with 200 000 rays. The 500k-ray
script in this folder has not been run.

| Quantity | OASYS workspace (Fcurv = 1.02) | Original script (Fcurv = 1.199) |
|---|---|---|
| Flux at screen (ph/s) | **1.04×10¹³** | 5.24×10¹² |
| Horizontal FWHM at screen | 0.50 mm | 0.275 mm |
| Vertical FWHM at screen | 0.123 mm | 0.126 mm |
| Energy centroid / peak | 39.993 / 39.994 keV | — |
| Energy FWHM | not output by the Runner | 46 eV |

- **Measurement:**
  - The workspace values come from the Runner's plot: 1024 bins, `plot.dx` and `plot.dy`.
  - The script's values come from `fwhm_from_hist` on its own plots.
- **Flux is twice the script's.** In the script's geometry the mirror sits 10 m farther
  from the source (see below). The horizontal beam is then wider at M1, and more of it is
  cut by the 50 mm mirror width.
- **The horizontal spot is larger.** With the source at its true distance, the
  horizontal magnification of the double-Laue focusing is about 0.75 instead of 0.57.
  The horizontal source is large (σx = 137 µm), so the image grows accordingly.
- **Reflectivity:** both use the internal ML Laue reflectivity, which is too high (see
  `../BUG_ML.md`).

## Differences from `XRT_BNL_SBDL.py`

1. **Geometry: global beam after the mask (the main difference).**
   - The script passes the local beam of `Screen1` (at 10 m) to the filters. That
     effectively moves the filters, crystals and mirror 10 m farther from the source
     (`../BUG_GEOMETRY.md`).
   - The OASYS Runner keeps the global beam, which is correct.
   - Effects:
     - **Flux about 2× higher**, because less is clipped at M1.
     - **Fcurv ≈ 1.02 instead of 1.199** (1.18–1.20 over the script's energy range).
   - In the correct geometry at 40 keV, the horizontal focus moves with Fcurv:

     | Fcurv | Horizontal focus |
     |---|---|
     | 1.0 | 55.6 m |
     | 1.02 | 56.5 m (the screen is at 56.4 m) |
     | 1.05 | 57.7 m |
     | 1.199 (the script's value) | 63.0 m, 9.6 mm FWHM at the screen |

   - The analytic M1 radius puts the vertical focus at 55.9–56.0 m, with 0.114 mm FWHM
     at the screen. That is essentially the minimum width, so the radius was not
     adjusted. The script's own R scan is replaced by the analytic value.
2. **One energy, no optimization.**
   - There is no energy scan, no Fcurv optimization, no M1 radius scan and no caustic
     search.
   - The Runner outputs flux, sizes and energy centroid and peak for the screen, but not
     the script's `res_*_energy.txt` table: no flux0, no divergences, no energy FWHM.
3. **Mask size.**
   - The script computes the mask as `d_slit1 * acch` without the `1e-3` used in the
     other scripts. Since `acch` is in mrad, this gives 11000 × 1000 mm, i.e. no mask.
   - The workspace keeps this size, to reproduce the script. The beam is limited by the
     source acceptance (±1.1 × ±0.1 mrad).
   - The intended mask is probably 11 × 1 mm (1.1 × 0.1 mrad at 10 m). It would cut the
     acceptance to ±0.55 × ±0.05 mrad and reduce the flux.
4. **Filter size:** the Plate widget has no `limPhysX`/`limPhysY` fields, so the filters
   are ±1000 mm (xrt default) instead of the mask size. Given the mask size, this makes no
   difference.
5. **M1:**
   - The Toroid Mirror Distorted widget, without a distortion file, generates xrt's
     `ToroidMirror` with the same material, radii, size and pitch.
   - `limOptX`/`limOptY` and `alarmLevel` are not passed. In the script the optical
     limits equal the physical ones, and `alarmLevel` only controls warnings.
6. **Elements not included:** the screen at the mask (`Screen1`).
7. **Bragg angle:** the widgets and the script, after the optlnls replacement, both use
   xrt's refraction-corrected Bragg angle. The script's results above were computed with
   optlnls, which differs by less than 1 µrad for Si(111).

## Verification

- **Settings.** Every workspace property matches an existing widget setting. The Runner
  script generated from the `.ows` builds the chain shine → propagate →
  double_refract × 2 → reflect × 3 → expose, passing the global beam at each step.
- **Geometry.** Positions, pitches, radii and asymmetries were taken from the script's
  own `build_beamline` and `energy_scan` at 40 keV.
- **Effect of the local beam.** Tracing the script's own elements at 40 keV with the
  local `Screen1` beam (script) and with the global beam (correct) gives the intensity
  and focus differences listed above.
- **Not tested in the OASYS GUI.** The workspace was checked by loading each widget
  offscreen with the `.ows` properties and generating and running the Runner script.
