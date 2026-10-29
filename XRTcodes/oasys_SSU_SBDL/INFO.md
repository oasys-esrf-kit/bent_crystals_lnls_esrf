# SSU_SBDL in OASYS2-XRT (50 keV)

This folder holds an OASYS2 workspace for the SSU_SBDL beamline of `XRT_SSU_SBDL.py` at one
energy, 50 keV, where the script's flux peaks (50–55 keV). The workspace uses OASYS2-XRT
widgets only.

## Files

| File | Content |
|---|---|
| `SSU_SBDL_50keV.ows` | OASYS2 workspace |
| `SSU_SBDL_50keV_runner_script.py` | Python script that the XRT Runner widget generates from this workspace (500k rays) |

## Workspace

Seven widgets in a chain, each sending XRTData to the next:

| # | Widget | Name | Settings |
|---|---|---|---|
| 0 | Bending Magnet Light Source | `SWLS` | B0 = 6.6 T, 3 GeV, 0.1 A; σx = 19.3 µm, σz = 1.9 µm; εx = 0.24897 nm·rad, εz = 0.00247 nm·rad; acceptance 1.6 × 0.66 mrad; E = 47.5–52.5 keV (±5%); 500 000 rays |
| 1 | Rectangular Aperture | `Mask1` | at [0, 10000, 0]; gap 16 × 6.6 mm |
| 2 | Plate | `Filter1` | diamond, `Material('C', rho=3.52, kind='plate')`, t = 0.05 mm, at [0, 10000, 0], pitch π/2 |
| 3 | Bent Crystal 2D | `LaueC1` | Si(111) Laue, nrm [1,0,0], sag [0,1,1], t = 1 mm; alpha `alpha_hkl_nrm` (35.26°); pitch `np.pi/2 + alpha + thetaB`; at [0, 31000, 0]; Rs = −1263.9563 mm, Rm = 40645.9958 mm, **Fcurv = 0.985**; **PP**, `ML_sign` +1 |
| 4 | Bent Crystal 2D | `LaueC2` | same crystal; alpha `-alpha_hkl_nrm`; pitch `np.pi/2 + alpha + thetaB`; at [0, 31630.6672381, 50]; same radii and Fcurv; PP, `ML_sign` −1 |
| 5 | Screen | `SCRexp` | at [0, 56000, 50]; plotted with limits x ±200 µm, z ±15000 µm |
| 6 | XRT Runner Python Script | | beamline name `SSU_SBDL_50keV`, 1 repetition, scores written to `SSU_SBDL_50keV/scores/` |

- **Where the values come from:**
  - Taken from the script's own `energy_scan` at 50 keV: θB = 39.56 mrad (xrt,
    refraction-corrected), C1–C2 distance 630.67 mm.
  - Rs = 4 sin θB sin α · fs, with p = 31 m and q = 25 m.
  - Rm = Rs / (c_fac · poL), with c_fac = 0.5.
- **Running it:** all optical elements have automatic execution switched on. Clicking
  "Send Data" on the Bending Magnet widget should run the whole chain.
- **Changing the energy:** update the following, using the formulas in
  `XRT_SSU_SBDL.py` (`energy_scan`):
  - `eMin` and `eMax`, and the two `alignment_energy` fields
  - Rs and Rm
  - the C2 position, 31000 + 50/tan(2θB)
  - possibly Fcurv; after the fixes the script's own Fcurv falls from 0.991 to 0.938 over
    35–200 keV

## Results (50 keV)

The generated Runner script was run outside the GUI with 200 000 rays, before and after
the ML fix. The 500k-ray script in this folder has not been run.

| Quantity | Workspace, current code (Fcurv = 0.985) | Script after the fixes (Fcurv = 0.986) | Workspace before the ML fix | Script before the fixes (Fcurv = 1.106) |
|---|---|---|---|---|
| Flux at screen (ph/s) | **4.85×10¹²** | 4.95×10¹² | 4.86×10¹² | 3.74×10¹² |
| Horizontal FWHM at screen | 65 µm | 34 µm | 68 µm | 37.5 µm |
| Vertical FWHM at screen | 14.3 mm | 12.4 mm | 13.1 mm | 16.1 mm (not focused vertically) |
| Energy FWHM | not output by the Runner | 89 eV | — | 115 eV |

- **The ML fix does not change this beamline:** both crystals use PP. The workspace and
  the script after the fixes agree within 2% in flux.

- **Measurement:**
  - The workspace values come from the Runner's plot: 1024 bins, `plot.dx` and `plot.dy`.
  - The script's values come from `fwhm_from_hist` on its own plots.
  - In a direct trace with 200-bin histograms over x = ±0.2 mm, the correct geometry
    gives 46 µm at the screen (34 µm at the focus, 55.95 m).
- **Flux is about 30% higher than the script's before the fixes.** The same difference appears when the
  script's own elements are traced with the global beam (intensity 742 vs 554 at
  Fcurv = 1.0 vs the script's geometry; see below).
- **Reflectivity:** both use the internal PP Laue reflectivity, which is within about
  10% of XOPPY (see `../RUN_XRT.md`, "Laue reflectivity check"). The ML bug
  (`../BUG_ML.md`) does not apply here.

## Differences from `XRT_SSU_SBDL.py`

1. **Geometry: global beam after the mask (the main difference).**
   - The script passes the mask's local beam to the filter. That effectively moves the
     filter and both crystals 10 m farther from the source (`../BUG_GEOMETRY.md`).
   - The OASYS Runner keeps the global beam, which is correct.
   - The table compares the two at 50 keV, with the same 200k source rays (horizontal
     FWHM from 200-bin histograms over ±1 mm, then ±0.2 mm for the fine scan):

     | | Fcurv | Σ intensity | Horizontal focus | x FWHM at 56 m |
     |---|---|---|---|---|
     | script (local beam) | 1.107 (the script found 1.106 at 50 keV) | 554 | 56.1 m | 0.05 mm |
     | correct (global beam) | 0.98 | 727 | 55.7 m | 0.24 mm |
     | **correct (global beam)** | **0.985** | **730** | **55.95 m** | **0.046 mm** |
     | correct (global beam) | 0.99 | 734 | 56.2 m | 0.16 mm |
     | correct (global beam) | 1.0 | 742 | 56.6 m | 0.58 mm |
     | correct (global beam) | 1.107 | 821 | 61.9 m | 1.9 mm |

   - **In the correct geometry the analytic double-Laue radius needs only a small
     correction (Fcurv ≈ 0.985).** The script needed about 1.11 to compensate for the
     bug.
   - **The focus is steep:** a change of 0.005 in Fcurv moves it by about 0.25 m, and
     the size at the screen grows by a factor of 3.
   - **Intensity at the screen is about 30% higher** in the correct geometry.
2. **One energy, no optimization.**
   - There is no energy scan, no Fcurv optimization and no caustic search.
   - The Runner outputs flux, sizes and energy centroid and peak for the screen, but not
     the script's `res_*_energy.txt` table.
3. **Filter size:** the Plate widget has no `limPhysX`/`limPhysY` fields, so the filter is
   ±1000 mm (xrt default) instead of the mask size. Since the mask is in front, this makes
   no difference.
4. **Elements not included:**
   - The screens at the mask and after the filter (`SCR1`, `SCR2`) are omitted.
   - The silicon filter, commented out in the script, is not included.
5. **Bragg angle:**
   - The widgets and the updated script both use xrt's refraction-corrected Bragg angle.
   - The script's results above were computed with optlnls, which differs by less than
     1 µrad for Si(111) at 50 keV.
6. **Crystal material asymmetry:** the widget passes the same alpha to the optical
   element and to the crystal material (+α for C1, −α for C2). That is exactly what the
   script does for this beamline.

## Verification

- **Settings.** Every workspace property matches an existing widget setting. The Runner
  script generated from the `.ows` builds the chain shine → propagate → double_refract →
  reflect (C1) → reflect (C2) → expose, passing the global beam at each step. Both
  crystals use `method='PP'`.
- **Geometry.** Positions, pitches, radii and asymmetries were taken from the script's
  own `build_beamline` and `energy_scan` at 50 keV.
- **Not tested in the OASYS GUI.** The workspace was checked by loading each widget
  offscreen with the `.ows` properties and generating and running the Runner script.
