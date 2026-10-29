# XRT runs: bent-crystal monochromators

The four scripts `XRT_APS_BLM.py`, `XRT_BNL_SBDL.py`, `XRT_SSU_BLM.py` and
`XRT_SSU_SBDL.py` (plus `XRT_SSU_BLM_PP.py`, SSU_BLM with the PP Laue model) were run with
xrt 2.0.0 (`~/venv/oasys2`), 500 000 rays, full energy lists. They have not been compared
with xrt 1.x (not installed in this environment).

There are two sets of runs:

- **Current runs (2026-10-28), in `run_XRT_*_fixed/`:** made after fixing the geometry bug
  in the scripts (`BUG_GEOMETRY.md`) and three bugs in the ML reflectivity of
  `BentCrys2D.py` (`BUG_ML.md`).
  - **The double-Laue beamlines (SSU_SBDL, BNL_SBDL) give valid results:** the simulated
    monochromator works (crystals aligned, focus found), so the numbers can be quoted.
  - **The Bragg+Laue beamlines (SSU_BLM, SSU_BLM_PP, APS_BLM) do not:** the flux falls to
    10⁻⁴–10⁻³ of the old values because the Bragg and Laue reflection windows no longer
    overlap (see [Bragg+Laue beamlines](#bragglaue-beamlines-the-reflection-windows-do-not-overlap)).
    This is an alignment artefact of the simulation, not a property of the beamlines.
  - Two test runs (`run_XRT_SSU_BLM_test_fc0_rock/`, `run_XRT_APS_BLM_test_fc0_rock/`) show
    what these beamlines give once the windows are aligned. They use changes that are not
    yet in the code; see the [TODO list](#todo-decisions-pending).
- **Old runs, in `run_XRT_*/`:** made before the fixes. They are kept for reference in
  [Results before the fixes](#results-before-the-fixes). Their geometry was wrong for all
  four beamlines, and their ML flux was wrong.

The comparison of old and new results, energy by energy, is in `RERUN_RESULTS.md`. The OASYS
workspaces (`oasys_*/INFO.md`) were rerun with the fixed crystal code and agree with the
scripts after the fixes.

## How the runs were made

- **xrt 2.0 compatibility:** `BentCrys2D.py` was adapted to xrt 2.0. The scripts were written for
  xrt 1.x with a patched `oes_base.py` that passed the ray directions and the
  normals to `get_amplitude`. On xrt 2.0, `BentCrys2D.py` now records these values
  during the reflection and passes them to `bentSi_Crys`, including for the Bragg
  crystal on xrt's `JohannCylinder`. The reflectivity physics is unchanged, and the
  file still works with xrt 1.x.
- **Separate folders:** each script ran in its own folder, `run_<script>/`. `XRT_APS_BLM` and
  `XRT_SSU_BLM` do not pass `ST=` to some crystals, so these re-read
  `S_tensor.txt` from the working folder. Parallel runs in one folder would
  overwrite each other's tensor.
- **No display:** the runs used `MPLBACKEND=Agg`. The only warnings in the logs are matplotlib saying
  `plt.show()` cannot open a window.
- **Folders of the current runs:** `run_XRT_<script>_fixed/`, run one after another with
  `/usr/bin/time -v` (`time.log`). The test runs are in `run_XRT_<script>_test_fc0_rock/`.
- **Outputs:** each folder contains `res_*_energy.txt`, the png images for each energy,
  `S_tensor.txt` and `run.log`.
- **Header fix:** in `XRT_SSU_BLM.py` and `XRT_SSU_SBDL.py`, the
  `dE` column was missing from the result-file header. It has now been added, and
  the headers of the existing result files were corrected.

## What each script computes

Each script makes an energy scan of a sagittally focusing bent-crystal
monochromator. At every energy:

1. **Crystal setup:**
   - The script computes the Bragg angle, refraction-corrected, with xrt:
     `get_Bragg_angle(E) - get_dtheta_symmetric_Bragg(E)`. The old runs used optlnls
     `calc_Bragg_angle(corrected=True)`, which differs by at most 1.4 µrad (Si(311) at
     35 keV).
   - It places the second crystal at the spacing that gives a fixed exit offset (±50 mm).
   - It computes the sagittal radius Rs from the thin-lens formula.
   - It computes the meridional radius Rm from the anticlastic ratio,
     Rm = Rs / (c_fac · Poisson ratio).
2. **Fcurv optimization:** a reduced beam of 50k rays is traced with the ML model.
   Fcurv, a scale factor on the bending, is adjusted until the horizontal focus
   lands on the sample screen. The search is broad at the first energy and fine after that.
3. **Full trace:** the full beam is traced with the final reflectivity model (ML,
   or PP for SSU_SBDL). The script then scans along the beam to find where the
   horizontal and vertical foci actually are.
4. **Output:** one row per energy is written to `res_*_energy.txt`, and the beam
   images are saved.

Columns of `res_*_energy.txt`:

| Column | Meaning |
|---|---|
| `flux0` | flux before the monochromator, within the simulated band of ±δE (ph/s) |
| `flux` | flux at the sample screen (ph/s) |
| `flux_density` | flux at the sample screen per mm² (ph/s/mm²) |
| `beam_x`, `beam_z` | FWHM beam size at the sample screen (mm) |
| `beam_dx`, `beam_dz` | FWHM divergence at the sample screen (mrad) |
| `dE` | FWHM energy bandwidth at the sample screen (keV) |
| last columns | the bending radii used (mm) |

## The four beamlines

| Script | Source | Monochromator | Sample screen | Energies |
|---|---|---|---|---|
| APS_BLM | 7 GeV bending magnet, 0.6 T | Si(311) Bragg (Johann cylinder) at 26.5 m + bent Laue, 0.6 mm thick | 52 m | 30–60 keV (16 points) |
| SSU_BLM | 3 GeV bending magnet, 6.6 T | Si(311) Bragg at 31 m + bent Laue, 1 mm thick | 56 m | 35–200 keV (34 points) |
| SSU_SBDL | same source as SSU_BLM, 50 µm diamond filter | two bent Si(111) Laue crystals at 31 m, PP model | 56 m | 35–200 keV (34 points) |
| BNL_SBDL | 3 GeV wiggler, 68 periods, K=16.5, 2 mm diamond + 4.8 mm SiC filters | two bent Si(111) Laue crystals at 32.3 m + Pt mirror M1 at 40.1 m for vertical focusing (its radius is scanned too) | 56.4 m | 30–70 keV (5 points) |

## Main results

| | Peak flux at sample (ph/s) | Horizontal size (µm) | Vertical size | Bandwidth ΔE/E | Fcurv | Status |
|---|---|---|---|---|---|---|
| SSU_SBDL | 4.95×10¹² at 50 keV, falling to 4.7×10⁹ at 200 keV | 34–84 | 7–16 mm (not focused vertically) | 0.16–4.3 % | 0.94–0.99 | current |
| BNL_SBDL | 3.92×10¹² at 50 keV | 404–463 | **0.11–0.13 mm** | 0.06–0.26 % | 1.005–1.030 | current |
| SSU_BLM | 6.7×10¹¹ at 70 keV | 26–60 | 0.5–1.0 mm | 0.03–0.8 % | 0.999–1.002 | test run only |
| APS_BLM | 3.5×10¹¹ at 54 keV | 167–234 | 0.45–1.15 mm | 0.01–0.18 % | 0.996–1.003 | test run only |

- **Fcurv is now close to 1 everywhere.** The analytic radii focus where they should. The
  old values of 1.05–1.20 came from the geometry bug.
- **SSU_SBDL:**
  - Flux is 1.16–1.52× the old values. Part of the gain comes before the crystals: the
    diamond filter is sized like the mask and clipped the misplaced beam of the old runs,
    so `flux0` (measured after the filter) is now 2.2× higher.
  - The horizontal focus is on the screen (55.6–56.3 m, screen at 56 m). It was about 1 m
    upstream before.
  - The bandwidth is 15–36% smaller above 50 keV.
- **BNL_SBDL:**
  - Flux is 0.25× the old value at 30–40 keV and 0.94–1.07× at 50–70 keV. The peak moved
    from 40 to 50 keV.
  - The geometry fix alone gives about 2× the old flux at 40 keV, and the ML fix ×0.12: the
    OASYS workspace, which always had the correct geometry, gives 1.04×10¹³ ph/s before the
    ML fix and 1.29×10¹² after it (`oasys_BNL_SBDL/INFO.md`). **This large drop is not understood yet.**
  - The horizontal spot is larger than before (0.40–0.46 mm instead of 0.24–0.30 mm). Both
    foci are within 0.3 m of the 56.4 m screen.
  - The M1 mirror still gives a 0.11–0.13 mm vertical spot. Peak flux density 8.1×10¹³
    ph/s/mm² at 50 keV (was 1.7×10¹⁴ at 40 keV).
- **SSU_BLM and APS_BLM:** the numbers in the table come from the test runs (see below), not
  from the scripts. They are 0.06–0.21× (SSU_BLM) and 0.16–0.36× (APS_BLM) the old flux.
  They still depend on the open decisions, in particular the absorption coefficient.

### Bragg+Laue beamlines: the reflection windows do not overlap

In `run_XRT_SSU_BLM_fixed/`, `run_XRT_SSU_BLM_PP_fixed/` and `run_XRT_APS_BLM_fixed/` the
flux is 3×10⁷–1.6×10⁹ ph/s, 10⁻⁴–10⁻³ of the old values. The few rays that reach the
screen give meaningless sizes and bandwidths, and the Fcurv search fails at low energy
for SSU_BLM (it stays at its 1.15 starting value).

- **Before the ML fix** the bent Bragg crystal reflected R ≈ 1 over about 330 µrad, because
  absorption had no effect (`BUG_ML.md`). Any Laue alignment fell inside that window.
- **Now** the Bragg reflectivity is narrow, and limited by absorption with depth.
- **The internal model takes the crystal mid-plane as the depth reference.** The
  `PP_center_fraction` argument is stored but ignored. In the thick (10 mm) Bragg crystal
  this shifts the reflection profile by about 160 µrad. XOPPY, and the geometry of the
  ray tracing, use the entrance surface.
- **The scripts align the Laue crystal at the nominal angle and never rock it.** Its
  window therefore misses the band the Bragg crystal reflects.
- The double-Laue beamlines depend only weakly on the depth reference. At one energy each,
  the entrance-surface reference changes the intensity by −20% (SSU_SBDL) and +4%
  (BNL_SBDL).

### Test runs: entrance-surface reference and Laue rocking

The wrapper `test_fc0_rock/run_fc0_rock.py` runs the unmodified `XRT_SSU_BLM.py` and
`XRT_APS_BLM.py` with two changes:

1. **Entrance-surface reference:** `test_fc0_rock/BentCrys2D.py` honours
   `PP_center_fraction` (`test_fc0_rock/PP_center_fraction.diff`), and the wrapper sets it
   to 0 for every crystal. The internal profiles then sit where XOPPY's do, within 0.1 µrad
   (Laue ML at 80 keV: [−0.6, 43.4] µrad against XOPPY's [−0.5, 43.3]).
2. **Rocking the Laue crystal:** before each energy, the Laue pitch is scanned from −60 to
   +160 µrad in 5 µrad steps (100k rays, Fcurv = 1) and set to the offset that maximizes the
   intensity. The script's own Fcurv optimization and full run then follow.

Results:

- **SSU_BLM:** the best offset grows with energy, from 0 at 35 keV to +120 µrad at 200 keV.
  Rocking gains ×1–4 below 60 keV, ×8–14 at 60–85 keV and ×20–30 from 90 keV.
- **APS_BLM:** best offset +5 to +30 µrad, gain ×1–4.
- **Focus:** Fcurv 0.996–1.003. Horizontal focus at 55.9–56.7 m (SSU_BLM, screen 56 m) and
  52.6–53.0 m (APS_BLM, screen 52 m).
- **Vertical size:** 0.5–1.0 mm for SSU_BLM, against 2.3–8.2 mm in the old runs.
- **Flux:** still well below the old values. Most of the remaining factor comes from the
  absorption coefficient in the Bragg ML reflectivity: with xrt's total μ (current) the
  integrated reflectivity at 80 keV is 17 µrad, against 67 µrad in XOPPY (photoelectric μ).
  With XOPPY's μ the internal code gives 71 µrad (`BUG_ML.md`).

`flux/flux0` is mostly the monochromator bandwidth divided by the simulated ±δE band, so
it is not the crystal reflectivity.

### SSU: double Laue (SBDL) vs Bragg + Laue (BLM)

Both use the same source (3 GeV BM, 6.6 T, 1.6 × 0.66 mrad) and the same 56 m sample
screen. SSU_SBDL is the run after the fixes; SSU_BLM is the test run (the only working one),
so its numbers are provisional.

| E (keV) | Flux SBDL | Flux BLM | BLM/SBDL | Flux density BLM/SBDL | ΔE/E SBDL | ΔE/E BLM | z FWHM SBDL | z FWHM BLM |
|---|---|---|---|---|---|---|---|---|
| 35 | 2.68e12 | 6.77e10 | 0.03 | 0.4 | 0.17% | 0.05% | 15.9 mm | 0.8 mm |
| 50 | 4.95e12 | 2.55e11 | 0.05 | 0.6 | 0.18% | 0.08% | 12.4 mm | 0.8 mm |
| 80 | 1.86e12 | 5.90e11 | 0.32 | 4.5 | 0.55% | 0.03% | 9.6 mm | 0.6 mm |
| 100 | 6.48e11 | 3.53e11 | 0.54 | 8.0 | 0.98% | 0.12% | 8.5 mm | 0.7 mm |
| 130 | 1.29e11 | 1.30e11 | 1.01 | 11.7 | 1.74% | 0.30% | 8.4 mm | 0.9 mm |
| 160 | 2.91e10 | 4.60e10 | 1.58 | 19.1 | 2.59% | 0.48% | 7.3 mm | 0.8 mm |
| 200 | 4.75e9 | 1.12e10 | 2.36 | 23.6 | 4.27% | 0.80% | 7.1 mm | 0.9 mm |

- **Flux:** SSU_SBDL is higher up to 130 keV (×14–40 at 35–55 keV, ×3 at 80 keV); SSU_BLM
  is higher above, ×2.4 at 200 keV.
- **Flux density:** SSU_BLM is higher from 60 keV on (×2.5–25), because its Johann Bragg
  crystal focuses vertically (0.5–1 mm against 7–16 mm). Horizontal sizes are similar
  (34–84 µm and 26–60 µm).
- **The bandwidth does not measure the throughput.** `dE` is the energy spread over the
  whole beam at the screen. In SSU_SBDL each vertical angle in the 0.66 mrad acceptance
  selects a different energy, so the spread approaches 0.66 mrad × cot θB (6.7% at 200 keV
  for Si(111); 4.3% FWHM measured). The flux instead follows the throughput per ray: at
  200 keV only about 0.05% of the 4.3% band reaches the screen on average for SSU_SBDL,
  against 0.6% of the 0.8% band for SSU_BLM. So a wider bandwidth with less flux is not a
  contradiction.

**XOPPY integrated reflectivities.** `XOPPY/integrated_SSU.py` computes with `diff_pat` the
integrated s-reflectivity of the three SSU crystals, with the radii of the runs at each
energy (results in `XOPPY/SSU_integrated/integrated_SSU.dat`):

| E (keV) | Bragg Si(311) ML, 10 mm (µrad) | Laue Si(311) PP, 1 mm | Laue Si(111) PP, 1 mm | (L·B)/(L·L) | Flux BLM/SBDL (ray tracing) |
|---|---|---|---|---|---|
| 35 | 9.3 | 7.8 | 33.3 | 0.07 | 0.03 |
| 50 | 18.8 | 20.2 | 62.4 | 0.10 | 0.05 |
| 65 | 37.2 | 30.3 | 80.4 | 0.17 | 0.18 |
| 80 | 67.8 | 33.0 | 80.1 | 0.35 | 0.32 |
| 100 | 128.9 | 29.7 | 67.4 | 0.84 | 0.54 |
| 120 | 152.4 | 24.3 | 53.3 | 1.31 | 0.85 |
| 140 | 174.9 | 19.5 | 41.9 | 1.94 | 1.23 |
| 160 | 194.8 | 15.8 | 33.4 | 2.75 | 1.58 |
| 180 | 211.6 | 12.9 | 27.1 | 3.72 | 1.95 |
| 200 | 224.1 | 10.7 | 22.3 | 4.82 | 2.36 |

- The 1 mm Laue crystals lose reflectivity above 65–80 keV (peak from about 0.9 to
  0.08–0.09 at 200 keV), while the thick bent Bragg crystal keeps a peak near 1 and its
  integrated reflectivity grows with energy. SSU_SBDL pays the Laue loss twice, SSU_BLM
  once.
- The product (L·B)/(L·L) crosses 1 at about 110 keV, the ray tracing at about 130 keV, and
  they agree within a factor of about 2. The product is only a figure of merit: it ignores
  how the two reflectivity curves overlap for each ray.
- XOPPY uses the photoelectric μ. The ray tracing uses xrt's total μ, which lowers the
  Bragg reflectivity (17 instead of 68 µrad at 80 keV in the internal code), consistent
  with the ray-tracing ratio being lower above 100 keV.
- The same ratio in energy units (each integral × cot θB) crosses 1 only at about
  180 keV and agrees less well with the ray tracing.

## Per-energy results

The numbers below are taken from `res_*_energy.txt` in the folders named in each heading.
Sizes are given in µm here (the files use mm). Fcurv and the H/V foci come from `run.log`.
The fixed runs of the Bragg+Laue beamlines are not listed here; their fluxes are in
`RERUN_RESULTS.md`.

### SSU_SBDL (Si111 double Laue, PP): `run_XRT_SSU_SBDL_fixed/`

The script has no vertical focusing and does not search for a vertical focus.

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] | Fcurv | H / V focus [m] |
|---|---|---|---|---|---|---|---|---|---|---|
| 35 | 8.29e14 | 2.68e12 | 3.66e12 | 44 | 15926 | 1.223 | 0.236 | 0.061 | 0.9910 | 55.98 / — |
| 40 | 7.60e14 | 3.79e12 | 6.64e12 | 38 | 13628 | 1.118 | 0.255 | 0.064 | 0.9885 | 55.96 / — |
| 45 | 6.95e14 | 4.60e12 | 8.49e12 | 34 | 14493 | 1.010 | 0.213 | 0.071 | 0.9872 | 55.93 / — |
| 50 | 6.35e14 | 4.95e12 | 1.07e13 | 34 | 12424 | 1.012 | 0.209 | 0.089 | 0.9860 | 55.95 / — |
| 55 | 5.77e14 | 4.74e12 | 1.07e13 | 34 | 11589 | 0.885 | 0.198 | 0.123 | 0.9848 | 55.99 / — |
| 60 | 5.23e14 | 4.31e12 | 7.71e12 | 41 | 12924 | 0.886 | 0.209 | 0.157 | 0.9835 | 55.95 / — |
| 65 | 4.74e14 | 3.65e12 | 7.82e12 | 41 | 10847 | 0.870 | 0.183 | 0.218 | 0.9816 | 55.95 / — |
| 70 | 4.28e14 | 2.95e12 | 5.05e12 | 52 | 10168 | 0.817 | 0.172 | 0.289 | 0.9791 | 55.92 / — |
| 75 | 3.87e14 | 2.41e12 | 5.48e12 | 41 | 9646 | 0.780 | 0.184 | 0.381 | 0.9791 | 55.95 / — |
| 80 | 3.49e14 | 1.86e12 | 4.03e12 | 45 | 9598 | 0.819 | 0.169 | 0.438 | 0.9785 | 55.98 / — |
| 85 | 3.14e14 | 1.45e12 | 2.83e12 | 48 | 9359 | 0.801 | 0.169 | 0.568 | 0.9769 | 55.99 / — |
| 90 | 2.83e14 | 1.10e12 | 2.56e12 | 47 | 8630 | 0.804 | 0.149 | 0.654 | 0.9750 | 55.99 / — |
| 95 | 2.54e14 | 8.65e11 | 1.65e12 | 55 | 9119 | 0.713 | 0.148 | 0.839 | 0.9725 | 55.88 / — |
| 100 | 2.29e14 | 6.48e11 | 1.22e12 | 59 | 8524 | 0.770 | 0.153 | 0.976 | 0.9719 | 55.91 / — |
| 105 | 2.05e14 | 4.90e11 | 1.00e12 | 52 | 8799 | 0.735 | 0.152 | 1.145 | 0.9702 | 55.90 / — |
| 110 | 1.83e14 | 3.70e11 | 6.39e11 | 62 | 8596 | 0.681 | 0.139 | 1.364 | 0.9671 | 55.83 / — |
| 115 | 1.64e14 | 2.88e11 | 5.69e11 | 58 | 8175 | 0.699 | 0.136 | 1.517 | 0.9675 | 55.92 / — |
| 120 | 1.47e14 | 2.17e11 | 3.81e11 | 62 | 8672 | 0.680 | 0.142 | 1.564 | 0.9675 | 55.98 / — |
| 125 | 1.32e14 | 1.70e11 | 2.36e11 | 83 | 8389 | 0.674 | 0.136 | 2.026 | 0.9725 | 56.31 / — |
| 130 | 1.18e14 | 1.29e11 | 2.21e11 | 67 | 8402 | 0.645 | 0.127 | 2.256 | 0.9644 | 56.13 / — |
| 135 | 1.06e14 | 1.00e11 | 2.02e11 | 61 | 7660 | 0.686 | 0.125 | 2.660 | 0.9612 | 55.98 / — |
| 140 | 9.39e13 | 7.62e10 | 1.28e11 | 59 | 8405 | 0.650 | 0.124 | 2.973 | 0.9593 | 55.87 / — |
| 145 | 8.41e13 | 6.05e10 | 1.07e11 | 72 | 7606 | 0.652 | 0.121 | 3.185 | 0.9568 | 55.90 / — |
| 150 | 7.54e13 | 4.77e10 | 7.47e10 | 75 | 8199 | 0.657 | 0.122 | 3.562 | 0.9528 | 55.72 / — |
| 155 | 6.71e13 | 3.67e10 | 6.65e10 | 70 | 7493 | 0.638 | 0.119 | 3.807 | 0.9528 | 55.79 / — |
| 160 | 6.00e13 | 2.91e10 | 4.81e10 | 80 | 7280 | 0.642 | 0.107 | 4.142 | 0.9478 | 55.64 / — |
| 165 | 5.35e13 | 2.32e10 | 4.23e10 | 69 | 7500 | 0.639 | 0.113 | 4.919 | 0.9510 | 55.70 / — |
| 170 | 4.79e13 | 1.82e10 | 3.19e10 | 73 | 7492 | 0.663 | 0.108 | 5.448 | 0.9465 | 55.70 / — |
| 175 | 4.25e13 | 1.45e10 | 2.61e10 | 72 | 7472 | 0.621 | 0.108 | 5.493 | 0.9497 | 55.89 / — |
| 180 | 3.79e13 | 1.16e10 | 1.81e10 | 84 | 7383 | 0.638 | 0.108 | 6.347 | 0.9510 | 56.02 / — |
| 185 | 3.38e13 | 9.34e9 | 1.71e10 | 70 | 7290 | 0.663 | 0.110 | 6.849 | 0.9453 | 55.86 / — |
| 190 | 3.01e13 | 7.38e9 | 1.21e10 | 81 | 7276 | 0.643 | 0.103 | 7.225 | 0.9466 | 55.76 / — |
| 195 | 2.68e13 | 5.99e9 | 1.13e10 | 69 | 7401 | 0.665 | 0.102 | 8.143 | 0.9428 | 55.87 / — |
| 200 | 2.39e13 | 4.75e9 | 8.09e9 | 78 | 7064 | 0.651 | 0.099 | 8.540 | 0.9379 | 55.79 / — |

### BNL_SBDL (Si111 double Laue, ML, + M1): `run_XRT_BNL_SBDL_fixed/`

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] | Fcurv | H / V focus [m] |
|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 1.13e15 | 9.36e11 | 1.77e13 | 463 | 107 | 1.352 | 0.066 | 0.017 | 1.0300 | 56.59 / 56.71 |
| 40 | 4.99e14 | 1.30e12 | 2.82e13 | 440 | 108 | 1.131 | 0.067 | 0.035 | 1.0194 | 56.35 / 56.52 |
| 50 | 2.11e14 | 3.92e12 | 8.14e13 | 442 | 112 | 0.947 | 0.068 | 0.070 | 1.0094 | 56.22 / 56.50 |
| 60 | 8.73e13 | 1.64e12 | 3.09e13 | 429 | 120 | 0.847 | 0.062 | 0.116 | 1.0194 | 56.40 / 56.54 |
| 70 | 3.55e13 | 4.47e11 | 8.80e12 | 404 | 129 | 0.814 | 0.054 | 0.180 | 1.0048 | 56.39 / 56.40 |

### SSU_BLM (Si311, ML), test run: `run_XRT_SSU_BLM_test_fc0_rock/`

Not produced by the script: entrance-surface reference and Laue rocking (see above).

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] | Fcurv | H / V focus [m] | Laue offset [µrad] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 35 | 1.66e14 | 6.77e10 | 1.61e12 | 46 | 800 | 1.886 | 0.096 | 0.016 | 0.9992 | 56.13 / 56.15 | +0 |
| 40 | 1.52e14 | 1.50e11 | 4.67e12 | 30 | 900 | 1.950 | 0.128 | 0.024 | 0.9999 | 56.11 / 55.94 | +15 |
| 45 | 1.39e14 | 3.03e11 | 1.08e13 | 26 | 900 | 1.961 | 0.160 | 0.030 | 0.9999 | 56.30 / 55.93 | +20 |
| 50 | 1.27e14 | 2.55e11 | 6.06e12 | 42 | 800 | 1.703 | 0.192 | 0.038 | 0.9992 | 56.34 / 56.13 | +20 |
| 55 | 1.15e14 | 2.90e11 | 7.10e12 | 44 | 800 | 1.832 | 0.192 | 0.042 | 1.0005 | 56.45 / 56.15 | +25 |
| 60 | 1.05e14 | 5.55e11 | 1.93e13 | 34 | 700 | 1.972 | 0.192 | 0.043 | 0.9999 | 56.65 / 56.25 | +30 |
| 65 | 9.48e13 | 6.60e11 | 2.44e13 | 34 | 600 | 1.713 | 0.160 | 0.038 | 0.9999 | 56.70 / 56.45 | +35 |
| 70 | 8.52e13 | 6.70e11 | 1.95e13 | 38 | 800 | 1.832 | 0.160 | 0.034 | 0.9999 | 56.51 / 56.58 | +40 |
| 75 | 7.73e13 | 6.28e11 | 2.65e13 | 34 | 500 | 1.541 | 0.160 | 0.027 | 0.9999 | 56.57 / 56.35 | +40 |
| 80 | 6.97e13 | 5.90e11 | 1.80e13 | 42 | 600 | 1.627 | 0.160 | 0.023 | 1.0005 | 56.51 / 56.42 | +45 |
| 85 | 6.28e13 | 5.16e11 | 1.63e13 | 36 | 700 | 1.616 | 0.160 | 0.045 | 0.9999 | 56.60 / 56.52 | +50 |
| 90 | 5.64e13 | 4.65e11 | 1.47e13 | 42 | 600 | 1.498 | 0.160 | 0.061 | 1.0005 | 56.55 / 56.41 | +50 |
| 95 | 5.08e13 | 4.03e11 | 1.40e13 | 40 | 600 | 1.509 | 0.160 | 0.093 | 1.0005 | 56.38 / 56.37 | +55 |
| 100 | 4.55e13 | 3.53e11 | 9.75e12 | 40 | 700 | 1.541 | 0.128 | 0.122 | 1.0005 | 56.39 / 56.58 | +60 |
| 105 | 4.08e13 | 3.03e11 | 9.83e12 | 36 | 700 | 1.422 | 0.128 | 0.157 | 1.0005 | 56.19 / 56.46 | +60 |
| 110 | 3.66e13 | 2.60e11 | 7.40e12 | 40 | 700 | 1.390 | 0.128 | 0.180 | 0.9999 | 56.00 / 55.99 | +65 |
| 115 | 3.27e13 | 2.19e11 | 4.79e12 | 42 | 900 | 1.369 | 0.128 | 0.241 | 1.0005 | 56.00 / 56.01 | +70 |
| 120 | 2.93e13 | 1.84e11 | 4.19e12 | 46 | 800 | 1.336 | 0.128 | 0.280 | 1.0005 | 56.01 / 56.16 | +70 |
| 125 | 2.63e13 | 1.54e11 | 3.32e12 | 44 | 900 | 1.347 | 0.128 | 0.341 | 0.9999 | 56.00 / 55.95 | +75 |
| 130 | 2.36e13 | 1.30e11 | 2.59e12 | 50 | 900 | 1.347 | 0.128 | 0.391 | 1.0005 | 56.01 / 55.80 | +80 |
| 135 | 2.10e13 | 1.12e11 | 2.27e12 | 50 | 900 | 1.272 | 0.128 | 0.454 | 1.0011 | 56.07 / 55.96 | +80 |
| 140 | 1.87e13 | 9.37e10 | 1.42e12 | 60 | 1000 | 1.272 | 0.128 | 0.522 | 1.0017 | 56.04 / 56.23 | +85 |
| 145 | 1.67e13 | 7.81e10 | 1.75e12 | 44 | 900 | 1.185 | 0.128 | 0.537 | 1.0005 | 56.02 / 56.13 | +85 |
| 150 | 1.50e13 | 6.66e10 | 1.15e12 | 56 | 900 | 1.185 | 0.128 | 0.658 | 0.9992 | 55.93 / 55.85 | +90 |
| 155 | 1.34e13 | 5.50e10 | 9.80e11 | 54 | 900 | 1.228 | 0.128 | 0.750 | 0.9992 | 55.98 / 56.05 | +95 |
| 160 | 1.19e13 | 4.60e10 | 9.18e11 | 54 | 800 | 1.164 | 0.128 | 0.769 | 1.0024 | 56.10 / 56.17 | +95 |
| 165 | 1.07e13 | 3.90e10 | 6.64e11 | 56 | 900 | 1.207 | 0.128 | 0.911 | 1.0024 | 56.13 / 56.06 | +100 |
| 170 | 9.47e12 | 3.28e10 | 4.57e11 | 56 | 1000 | 1.196 | 0.128 | 0.950 | 1.0024 | 56.10 / 56.05 | +105 |
| 175 | 8.43e12 | 2.72e10 | 5.77e11 | 52 | 800 | 1.153 | 0.128 | 1.053 | 1.0011 | 56.05 / 55.96 | +105 |
| 180 | 7.53e12 | 2.27e10 | 4.65e11 | 54 | 800 | 1.153 | 0.096 | 1.126 | 0.9999 | 56.00 / 55.95 | +110 |
| 185 | 6.67e12 | 1.91e10 | 2.93e11 | 58 | 1000 | 1.142 | 0.128 | 1.206 | 1.0005 | 55.96 / 55.95 | +115 |
| 190 | 5.96e12 | 1.59e10 | 3.07e11 | 52 | 900 | 1.078 | 0.128 | 1.422 | 1.0011 | 55.98 / 55.78 | +115 |
| 195 | 5.34e12 | 1.33e10 | 2.14e11 | 56 | 1000 | 1.142 | 0.096 | 1.422 | 0.9999 | 55.93 / 55.84 | +120 |
| 200 | 4.74e12 | 1.12e10 | 1.91e11 | 56 | 900 | 1.078 | 0.128 | 1.595 | 0.9992 | 55.98 / 55.73 | +120 |

### APS_BLM (Si311, ML), test run: `run_XRT_APS_BLM_test_fc0_rock/`

Not produced by the script: entrance-surface reference and Laue rocking (see above).

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] | Fcurv | H / V focus [m] | Laue offset [µrad] |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 30 | 7.54e13 | 1.76e11 | 1.85e12 | 195 | 450 | 1.008 | 0.068 | 0.002 | 1.0010 | 53.03 / 52.90 | +5 |
| 32 | 6.98e13 | 1.47e11 | 1.67e12 | 178 | 450 | 0.924 | 0.085 | 0.003 | 0.9985 | 52.96 / 52.83 | +5 |
| 34 | 6.42e13 | 1.62e11 | 1.38e12 | 189 | 550 | 1.047 | 0.085 | 0.006 | 0.9972 | 52.77 / 52.76 | +10 |
| 36 | 5.90e13 | 2.63e11 | 3.17e12 | 167 | 450 | 1.042 | 0.068 | 0.008 | 0.9997 | 52.84 / 52.57 | +10 |
| 38 | 5.41e13 | 2.14e11 | 2.15e12 | 189 | 500 | 1.008 | 0.068 | 0.011 | 1.0022 | 52.71 / 52.64 | +10 |
| 40 | 5.00e13 | 1.86e11 | 1.33e12 | 234 | 550 | 0.874 | 0.068 | 0.014 | 0.9957 | 52.58 / 52.64 | +10 |
| 42 | 4.58e13 | 1.70e11 | 1.35e12 | 167 | 600 | 1.036 | 0.068 | 0.018 | 1.0007 | 52.73 / 52.48 | +15 |
| 44 | 4.21e13 | 1.90e11 | 1.20e12 | 206 | 700 | 0.946 | 0.085 | 0.024 | 1.0032 | 52.67 / 52.53 | +15 |
| 46 | 3.86e13 | 1.97e11 | 1.22e12 | 195 | 750 | 1.019 | 0.068 | 0.035 | 1.0007 | 52.63 / 52.55 | +20 |
| 48 | 3.55e13 | 2.91e11 | 1.89e12 | 195 | 700 | 1.042 | 0.068 | 0.042 | 1.0013 | 52.80 / 52.43 | +20 |
| 50 | 3.25e13 | 3.35e11 | 2.19e12 | 189 | 750 | 1.042 | 0.068 | 0.049 | 0.9997 | 52.71 / 52.60 | +20 |
| 52 | 2.97e13 | 3.43e11 | 2.31e12 | 189 | 750 | 1.030 | 0.068 | 0.050 | 1.0003 | 52.69 / 52.51 | +20 |
| 54 | 2.72e13 | 3.46e11 | 1.75e12 | 200 | 900 | 1.058 | 0.068 | 0.068 | 0.9978 | 52.75 / 52.49 | +25 |
| 56 | 2.49e13 | 3.38e11 | 1.84e12 | 184 | 900 | 1.030 | 0.068 | 0.075 | 1.0019 | 52.84 / 52.55 | +25 |
| 58 | 2.28e13 | 3.25e11 | 1.43e12 | 200 | 1000 | 1.058 | 0.068 | 0.091 | 1.0019 | 52.78 / 52.90 | +30 |
| 60 | 2.09e13 | 3.16e11 | 1.37e12 | 184 | 1150 | 1.047 | 0.085 | 0.107 | 1.0000 | 52.91 / 52.59 | +30 |

## Results before the fixes

These are the original runs, in `run_XRT_*/`. **They are superseded:**

- All four beamlines used the wrong geometry: the optics sat effectively 10 m farther
  from the source (`BUG_GEOMETRY.md`). Fcurv, focus positions and sizes are affected.
- The ML reflectivity was wrong: Laue too high, Bragg without absorption (`BUG_ML.md`).
  This affects the flux of SSU_BLM, APS_BLM and BNL_SBDL.

### Main results (old)

| | Peak flux at sample (ph/s) | Horizontal size (µm) | Vertical size | Bandwidth ΔE/E | Fcurv |
|---|---|---|---|---|---|
| APS_BLM | ~1.0×10¹², flat from 38 to 54 keV | 130–160 | 1.25 → 2.5 mm | 0.05–0.35 % | 1.149–1.155 |
| SSU_BLM | 3.4×10¹² at 75–85 keV (ML, suspect; 2.2×10¹² at 70 keV with PP) | 26–64 | 2.3 → 7.9 mm | 0.03–1.4 % | 1.110–1.120 |
| SSU_SBDL | 3.7×10¹² at 50–55 keV, falling to 4×10⁹ at 200 keV | 25–90 | 8–20 mm (not focused vertically) | 0.13–5 % | 1.05–1.11 |
| BNL_SBDL | 5.2×10¹² at 40 keV | 240–295 | **0.11–0.13 mm** | 0.07–0.24 % | 1.18–1.20 |

- **APS_BLM:**
  - Both foci sit about 1 m beyond the sample screen (≈53.0–53.2 m, against 52 m).
  - At the true focus the horizontal width is 12–80 µm. On the screen it is ~140 µm.
- **SSU_BLM:**
  - The horizontal focus is on the screen (55.9–56.1 m), so the spot is narrow (≈30–60 µm).
  - Of the two SSU options, it gives the most flux at high energy.
- **SSU_SBDL:**
  - The horizontal focus falls about 1 m before the screen (≈54.8–55.2 m).
  - The beam is spread vertically over 8–20 mm, because nothing focuses it vertically.
  - Above about 100 keV the flux drops much faster than in SSU_BLM.
- **BNL_SBDL:**
  - The M1 mirror gives the only sub-mm vertical spot, so this beamline has by far the
    highest flux density, up to 1.7×10¹⁴ ph/s/mm².
  - The two Laue crystals together transmit 9–15% of the beam reaching them.
- **Fcurv:** in every case the optimum is above 1. The analytic sagittal radius
  has to be corrected by 5–20% to put the focus at the target distance.

`flux/flux0` ranges from 0.03% to 5%. It is mostly the monochromator bandwidth
divided by the simulated ±δE band, so it is not the crystal reflectivity.

### Per-energy results (old)

The numbers below are taken from `run_XRT_*/res_*_energy.txt` (old runs). Sizes are given in µm here (the files use mm).

#### APS_BLM (Si311, ML)

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] |
|---|---|---|---|---|---|---|---|---|
| 30 | 7.55e13 | 8.06e11 | 3.68e12 | 163 | 1250 | 1.424 | 0.094 | 0.014 |
| 32 | 6.95e13 | 8.89e11 | 4.56e12 | 129 | 1400 | 1.432 | 0.094 | 0.020 |
| 34 | 6.41e13 | 9.30e11 | 4.25e12 | 144 | 1450 | 1.432 | 0.075 | 0.026 |
| 36 | 5.90e13 | 9.68e11 | 3.49e12 | 144 | 1800 | 1.439 | 0.075 | 0.031 |
| 38 | 5.44e13 | 1.02e12 | 3.59e12 | 158 | 1700 | 1.432 | 0.075 | 0.038 |
| 40 | 5.00e13 | 1.01e12 | 3.31e12 | 144 | 1950 | 1.424 | 0.075 | 0.053 |
| 42 | 4.59e13 | 1.04e12 | 3.34e12 | 139 | 2100 | 1.432 | 0.075 | 0.059 |
| 44 | 4.22e13 | 1.05e12 | 3.47e12 | 148 | 1950 | 1.432 | 0.075 | 0.064 |
| 46 | 3.87e13 | 1.02e12 | 3.39e12 | 134 | 2150 | 1.439 | 0.075 | 0.080 |
| 48 | 3.54e13 | 1.03e12 | 2.98e12 | 144 | 2300 | 1.455 | 0.075 | 0.091 |
| 50 | 3.24e13 | 1.02e12 | 2.66e12 | 139 | 2500 | 1.439 | 0.075 | 0.107 |
| 52 | 2.96e13 | 1.00e12 | 2.80e12 | 148 | 2350 | 1.439 | 0.056 | 0.130 |
| 54 | 2.73e13 | 9.95e11 | 2.83e12 | 134 | 2250 | 1.447 | 0.056 | 0.137 |
| 56 | 2.48e13 | 9.56e11 | 2.57e12 | 153 | 2300 | 1.447 | 0.056 | 0.188 |
| 58 | 2.28e13 | 9.38e11 | 2.58e12 | 144 | 2400 | 1.447 | 0.056 | 0.183 |
| 60 | 2.08e13 | 8.75e11 | 2.24e12 | 148 | 2500 | 1.447 | 0.075 | 0.211 |

#### BNL_SBDL (Si111 double Laue, ML, + M1)

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] | Laue transmission |
|---|---|---|---|---|---|---|---|---|---|
| 30 | 1.13e15 | 3.73e12 | 9.88e13 | 295 | 128 | 1.960 | 0.073 | 0.021 | 0.093 |
| 40 | 4.99e14 | 5.24e12 | 1.72e14 | 275 | 126 | 1.599 | 0.071 | 0.046 | 0.123 |
| 50 | 2.11e14 | 3.66e12 | 1.22e14 | 255 | 126 | 1.372 | 0.065 | 0.084 | 0.151 |
| 60 | 8.72e13 | 1.65e12 | 5.59e13 | 260 | 121 | 1.187 | 0.059 | 0.127 | 0.148 |
| 70 | 3.55e13 | 4.75e11 | 2.09e13 | 239 | 109 | 1.052 | 0.050 | 0.166 | 0.102 |

#### SSU_BLM (Si311, ML)

The flux below is suspect, because the internal ML overestimates the Laue reflectivity. The
same run with PP is in [Laue reflectivity check: PP vs ML](#laue-reflectivity-check-pp-vs-ml) (also an old run).

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] |
|---|---|---|---|---|---|---|---|---|
| 35 | 1.66e14 | 1.21e12 | 1.43e13 | 28 | 2400 | 2.642 | 0.130 | 0.010 |
| 40 | 1.52e14 | 1.59e12 | 1.67e13 | 38 | 2300 | 2.642 | 0.130 | 0.013 |
| 45 | 1.39e14 | 1.98e12 | 3.02e13 | 26 | 2300 | 2.656 | 0.163 | 0.016 |
| 50 | 1.27e14 | 2.37e12 | 3.11e13 | 30 | 2400 | 2.656 | 0.163 | 0.022 |
| 55 | 1.15e14 | 2.66e12 | 3.48e13 | 26 | 2500 | 2.656 | 0.163 | 0.031 |
| 60 | 1.05e14 | 2.88e12 | 2.96e13 | 32 | 2800 | 2.670 | 0.163 | 0.047 |
| 65 | 9.48e13 | 3.11e12 | 3.08e13 | 32 | 2900 | 2.670 | 0.163 | 0.067 |
| 70 | 8.55e13 | 3.25e12 | 2.46e13 | 34 | 3700 | 2.684 | 0.195 | 0.090 |
| 75 | 7.73e13 | 3.38e12 | 2.42e13 | 34 | 3800 | 2.684 | 0.195 | 0.120 |
| 80 | 6.98e13 | 3.38e12 | 2.32e13 | 34 | 3900 | 2.684 | 0.195 | 0.163 |
| 85 | 6.28e13 | 3.38e12 | 2.06e13 | 38 | 4100 | 2.684 | 0.228 | 0.190 |
| 90 | 5.65e13 | 3.11e12 | 1.61e13 | 36 | 4900 | 2.699 | 0.228 | 0.254 |
| 95 | 5.07e13 | 2.67e12 | 1.14e13 | 44 | 5000 | 2.699 | 0.260 | 0.329 |
| 100 | 4.55e13 | 2.30e12 | 1.16e13 | 38 | 4800 | 2.713 | 0.260 | 0.387 |
| 105 | 4.09e13 | 2.00e12 | 8.56e12 | 42 | 5200 | 2.699 | 0.293 | 0.423 |
| 110 | 3.66e13 | 1.71e12 | 6.12e12 | 40 | 6200 | 2.713 | 0.293 | 0.558 |
| 115 | 3.28e13 | 1.47e12 | 6.07e12 | 36 | 6100 | 2.713 | 0.325 | 0.646 |
| 120 | 2.93e13 | 1.26e12 | 4.32e12 | 42 | 6400 | 2.713 | 0.325 | 0.780 |
| 125 | 2.62e13 | 1.08e12 | 3.34e12 | 46 | 6500 | 2.713 | 0.358 | 0.858 |
| 130 | 2.34e13 | 9.34e11 | 2.88e12 | 44 | 6700 | 2.713 | 0.358 | 0.999 |
| 135 | 2.10e13 | 7.92e11 | 1.96e12 | 50 | 7000 | 2.713 | 0.358 | 1.141 |
| 140 | 1.87e13 | 6.66e11 | 2.00e12 | 46 | 6500 | 2.670 | 0.358 | 1.308 |
| 145 | 1.67e13 | 5.49e11 | 1.34e12 | 46 | 8200 | 2.500 | 0.358 | 1.375 |
| 150 | 1.50e13 | 4.55e11 | 1.09e12 | 52 | 7600 | 2.414 | 0.325 | 1.664 |
| 155 | 1.33e13 | 3.75e11 | 9.30e11 | 52 | 6500 | 2.358 | 0.325 | 1.809 |
| 160 | 1.19e13 | 3.14e11 | 7.63e11 | 50 | 7800 | 2.244 | 0.325 | 2.000 |
| 165 | 1.06e13 | 2.58e11 | 6.04e11 | 54 | 7400 | 2.102 | 0.325 | 2.145 |
| 170 | 9.47e12 | 2.14e11 | 4.73e11 | 58 | 7500 | 2.074 | 0.325 | 2.159 |
| 175 | 8.45e12 | 1.75e11 | 4.03e11 | 56 | 7200 | 1.946 | 0.325 | 2.292 |
| 180 | 7.53e12 | 1.43e11 | 2.99e11 | 60 | 7800 | 1.889 | 0.325 | 2.376 |
| 185 | 6.69e12 | 1.18e11 | 2.57e11 | 58 | 7600 | 1.818 | 0.325 | 2.461 |
| 190 | 5.98e12 | 9.81e10 | 2.38e11 | 52 | 7600 | 1.761 | 0.325 | 2.546 |
| 195 | 5.32e12 | 8.02e10 | 1.74e11 | 58 | 7700 | 1.704 | 0.293 | 2.691 |
| 200 | 4.73e12 | 6.59e10 | 1.27e11 | 64 | 7900 | 1.591 | 0.358 | 2.860 |

#### SSU_SBDL (Si111 double Laue, PP)

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] |
|---|---|---|---|---|---|---|---|---|
| 35 | 3.69e14 | 1.77e12 | 2.71e12 | 30 | 19839 | 1.203 | 0.230 | 0.044 |
| 40 | 3.49e14 | 2.62e12 | 5.65e12 | 30 | 14535 | 1.088 | 0.194 | 0.057 |
| 45 | 3.25e14 | 3.33e12 | 6.84e12 | 25 | 14357 | 1.161 | 0.217 | 0.079 |
| 50 | 3.01e14 | 3.74e12 | 6.02e12 | 38 | 16061 | 1.186 | 0.212 | 0.115 |
| 55 | 2.77e14 | 3.73e12 | 7.95e12 | 34 | 12595 | 1.057 | 0.186 | 0.178 |
| 60 | 2.53e14 | 3.42e12 | 4.70e12 | 45 | 15080 | 0.968 | 0.203 | 0.247 |
| 65 | 2.32e14 | 3.07e12 | 5.89e12 | 36 | 13481 | 1.015 | 0.192 | 0.330 |
| 70 | 2.10e14 | 2.51e12 | 4.08e12 | 45 | 13096 | 1.003 | 0.184 | 0.413 |
| 75 | 1.90e14 | 1.97e12 | 3.16e12 | 53 | 11366 | 0.858 | 0.167 | 0.554 |
| 80 | 1.72e14 | 1.56e12 | 3.06e12 | 42 | 11608 | 0.956 | 0.157 | 0.660 |
| 85 | 1.55e14 | 1.22e12 | 1.88e12 | 58 | 10757 | 0.794 | 0.163 | 0.812 |
| 90 | 1.40e14 | 9.16e11 | 1.85e12 | 50 | 9365 | 0.896 | 0.166 | 1.014 |
| 95 | 1.26e14 | 7.27e11 | 1.30e12 | 48 | 10796 | 0.819 | 0.153 | 1.149 |
| 100 | 1.13e14 | 5.35e11 | 1.04e12 | 48 | 9995 | 0.858 | 0.138 | 1.408 |
| 105 | 1.02e14 | 4.21e11 | 9.77e11 | 41 | 9751 | 0.835 | 0.133 | 1.647 |
| 110 | 9.16e13 | 3.09e11 | 5.07e11 | 56 | 10550 | 0.780 | 0.133 | 1.859 |
| 115 | 8.20e13 | 2.42e11 | 5.01e11 | 45 | 9997 | 0.801 | 0.136 | 2.268 |
| 120 | 7.33e13 | 1.85e11 | 2.90e11 | 64 | 9679 | 0.760 | 0.134 | 2.461 |
| 125 | 6.58e13 | 1.41e11 | 2.77e11 | 55 | 8773 | 0.752 | 0.127 | 2.831 |
| 130 | 5.87e13 | 1.11e11 | 1.79e11 | 69 | 8775 | 0.751 | 0.127 | 2.916 |
| 135 | 5.28e13 | 8.43e10 | 1.67e11 | 53 | 9047 | 0.775 | 0.121 | 3.433 |
| 140 | 4.71e13 | 6.53e10 | 1.17e11 | 61 | 7730 | 0.753 | 0.109 | 3.680 |
| 145 | 4.19e13 | 5.16e10 | 7.67e10 | 73 | 8745 | 0.753 | 0.116 | 4.146 |
| 150 | 3.76e13 | 4.00e10 | 7.35e10 | 61 | 8483 | 0.730 | 0.109 | 4.411 |
| 155 | 3.37e13 | 3.11e10 | 7.77e10 | 47 | 8045 | 0.766 | 0.098 | 5.030 |
| 160 | 3.00e13 | 2.43e10 | 3.38e10 | 83 | 8300 | 0.751 | 0.107 | 5.105 |
| 165 | 2.67e13 | 1.93e10 | 2.76e10 | 75 | 8912 | 0.757 | 0.106 | 5.625 |
| 170 | 2.38e13 | 1.52e10 | 2.59e10 | 45 | 8956 | 0.761 | 0.112 | 6.662 |
| 175 | 2.12e13 | 1.22e10 | 2.10e10 | 73 | 7660 | 0.781 | 0.095 | 6.905 |
| 180 | 1.89e13 | 9.53e9 | 1.24e10 | 80 | 8581 | 0.733 | 0.105 | 7.965 |
| 185 | 1.69e13 | 7.62e9 | 1.06e10 | 83 | 8228 | 0.788 | 0.097 | 7.652 |
| 190 | 1.51e13 | 6.08e9 | 1.04e10 | 72 | 7901 | 0.752 | 0.090 | 8.320 |
| 195 | 1.34e13 | 4.76e9 | 6.25e9 | 89 | 8235 | 0.796 | 0.093 | 8.877 |
| 200 | 1.20e13 | 3.83e9 | 5.87e9 | 75 | 8576 | 0.763 | 0.083 | 9.986 |

## How the bent-crystal model works

This section describes `BentCrys2D.py` (copied into OASYS2-XRT as
`orangecontrib/xrt/util/bent_crys_2d.py`) together with the xrt 2.0 code it relies on,
`xrt/backends/raycing/oes/reflect.py`.

### 1) Output direction: dispersive, not specular on the Bragg planes

`BentCrys2D.local_n` returns two normals for each hit point: the local Bragg-plane normal
and the local surface normal. Both are tilted by the bending: the surface slope is
x/(Rs·Fcurv) and y/Rm, and the planes are rotated by the asymmetry angle α. When xrt
receives two normals for a crystal, it computes the outgoing direction like a grating
(`_get_asymmetric_reflection_grating` → `_grating_deflection`):

- **Grating equation:** the part of the wavevector along the surface changes by the part
  of the reciprocal vector H (length 1/d) that lies along the surface. The component
  normal to the surface is then fixed by |k_out| = |k_in|.
- **Each ray is diffracted individually.** A ray slightly off the Bragg angle comes out in
  the direction the grating equation gives for its own energy and angle, so energy is
  spread in angle the way a real crystal does it. In the symmetric Bragg case H is normal
  to the surface, and this reduces to specular reflection.
- **No refraction correction:** xrt sets `wH = 0`, so the direction leaves out the
  refraction shift. The scripts compensate for it in the alignment by using the corrected
  Bragg angle.
- **The direction does not depend on the reflectivity.** Every ray goes to its geometric
  direction and is then weighted by the reflectivity. There is no angular broadening
  inside the reflectivity curve.
- **Where the focusing comes from:** the curvature of the planes and of the surface
  (Rm, Rs·Fcurv). Fcurv changes only the normals, i.e. the geometry. The reflectivity uses
  the crystal's own Rs, without Fcurv.
- **Mismatch in the code:** `local_z`, which defines where rays hit, uses `Rs` without
  Fcurv, while `local_n` uses `Rs·Fcurv`. When Fcurv ≠ 1 the surface shape and the surface
  normals are therefore not quite consistent. The effect is second-order.

The specular alternative, mirroring on the Bragg planes, is only used with mosaicity or
with `volumetricDiffraction=True`. The scripts use neither.

### 2) Penetration: no, rays reflect at the front surface

The crystals are created with `volumetricDiffraction=False`. Each ray is diffracted at the
point where it meets the curved surface, and leaves from that same point:

- **Laue case:** the real beam leaves from the back face, spread over the Borrmann fan.
  For symmetric Laue the fan width is 2t·tanθB, about 0.11 mm for t = 1 mm at 35 keV, and
  less at higher energy. Neither the shift nor the spread is in the ray tracing, and that
  is comparable to the 25–60 µm foci in the SSU runs.
- **Where the thickness t does enter:** only the reflectivity uses it (absorption, the PP
  integration range, the number of ML lamellae).
- **The xrt option is not a fix:** turning on `volumetricDiffraction` would place the
  diffraction point at a random depth. It would also switch the direction to specular on
  the planes, and the code says `local_n_depth` "doesn't work for bent crystal using PP
  or ML".

### 3) Reflectivity: per-ray dynamical theory with bending strain

`bentSi_Crys.get_amplitude` computes a reflectivity for each ray from its own energy and
directions. The steps are:

1. **Local deviation from the Bragg condition** (`deviation_common`). It uses
   α_Z = (λ/d)² + 2(λ/d)(k̂_in·Ĥ), with the asymmetry factor b = γ₀/γ_H taken from the
   actual incoming and outgoing directions, plus χ₀ and χ_H. Because the plane normal is
   local, every point of the bent crystal has its own deviation.
2. **Strain gradient from the bending.**
   - The bending moments M1 and M2 come from Rm and Rs with the anisotropic compliance
     tensor S (`stiffness_tensor` → `ST`).
   - For PP they are projected onto the incoming direction, the outgoing direction and H
     to give the strain-gradient parameter β (`G`, `beta_s`, `beta_p`).
   - For ML they give the gradient of the deviation through the depth (`dalpha_over_dt`).
3. **One of two models:**
   - **PP (Penning–Polder), Laue only.** The Bragg-condition point is placed at the crystal
     mid-plane. `PP_center_fraction` is stored but ignored; the mid-plane is hard-coded
     in `build_deviation_reference`, for PP and ML alike. The deviation is integrated between the
     entrance and exit faces, η_in → η_out, which differ by −βt. The intensity is the sum
     of the two wavefield branches, each with its own absorption,
     exp(−μt(1+B+D)/|γ₀|). It is multiplied by (1 − e^{−2πβc/|β|}), the loss to newly
     created wavefields when the strain gradient is large.
   - **ML (multilamellar).** The thickness is divided into N lamellae, each a perfect thin
     crystal. The deviation step between lamellae follows from the strain gradient: Δy = π/2
     (Laue) or 2 (Bragg), as in `crystal3.f90` (`ML_refinement=1`; the old runs used 4). Each lamella's reflectivity comes from the standard thin-crystal
     formulas (separate Bragg and Laue versions). The lamellae are summed incoherently,
     weighted by the transmission through the lamellae before it and by the absorption
     along the path. The susceptibilities from xrt are complex-conjugated first, to match
     the sign convention of the formulas.
4. **Polarization:** σ and π are computed separately, with a polarization factor of 1 and
   |cos 2θB|. The π factor uses the median over the rays.
5. **Return value:** √R for σ and π, with no phase. That is fine for incoherent ray
   tracing, but there is no phase information.

`useTT=False`, so xrt's own Takagi–Taupin solver (pyTTE) is never used.

## Laue reflectivity check: PP vs ML

> **This section describes the code and runs before the ML fixes.** With the fixes
> (`BUG_ML.md`), the internal ML gives 30.42 µrad for the single crystal below (XOPPY ML
> 34.02, internal PP 30.05), and the internal PP/ML ratio is 0.988 (XOPPY 0.969). The
> comparison of the full beamline below used the old geometry and the old ML code. The
> rerun `run_XRT_SSU_BLM_PP_fixed/` fails like the ML one (see
> [Bragg+Laue beamlines](#bragglaue-beamlines-the-reflection-windows-do-not-overlap)).

The internal reflectivity code (`bentSi_Crys.get_amplitude`) was compared with XOPPY
(`diff_pat`, from xoppylib) for the SSU_BLM Laue crystal. The same rocking-curve
calculation was also made in the full SSU_BLM beamline.

### Single crystal: internal code vs XOPPY

Script: `XOPPY/comparison_SSU_BLM.py`. It writes its output to `XOPPY/SSU_BLM_radii/`:
`diff_pat_PP.dat`, `diff_pat_ML.dat` and `comparison_internal.dat`.

- **Crystal:**
  - Si(311) Laue, asymmetry 20.98°, t = 1 mm.
  - E = 80 keV, θB = 47.34 mrad.
  - Compliance tensor `/home/srio/Oasys2/S_tensor.txt`, identical to the SSU_BLM one.
- **Radii, as SSU_BLM uses them at 80 keV:**
  - XOPPY: RSAG = 121.76 cm, RMER = −7603.94 cm.
  - Internal: Rs = 1217.6 mm, Rm = −76039 mm.
- **Matching the two codes:**
  - XOPPY RMER corresponds to the internal Rm, and RSAG to the internal Rs. This was
    checked through M1/I and M2/I.
  - The internal code builds the in and out directions as xrt does. Its out-angle matches
    XOPPY's to 0.02 µrad.

| s-polarization | peak | integrated (µrad) | half-max range (µrad) |
|---|---|---|---|
| XOPPY PP (MOSAIC=3) | 0.760 | 32.96 | [0.3, 43.5] |
| XOPPY ML (MOSAIC=2) | 0.753 | 34.02 | [−0.5, 43.3] |
| internal PP | 0.699 | 30.05 | [−21.5, 21.8] |
| internal ML | **1.000 (capped)** | **48.77** | [−22.5, 23.4] |

| Integrated reflectivity ratio PP/ML | s-pol | p-pol |
|---|---|---|
| XOPPY | 0.969 | 0.973 |
| internal | 0.616 | 0.610 |

**Findings:**

- **XOPPY's PP and ML agree within 3%.** In the internal code they differ by a factor 1.6.
- **Internal PP:** 9% below XOPPY's PP.
  - The difference is in the absorption: with μ = 0, the internal PP gives 1 − e^{−2πβc/|β|}.
  - The internal absorption exponent is about 3× XOPPY's.
- **Internal ML: 43% above XOPPY's ML. This is the suspect result.**
  - The summed lamellae go above 1 and are capped at 1 over about 40 µrad.
  - XOPPY's ML peaks at 0.75 with the same lamella structure (47 lamellae of 21.2 µm,
    c = 2.06, β = 740 cm⁻¹).
  - Things to check: the per-lamella reflectivity, the transmission and absorption
    weighting between lamellae, and the ×4 refinement of N in `get_amplitude`.
- **The internal profiles are centred 22 µrad lower than XOPPY's.**
  - The internal code puts the Bragg-condition point at the crystal mid-plane, while XOPPY
    puts it at the entrance surface. The shift is half the profile width.
  - The `PP_center_fraction` argument is stored but not used:
    `build_deviation_reference` always uses the mid-plane.
- **Low-reflectivity check:** a run of `XOPPY/comparison.py` with Rm = 1000 mm and
  Rs = 1×10⁶ mm (peak R ≈ 0.02) gave integrated reflectivities of 62.80 µrad for XOPPY PP,
  58.65 for internal PP and 60.12 for internal ML. The internal ML only goes wrong when
  the reflectivity is high and the lamella sum exceeds 1.
- **XOPPY ML needs the radii signs flipped here.** With the SSU_BLM radii the strain
  gradient is negative, and `diff_pat` then gives a negative number of lamellae (−47) and a
  profile that is zero everywhere.
  - The script flips the sign of both radii for the XOPPY ML run, the equivalent of
    `ML_sign=-1` in the internal code.
  - This mirrors the profile in angle but does not change the integrated reflectivity.

### Full beamline: SSU_BLM with PP

Script: `XRT_SSU_BLM_PP.py`, a copy of `XRT_SSU_BLM.py`. Results are in
`run_XRT_SSU_BLM_PP/`.

- **What changed:** the Laue crystal uses `method='PP'`, and the output names say `_PP_`.
- **What did not change:**
  - The Bragg crystal stays on ML, since PP is Laue-only.
  - The Fcurv optimization still traces with ML, so the radii (`RsL`, `RmL`) are identical
    in both runs.

| E [keV] | flux ML | flux PP | PP/ML | flux density PP | x PP [µm] | z PP [µm] | dE PP [keV] | dE ML [keV] |
|---|---|---|---|---|---|---|---|---|
| 35 | 1.21e12 | 5.29e11 | 0.44 | 7.04e12 | 30 | 2000 | 0.009 | 0.010 |
| 40 | 1.59e12 | 7.92e11 | 0.50 | 1.34e13 | 30 | 1800 | 0.012 | 0.013 |
| 45 | 1.98e12 | 1.13e12 | 0.57 | 1.68e13 | 28 | 2200 | 0.015 | 0.016 |
| 50 | 2.37e12 | 1.39e12 | 0.59 | 1.84e13 | 32 | 2100 | 0.020 | 0.022 |
| 55 | 2.66e12 | 1.67e12 | 0.63 | 2.19e13 | 30 | 2300 | 0.030 | 0.031 |
| 60 | 2.88e12 | 1.88e12 | 0.65 | 2.33e13 | 30 | 2400 | 0.042 | 0.047 |
| 65 | 3.11e12 | 2.04e12 | 0.66 | 2.09e13 | 32 | 2700 | 0.063 | 0.067 |
| 70 | 3.25e12 | 2.16e12 | 0.66 | 2.06e13 | 34 | 2900 | 0.083 | 0.090 |
| 75 | 3.38e12 | 2.14e12 | 0.63 | 1.78e13 | 32 | 3400 | 0.114 | 0.120 |
| 80 | 3.38e12 | 2.08e12 | 0.61 | 1.51e13 | 34 | 3800 | 0.160 | 0.163 |
| 85 | 3.38e12 | 1.98e12 | 0.59 | 1.36e13 | 34 | 4000 | 0.193 | 0.190 |
| 90 | 3.11e12 | 1.85e12 | 0.60 | 1.09e13 | 36 | 4300 | 0.246 | 0.254 |
| 95 | 2.67e12 | 1.76e12 | 0.66 | 8.41e12 | 36 | 4900 | 0.314 | 0.329 |
| 100 | 2.30e12 | 1.60e12 | 0.70 | 7.27e12 | 42 | 5000 | 0.390 | 0.387 |
| 105 | 2.00e12 | 1.43e12 | 0.72 | 7.01e12 | 38 | 5000 | 0.442 | 0.423 |
| 110 | 1.71e12 | 1.28e12 | 0.75 | 5.90e12 | 40 | 5100 | 0.537 | 0.558 |
| 115 | 1.47e12 | 1.15e12 | 0.78 | 4.63e12 | 42 | 5500 | 0.638 | 0.646 |
| 120 | 1.26e12 | 1.01e12 | 0.80 | 3.48e12 | 48 | 5800 | 0.735 | 0.780 |
| 125 | 1.08e12 | 8.96e11 | 0.83 | 3.32e12 | 46 | 5600 | 0.831 | 0.858 |
| 130 | 9.34e11 | 7.81e11 | 0.84 | 2.93e12 | 44 | 5800 | 1.009 | 0.999 |
| 135 | 7.92e11 | 6.73e11 | 0.85 | 2.08e12 | 52 | 5900 | 1.162 | 1.141 |
| 140 | 6.66e11 | 5.73e11 | 0.86 | 1.73e12 | 48 | 6600 | 1.413 | 1.308 |
| 145 | 5.49e11 | 4.86e11 | 0.89 | 1.56e12 | 48 | 6300 | 1.421 | 1.375 |
| 150 | 4.55e11 | 4.05e11 | 0.89 | 1.02e12 | 58 | 6800 | 1.599 | 1.664 |
| 155 | 3.75e11 | 3.42e11 | 0.91 | 8.85e11 | 52 | 7000 | 1.787 | 1.809 |
| 160 | 3.14e11 | 2.85e11 | 0.91 | 9.01e11 | 48 | 6200 | 2.048 | 2.000 |
| 165 | 2.58e11 | 2.35e11 | 0.91 | 7.40e11 | 46 | 6400 | 2.090 | 2.145 |
| 170 | 2.14e11 | 1.97e11 | 0.92 | 5.54e11 | 50 | 6600 | 2.125 | 2.159 |
| 175 | 1.75e11 | 1.62e11 | 0.92 | 4.17e11 | 60 | 6100 | 2.275 | 2.292 |
| 180 | 1.43e11 | 1.34e11 | 0.93 | 2.99e11 | 62 | 7000 | 2.412 | 2.376 |
| 185 | 1.18e11 | 1.12e11 | 0.95 | 2.68e11 | 62 | 6400 | 2.461 | 2.461 |
| 190 | 9.81e10 | 9.30e10 | 0.95 | 2.37e11 | 56 | 6600 | 2.584 | 2.546 |
| 195 | 8.02e10 | 7.69e10 | 0.96 | 1.95e11 | 58 | 6500 | 2.730 | 2.691 |
| 200 | 6.59e10 | 6.36e10 | 0.96 | 1.58e11 | 56 | 6900 | 2.780 | 2.860 |

- **Flux:**
  - PP gives less flux than ML at every energy: 0.44× at 35 keV, 0.61× at 80 keV and
    0.96× at 200 keV. Over the whole scan PP gives 0.67× the ML flux.
  - The gap is largest at low energy, where the reflectivity is high and the internal ML
    is capped at 1.
  - The PP flux peaks at 70 keV (2.16×10¹²) and stays nearly flat up to 80 keV.
- **The single-crystal and beamline results agree.** At 80 keV the flux ratio of the
  beamline (0.615) equals the single-crystal integrated-reflectivity ratio (0.616). The
  flux difference therefore comes entirely from the Laue reflectivity model.
- **Corrected SSU_BLM flux:**
  - XOPPY's ML and PP agree, which makes the internal ML the wrong one.
  - Taking XOPPY as the reference, the SSU_BLM flux at 80 keV is about 2.3×10¹² ph/s:
    the PP value raised by the 9% that PP is below XOPPY. The 3.38×10¹² from the ML run
    is too high.
- **Bandwidth, horizontal size and divergences:** unchanged within the ray statistics.
  dE is a few percent smaller with PP below 100 keV.
- **Vertical size:** typically about 10% smaller with PP, since PP puts less weight on the
  outer part of the rocking curve.

## TODO (decisions pending)

1. **Depth reference of the bent-crystal model.** Implement `PP_center_fraction` in
   `BentCrys2D.py` and in OASYS2-XRT `orangecontrib/xrt/util/bent_crys_2d.py`, as in
   `test_fc0_rock/PP_center_fraction.diff`. Decide whether to use 0 (entrance surface, as
   XOPPY) in the scripts, as the widget default and in the `.ows` files.
2. **Alignment of the second crystal.** Decide whether to add a per-energy rocking scan of
   the Laue crystal to `XRT_SSU_BLM.py`, `XRT_SSU_BLM_PP.py` and `XRT_APS_BLM.py`, as in
   `test_fc0_rock/run_fc0_rock.py`. Without it the Bragg+Laue beamlines do not transmit.
3. **Absorption coefficient.** Keep xrt's total μ (current, 0.521 cm⁻¹ for Si at 80 keV) or
   use the photoelectric μ as XOPPY (0.118 cm⁻¹). This changes the Bragg integrated
   reflectivity by about ×4 (17 vs 71 µrad at 80 keV); the effect on the beamline flux has
   not been traced.
4. **BNL_SBDL at 30–40 keV.** Understand the drop to 0.25× the old flux.
5. **After 1–3:** rerun SSU_BLM, SSU_BLM_PP and APS_BLM (and SSU_SBDL and BNL_SBDL if the
   depth reference changes), update this file, `RERUN_RESULTS.md` and the `INFO.md` of the
   OASYS workspaces. The Bragg+Laue workspaces (`oasys_SSU_BLM`, `oasys_APS_BLM`) need the
   same changes as the scripts.
