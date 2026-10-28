# XRT runs: bent-crystal monochromators

The four scripts `XRT_APS_BLM.py`, `XRT_BNL_SBDL.py`, `XRT_SSU_BLM.py` and
`XRT_SSU_SBDL.py` were run with xrt 2.0.0 (`~/venv/oasys2`), 500 000 rays, full
energy lists. All four finished with exit code 0. These results have not yet been
compared with xrt 1.x results (xrt 1.x is not installed in this environment).

> **Warning: the ML (multilamellar) reflectivity is suspect.** For the SSU_BLM Laue
> crystal at 80 keV, the internal ML gives an integrated reflectivity 1.43× XOPPY's ML,
> while the internal PP is within 9% of XOPPY's PP. The SSU_BLM flux computed with ML is
> therefore about 1.6× too high at 80 keV, and up to 2.3× too high at 35 keV. APS_BLM and
> BNL_SBDL also use ML and may be affected in the same way; they have not been checked.
> See [Laue reflectivity check: PP vs ML](#laue-reflectivity-check-pp-vs-ml).

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
- **Outputs:** each folder contains `res_*_energy.txt`, the png images for each energy,
  `S_tensor.txt` and `run.log`.
- **Header fix:** in `XRT_SSU_BLM.py` and `XRT_SSU_SBDL.py`, the
  `dE` column was missing from the result-file header. It has now been added, and
  the headers of the existing result files were corrected.

## What each script computes

Each script makes an energy scan of a sagittally focusing bent-crystal
monochromator. At every energy:

1. **Crystal setup:**
   - The script computes the Bragg angle, refraction-corrected. These runs used optlnls
     `calc_Bragg_angle(corrected=True)`. The scripts now use xrt instead,
     `get_Bragg_angle(E) - get_dtheta_symmetric_Bragg(E)`, which differs by at most
     1.4 µrad (Si(311) at 35 keV).
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

## Per-energy results

The numbers below are taken from `run_*/res_*_energy.txt`. Sizes are given in µm here (the files use mm).

### APS_BLM (Si311, ML)

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

### BNL_SBDL (Si111 double Laue, ML, + M1)

| E [keV] | flux0 | flux | flux density | x [µm] | z [µm] | x' [mrad] | z' [mrad] | dE [keV] | Laue transmission |
|---|---|---|---|---|---|---|---|---|---|
| 30 | 1.13e15 | 3.73e12 | 9.88e13 | 295 | 128 | 1.960 | 0.073 | 0.021 | 0.093 |
| 40 | 4.99e14 | 5.24e12 | 1.72e14 | 275 | 126 | 1.599 | 0.071 | 0.046 | 0.123 |
| 50 | 2.11e14 | 3.66e12 | 1.22e14 | 255 | 126 | 1.372 | 0.065 | 0.084 | 0.151 |
| 60 | 8.72e13 | 1.65e12 | 5.59e13 | 260 | 121 | 1.187 | 0.059 | 0.127 | 0.148 |
| 70 | 3.55e13 | 4.75e11 | 2.09e13 | 239 | 109 | 1.052 | 0.050 | 0.166 | 0.102 |

### SSU_BLM (Si311, ML)

The flux below is suspect, because the internal ML overestimates the Laue reflectivity. The
same run with PP is in [Laue reflectivity check: PP vs ML](#laue-reflectivity-check-pp-vs-ml).

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

### SSU_SBDL (Si111 double Laue, PP)

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
     mid-plane (`PP_center_fraction=0.5`). The deviation is integrated between the
     entrance and exit faces, η_in → η_out, which differ by −βt. The intensity is the sum
     of the two wavefield branches, each with its own absorption,
     exp(−μt(1+B+D)/|γ₀|). It is multiplied by (1 − e^{−2πβc/|β|}), the loss to newly
     created wavefields when the strain gradient is large.
   - **ML (multilamellar).** The thickness is divided into N lamellae, each a perfect thin
     crystal. The deviation step between lamellae follows from the strain gradient, and N
     is refined ×4. Each lamella's reflectivity comes from the standard thin-crystal
     formulas (separate Bragg and Laue versions). The lamellae are summed incoherently,
     weighted by the transmission through the lamellae before it and by the absorption
     along the path.
4. **Polarization:** σ and π are computed separately, with a polarization factor of 1 and
   |cos 2θB|. The π factor uses the median over the rays.
5. **Return value:** √R for σ and π, with no phase. That is fine for incoherent ray
   tracing, but there is no phase information.

`useTT=False`, so xrt's own Takagi–Taupin solver (pyTTE) is never used.

## Laue reflectivity check: PP vs ML

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
