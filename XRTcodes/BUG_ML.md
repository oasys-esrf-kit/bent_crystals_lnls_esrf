# Bug: internal ML reflectivity (fixed)

The multilamellar (ML) reflectivity in `bentSi_Crys.get_amplitude` had three bugs. In the
Laue case it was too high: its integrated reflectivity was 1.43× XOPPY's for the SSU_BLM
crystal at 80 keV. In the Bragg case absorption had no effect.

The three bugs are **fixed** in both copies of the code:

- `XRTcodes/BentCrys2D.py`
- `OASYS2-XRT/orangecontrib/xrt/util/bent_crys_2d.py`

The scans have been rerun with the fixes (`run_XRT_*_fixed/`; see `RERUN_RESULTS.md` and
`RUN_XRT.md`). The Bragg+Laue beamlines now need a decision on the depth reference and the
alignment of the second crystal.

Reference: `crystal3.f90`, `CRYSTAL_BENT_ML`, in https://github.com/srio/CRYSTAL (used by
XOPPY `diff_pat`, MOSAIC=2).

## 1. Missing transmission update in the Laue branch

- In `crystal3.f90`, `tsmall_old = tsmall` runs after both the Bragg and the Laue
  branches. This makes `prodt` the product of the transmissions of the lamellae in front.
- The Python Laue branch did not have this line, so `prodt` stayed at 1 and every lamella
  reflected the full incident beam.
- The sum then exceeded 1, and `np.clip` capped it.
- **Fix:** `tsmall_old = tsmall` after `Rtot += rsmall * prodt * att` in the Laue branch.

## 2. `ML_refinement=4`

- Fortran uses the lamella step Δy = π/2 (Laue) or 2 (Bragg). Python used 4× more
  lamellae.
- ML sums the lamellae incoherently, so a finer subdivision changes the result instead of
  converging. The step is part of the model.
- **Fix:** `ML_refinement=1` in both calls to `build_deviation_reference`.

## 3. Sign convention of the susceptibilities

- xrt returns Ψ₀, Ψ_H and Ψ_H̄ as the complex conjugates of the values the ML formulas
  expect, which `crystal3.f90` and XOPPY use. At 80 keV for Si(311):
  - xrt: Ψ₀ = −1.507e−7 **+** 1.29e−10 i
  - XOPPY: Ψ₀ = −1.506e−7 **−** 2.92e−11 i
- **Effect:**
  - Absorption turned into gain in the lamella transmission. `tsmall` came out above 1 (up
    to 1.05) for off-Bragg lamellae.
  - In a thick Bragg crystal the product over about 200 lamellae reached about 10⁴. That
    cancelled the exit attenuation, so the reflectivity stayed at about 1 over the whole
    bending range, independently of absorption.
- **Fix:** in `ml_pol`, use `np.conj(PsiH)`, `np.conj(PsiH_)` and `np.conj(Psi0).imag`
  for `qsmall` and `g`.

## Remaining difference: the absorption coefficient (a modelling choice)

- **The two values:**
  - The internal code takes μ from xrt (`get_absorption_coefficient`, from Im n): 0.521 cm⁻¹
    for Si at 80 keV.
  - XOPPY uses 0.118 cm⁻¹ (photoelectric only).
  - xrt's value presumably includes incoherent (Compton) scattering, which is large at
    these energies.
- **The internal code is unchanged here: it uses xrt's μ.**
- **Check:** with all three fixes and μ (and Im Ψ₀) scaled to XOPPY's value (a test only),
  the internal ML reproduces XOPPY. See the "+ XOPPY μ" column below.

## Test: integrated s-reflectivity at 80 keV (µrad)

Si(311), compliance tensor `/home/srio/Oasys2/S_tensor.txt`. Laue crystal: asymmetry
20.98°, t = 1 mm, `ML_sign` −1. Bragg crystal: symmetric, t = 10 mm, Rm = 647668 mm,
flat sagittally (the SSU_BLM Bragg crystal).

| Case | XOPPY ML | Before the fixes | Fixes 1+2 | Fixes 1+2+3 (current code) | Fixes 1+2+3 + XOPPY μ (test) |
|---|---|---|---|---|---|
| Laue, SSU_BLM radii (Rm = −76039 mm, Rs = 1217.6 mm) | 34.02 | 48.77 (capped at 1) | 34.82 | **30.42** | 34.04 |
| Laue, Rm = 1000 mm, Rs = 1e6 mm | 62.85 | 60.12 | 64.37 | **56.35** | 63.17 |
| Bragg, SSU_BLM | 67.45 | 315 (R ≈ 1 over 330 µrad) | 309 | **17.1** | 71.4 |

- **Laue (current code):** about 11% below XOPPY, the same as the internal PP (30.05 and
  58.65 µrad). The internal PP/ML ratio is now 0.988, against XOPPY's 0.969; before the
  fixes it was 0.616.
- **Bragg (current code):**
  - 17 µrad, against 67 µrad for XOPPY.
  - The thick bent Bragg crystal is limited by absorption with depth, so its integrated
    reflectivity scales roughly as 1/μ. This is why the choice of μ matters much more for
    Bragg than for Laue.
- **Profile position:**
  - The internal profiles are centred on the crystal mid-plane, while XOPPY's start at
    the entrance surface.
  - `PP_center_fraction` is not used.
  - This shifts the profiles but not the integrals.
- **XOPPY radii sign:** XOPPY ML needs radii of the sign that makes the strain gradient
  positive. Otherwise it gives a negative number of lamellae and a zero profile.

## Consequences for the runs

- **The flux of every beamline traced with ML is affected:** SSU_BLM, APS_BLM and
  BNL_SBDL. SSU_SBDL uses PP for its Laue crystals and has no Bragg crystal.
- **Laue crystals:** the old ML flux was too high, by about 1.6× at 80 keV for SSU_BLM
  (see "Laue reflectivity check" in `RUN_XRT.md`).
- **Bragg crystals (SSU_BLM, APS_BLM):** the old ML reflectivity ignored absorption. With
  the fix, the Bragg crystal reflects much less. In the rerun, the Bragg and Laue reflection windows no longer
  overlap and the Bragg+Laue flux falls to 10⁻⁴–10⁻³ of the old values (see `RUN_XRT.md`).
- **The OASYS workspaces** (`oasys_*`) use `bent_crys_2d.py` from OASYS2-XRT. They were
  rerun with the fix; their `INFO.md` show the results before and after.

## To do

1. Decide which absorption coefficient to use: xrt total (current) or photoelectric
   (as XOPPY).
2. Done: scripts rerun and `RUN_XRT.md` updated. See the TODO list there for what remains.
