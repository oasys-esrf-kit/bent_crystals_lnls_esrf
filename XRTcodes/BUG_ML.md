# Bug: internal ML reflectivity too high for Laue crystals

The multilamellar (ML) reflectivity in `bentSi_Crys.get_amplitude` is too high for Laue
crystals. The cause is one missing line in the Laue branch of `ml_pol`. A second
difference from the Fortran reference (4× more lamellae) adds a smaller error. With both
changed, the internal ML matches XOPPY's ML within 2.4% in both test cases.

The code is **not corrected yet**. The fixes below were tested on scratch copies only.

Reference: `crystal3.f90`, `CRYSTAL_BENT_ML`, in https://github.com/srio/CRYSTAL (used by
XOPPY `diff_pat`, MOSAIC=2).

Both copies of the code are affected:

- `XRTcodes/BentCrys2D.py`
- `OASYS2-XRT/orangecontrib/xrt/util/bent_crys_2d.py`

## 1. Missing transmission update in the Laue branch (main cause)

- In `crystal3.f90`, `tsmall_old = tsmall` runs after both the Bragg and the Laue
  branches. This makes `prodt` the product of the transmissions of the lamellae in
  front, so each lamella only sees what is left of the incident beam.
- The Python Bragg branch has this line (`BentCrys2D.py:523`), but the Laue branch does
  not. In Laue, `prodt` therefore stays at 1 and every lamella reflects the full beam.
- The incoherent sum then goes above 1, and `np.clip(Rtot, 0.0, 1.0)` hides this by
  capping at 1.
- The fix is one line after `Rtot += rsmall * prodt * att` in the Laue branch:

  ```python
  tsmall_old = tsmall
  ```

## 2. `ML_refinement=4` (smaller effect)

- Fortran uses the step Δy = π/2 between lamellae for Laue (2 for Bragg). For the
  SSU_BLM Laue crystal at 80 keV that gives 47 lamellae of 21 µm.
- Python multiplies the number of lamellae by 4 (`BentCrys2D.py:359-360`,
  `ML_refinement=4`), giving a step of π/8.
- The ML model sums lamellae incoherently, so it is not a discretization that converges
  as the lamellae get thinner. The π/2 step is part of the model's definition, and 4×
  thinner lamellae change the answer.
- Setting `ML_refinement=1` matches Fortran.

## Test: integrated s-reflectivity at 80 keV (µrad)

Si(311) Laue, asymmetry 20.98°, t = 1 mm, compliance tensor
`/home/srio/Oasys2/S_tensor.txt`. The internal code is called as in
`XOPPY/comparison_SSU_BLM.py` with `method='ML'` and `ML_sign=-1`.

| | SSU_BLM radii (Rm = −76039 mm, Rs = 1217.6 mm) | Rm = 1000 mm, Rs = 1e6 mm |
|---|---|---|
| XOPPY ML (reference) | 34.02 | 62.85 |
| internal ML, current code | 48.77 (+43%, capped at 1) | 60.12 (−4%) |
| + `tsmall_old` fix | 37.09 (+9%) | 63.63 (+1%) |
| + `ML_refinement=1` | 44.79 (still capped) | 60.89 |
| both fixes | **34.82 (+2.4%)**, peak 0.86 | **64.37 (+2.4%)** |

- **Why the low-reflectivity case looked fine:** its peak reflectivity is only about
  0.02, so beam depletion barely matters. That case hid the bug. It only shows up when
  the reflectivity is high.
- **Bragg crystals:** the Bragg branch already has the transmission update. But
  `ML_refinement=4` also applies to Bragg crystals (SSU_BLM, APS_BLM), so their
  reflectivity probably changes too. This has not been checked against XOPPY.
- **XOPPY radii sign:** XOPPY ML runs the Rm = 1000 mm case with positive radii, and gives
  −2899 lamellae and a zero profile with negative radii. For the SSU_BLM radii it is the
  other way round. The ML step sign follows the strain gradient: with the wrong sign the
  lamella count goes negative and the profile is zero.

## Consequences for the runs

- All beamlines traced with ML are affected: SSU_BLM, APS_BLM and BNL_SBDL. SSU_SBDL uses
  PP.
- For SSU_BLM the ML flux is about 1.6× too high at 80 keV and 2.3× too high at 35 keV.
  See "Laue reflectivity check: PP vs ML" in `RUN_XRT.md`.
- The Fcurv optimization traces with ML in every script. Its result (the focus position)
  depends on geometry, not on the reflectivity magnitude, so it should be largely
  unaffected.

## To do after the fix

1. Apply both fixes to the two copies of the code.
2. Rerun `XOPPY/comparison.py` and `XOPPY/comparison_SSU_BLM.py`.
3. Check a Bragg ML case against XOPPY.
4. Rerun SSU_BLM, APS_BLM and BNL_SBDL, and update `RUN_XRT.md`.
