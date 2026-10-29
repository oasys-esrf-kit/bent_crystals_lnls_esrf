#
# Integrated reflectivities with XOPPY (diff_pat) for the SSU crystals, with the radii used
# in the runs at each energy:
#   SSU_BLM Bragg: Si(311) symmetric, t = 10 mm, ML (MOSAIC=2), Rm = RmB, flat sagittally
#   SSU_BLM Laue:  Si(311), asymmetry 20.98 deg, t = 1 mm, PP (MOSAIC=3), Rs = RsL, Rm = -RmL
#   SSU_SBDL Laue: Si(111), asymmetry 35.26 deg, t = 1 mm, PP (MOSAIC=3), Rs = -Rs, Rm = Rm
# radii from ../run_XRT_SSU_BLM_test_fc0_rock and ../run_XRT_SSU_SBDL_fixed (res_*_energy.txt);
# compliance tensors from the S_tensor.txt the scripts write in their run folders.
# s-polarization. Results in SSU_integrated/integrated_SSU.dat
#
import glob, os, shutil
import numpy
from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
from dabax.dabax_xraylib import DabaxXraylib

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "SSU_integrated")
os.makedirs(OUT, exist_ok=True)
os.chdir(OUT)
res_blm = numpy.loadtxt(glob.glob(os.path.join(HERE, "..", "run_XRT_SSU_BLM_test_fc0_rock", "res_*_energy.txt"))[0], skiprows=1)
res_sbdl = numpy.loadtxt(glob.glob(os.path.join(HERE, "..", "run_XRT_SSU_SBDL_fixed", "res_*_energy.txt"))[0], skiprows=1)
ST_BLM = os.path.join(HERE, "..", "run_XRT_SSU_BLM_test_fc0_rock", "S_tensor.txt")
ST_SBDL = os.path.join(HERE, "..", "run_XRT_SSU_SBDL_fixed", "S_tensor.txt")
energies = [35e3, 50e3, 65e3, 80e3, 100e3, 120e3, 140e3, 160e3, 180e3, 200e3]

def integrated(hkl, energy, mosaic, geometry, asymmetry, thickness_mm, rsag_mm, rmer_mm, compliance, scan=400.0, points=4001):
    bd = bragg_calc2(descriptor="Si", hh=hkl[0], kk=hkl[1], ll=hkl[2], temper=1.0,
                     emin=energy - 100, emax=energy + 100, estep=0.4, ANISO_SEL=0, fileout="xcrystal.bra",
                     do_not_prototype=0, verbose=False,
                     material_constants_library=DabaxXraylib(file_f1f2="f1f2_Windt.dat"))
    best = None
    for sign in (1.0, -1.0):  # diff_pat gives a zero profile for one sign of the strain gradient (ML)
        run_diff_pat(bd, preprocessor_file="xcrystal.bra", descriptor="Si", MOSAIC=mosaic, GEOMETRY=geometry,
                     SCAN=2, UNIT=1, SCANFROM=-scan, SCANTO=scan, SCANPOINTS=points, ENERGY=energy,
                     ASYMMETRY_ANGLE=asymmetry, THICKNESS=thickness_mm / 10, MOSAIC_FWHM=0.2,
                     RSAG=sign * rsag_mm / 10, RMER=sign * rmer_mm / 10, ANISOTROPY=3, POISSON=0.2201,
                     CUT="2 -1 -1 ; 1 1 1 ; 0 0 0", FILECOMPLIANCE=compliance)
        d = numpy.loadtxt("diff_pat.dat", skiprows=5)
        x, r = d[:, 0], d[:, -1]
        out = (numpy.trapezoid(r, x), r.max(), max(r[0], r[-1]), sign)
        if best is None or out[0] > best[0]:
            best = out
        if mosaic == 3:  # PP: the sign only mirrors the profile
            break
    return best

rows = []
for E in energies:
    i = numpy.argmin(abs(res_blm[:, 0] - E)); j = numpy.argmin(abs(res_sbdl[:, 0] - E))
    RmB, RsL, RmL = res_blm[i, 9:12]
    Rs, Rm = res_sbdl[j, 9:11]
    B = integrated([3, 1, 1], E, 2, 0, 0.0, 10.0, 1e32, RmB, ST_BLM, scan=2000.0, points=16001)
    L3 = integrated([3, 1, 1], E, 3, 1, 20.98, 1.0, RsL, -RmL, ST_BLM)
    L1 = integrated([1, 1, 1], E, 3, 1, 35.26, 1.0, -Rs, Rm, ST_SBDL)
    for name, v in (("Bragg311 ML", B), ("Laue311 PP", L3), ("Laue111 PP", L1)):
        if v[2] > 0.01 * v[1]:
            print("WARNING %s at %.0f keV: profile not zero at the scan edge (%.3g)" % (name, E / 1e3, v[2]))
    rows.append((E, B[0], L3[0], L1[0], B[1], L3[1], L1[1]))
    print("%5.0f keV  Bragg311 %8.2f urad (peak %.3f)  Laue311 %7.2f (peak %.3f)  Laue111 %7.2f (peak %.3f)" %
          (E / 1e3, B[0], B[1], L3[0], L3[1], L1[0], L1[1]), flush=True)
numpy.savetxt("integrated_SSU.dat", numpy.array(rows),
              header="E_eV Int_Bragg311_ML Int_Laue311_PP Int_Laue111_PP [urad] peak_Bragg311 peak_Laue311 peak_Laue111")
