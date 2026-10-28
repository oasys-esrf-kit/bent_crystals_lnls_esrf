
#
# script to calculate crystal diffraction profiles (created by XOPPY:crystal)
# SSU_BLM Laue crystal radii at 80 keV (RsL = 1217.64 mm, RmL = 76039.4 mm; xrt uses Rs=RsL, Rm=-RmL)
# run in SSU_BLM_radii/ so the files of comparison.py are not overwritten
#

import numpy
import os
os.chdir(os.path.join(os.path.dirname(os.path.abspath(__file__)), "SSU_BLM_radii"))
from xoppylib.crystals.tools import bragg_calc2, run_diff_pat
try: import xraylib
except: print("xraylib not available")
from dabax.dabax_xraylib import DabaxXraylib

#
# run bragg_calc (preprocessor) and create file xcrystal.bra
#
bragg_dictionary = bragg_calc2(
    descriptor = "Si",
    hh         = 3,
    kk         = 1,
    ll         = 1,
    temper     = 1.0,
    emin       = 79900.0,
    emax       = 80100.0,
    estep      = 0.4,
    ANISO_SEL  = 0,
    fileout    = "xcrystal.bra",
    do_not_prototype = 0,  # 0=use site groups (recommended), 1=use all individual sites
    verbose = False,
    material_constants_library = DabaxXraylib(file_f1f2="f1f2_Windt.dat"),
    )

#
# run external (fortran) diff_pat (note that some parameters may not be used)
# twice: MOSAIC=3 (bent Penning-Polder) and MOSAIC=2 (bent multilamellar); diff_pat.dat is copied after each run
# diff_pat ML fails (negative number of lamellae, zero profile) when the strain gradient is negative, as here;
# for ML the radii signs are flipped (same as ML_sign=-1 in the internal code): the profile is mirrored in angle,
# the integrated reflectivity is unchanged
#
import shutil
for mosaic, label, rsign in ((3, "PP", 1.0), (2, "ML", -1.0)):
    run_diff_pat(
        bragg_dictionary,
        preprocessor_file  = "xcrystal.bra",
        descriptor         = "Si",
        MOSAIC             = mosaic,
        GEOMETRY           = 1,
        SCAN               = 2,
        UNIT               = 1,
        SCANFROM           = -100.0,
        SCANTO             = 100.0,
        SCANPOINTS         = 1000,
        ENERGY             = 80000.0,
        ASYMMETRY_ANGLE    = 20.98,
        THICKNESS          = 0.1,
        MOSAIC_FWHM        = 0.2,
        RSAG               = rsign * 121.7636666,
        RMER               = rsign * -7603.9407473,
        ANISOTROPY         = 3,
        POISSON            = 0.2201,
        CUT                = "2 -1 -1 ; 1 1 1 ; 0 0 0",
        FILECOMPLIANCE     = "/home/srio/Oasys2/S_tensor.txt",
        )
    shutil.copy("diff_pat.dat", "diff_pat_%s.dat" % label)

data_pp = numpy.loadtxt("diff_pat_PP.dat", skiprows=5)
data_ml = numpy.loadtxt("diff_pat_ML.dat", skiprows=5)

#
# example plot
#
if True:
    from srxraylib.plot.gol import plot
    plot(data_pp[:,0], data_pp[:,-1], data_pp[:,0], data_pp[:,-2], data_ml[:,0], data_ml[:,-1], data_ml[:,0], data_ml[:,-2],
         ytitle='Crystal reflectivity', legend=['PP s-polarization','PP p-polarization','ML s-polarization','ML p-polarization'])

#
# same crystal with the internal reflectivity code (bentSi_Crys.get_amplitude in ../BentCrys2D.py)
#
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import BentCrys2D as BC2D
from xrt.backends.raycing.physconsts import CH

# same inputs as in run_diff_pat above (XOPPY lengths in cm, here in mm)
energy = 80000.0                                      # ENERGY [eV]
asymmetry = numpy.radians(20.98)                      # ASYMMETRY_ANGLE: angle between H and surface normal
thickness = 0.1 * 10                                  # THICKNESS
Rm = -7603.9407473 * 10                                # RMER
Rs = 121.7636666 * 10                                   # RSAG
ST = numpy.loadtxt("/home/srio/Oasys2/S_tensor.txt")  # FILECOMPLIANCE
dtheta = numpy.arange(-100.0, 100.01, 0.05)           # Th-ThB{in} [microrad], wider than XOPPY scan

def internal_reflectivity(method, ML_sign=1.0):
    crystal = BC2D.bentSi_Crys(hkl=[3, 1, 1], geom='Laue reflected', t=thickness, aLrad=asymmetry,
                               volumetricDiffraction=False, useTT=False, Rs=Rs, Rm=Rm, method=method,
                               ML_sign=ML_sign, PP_center_fraction=0.5, ST=ST)
    # local crystal frame as in xrt: surface normal along z, diffraction plane y-z
    n = numpy.array([0.0, 0.0, 1.0])
    H = numpy.array([0.0, numpy.sin(asymmetry), numpy.cos(asymmetry)])
    wavelength = CH / energy
    theta = numpy.arcsin(wavelength / 2 / crystal.d) + dtheta * 1e-6  # glancing angle to the planes
    phi = numpy.arctan2(H[2], H[1]) - (numpy.pi / 2 + theta)
    v_in = numpy.array([0 * phi, numpy.cos(phi), numpy.sin(phi)])     # v_in . H = -sin(theta)
    # outgoing direction as in xrt: momentum along the surface conserved, |k_out| = |k_in|
    H_t = H - n * (n @ H)
    v_out = v_in - numpy.outer(n, n @ v_in) + (wavelength / crystal.d) * H_t[:, None]
    v_out[2] = -numpy.sqrt(1 - v_out[0]**2 - v_out[1]**2)            # Laue: stays inside the crystal
    one = numpy.ones_like(theta)
    oe_normal = [H[0] * one, H[1] * one, H[2] * one, n[0] * one, n[1] * one, n[2] * one]
    amp_s, amp_p = crystal.get_amplitude(energy * one, v_in[2], v_out[2], v_in.T @ H,
                                         xd=0 * one, yd=0 * one, ad=v_in[0], bd=v_in[1], cd=v_in[2],
                                         oeNormal=oe_normal, outBN=list(v_out))
    return numpy.abs(amp_s)**2, numpy.abs(amp_p)**2

pp_s, pp_p = internal_reflectivity('PP')
ml_s, ml_p = internal_reflectivity('ML', ML_sign=-1.0)  # ML_sign as for the Laue crystal in XRT_SSU_BLM.py

def profile_summary(name, x, y):
    integral = numpy.trapezoid(y, x)
    half = x[y >= y.max() / 2]
    print("%-14s peak %.5f  integral %7.2f urad  centroid %8.1f urad  half-max [%7.1f, %7.1f] urad" %
          (name, y.max(), integral, numpy.trapezoid(y * x, x) / integral, half.min(), half.max()))
    return integral

print("\nComparison XOPPY (diff_pat) vs internal code (bentSi_Crys):")
integrals = {}
integrals["XOPPY PP s"] = profile_summary("XOPPY PP s", data_pp[:, 0], data_pp[:, -1])
integrals["XOPPY PP p"] = profile_summary("XOPPY PP p", data_pp[:, 0], data_pp[:, -2])
integrals["XOPPY ML s"] = profile_summary("XOPPY ML s", data_ml[:, 0], data_ml[:, -1])
integrals["XOPPY ML p"] = profile_summary("XOPPY ML p", data_ml[:, 0], data_ml[:, -2])
integrals["internal PP s"] = profile_summary("internal PP s", dtheta, pp_s)
integrals["internal PP p"] = profile_summary("internal PP p", dtheta, pp_p)
integrals["internal ML s"] = profile_summary("internal ML s", dtheta, ml_s)
integrals["internal ML p"] = profile_summary("internal ML p", dtheta, ml_p)

print("\nIntegrated reflectivity ratio PP/ML:")
for code in ("XOPPY", "internal"):
    for pol in ("s", "p"):
        print("%-9s %s-pol  PP/ML = %.4f" % (code, pol, integrals["%s PP %s" % (code, pol)] / integrals["%s ML %s" % (code, pol)]))

numpy.savetxt("comparison_internal.dat", numpy.column_stack((dtheta, pp_s, pp_p, ml_s, ml_p)),
              header="Th-ThB{in}[microrad] PP_s PP_p ML_s ML_p")

if True:
    from srxraylib.plot.gol import plot
    plot(data_pp[:, 0], data_pp[:, -1], data_ml[:, 0], data_ml[:, -1], dtheta, pp_s, dtheta, ml_s,
         xtitle='Th-ThB{in} [microrad]', ytitle='Crystal reflectivity (s-polarization)',
         legend=['XOPPY diff_pat (PP)', 'XOPPY diff_pat (ML)', 'internal PP', 'internal ML'])

#
# end script
#


