
import os
import time
import numpy as np
from scipy.signal import savgol_filter

import sys; sys.path.append(os.path.join('..', '..', '..'))  # analysis:ignore
import xrt

from xrt.backends.raycing import BeamLine
from xrt.backends.raycing.sources import BeamProxy
from xrt.backends.raycing.materials import Material

from xrt.plotter import XYCAxis, XYCPlot
from xrt.runner import run_ray_tracing

showIn3D = False

#
# build_beamline_code
#
def build_beamline(name=''):

    bl = BeamLine()
    bl.name = name


    #
    # Component index: 0 (DW)
    #    
    from xrt.backends.raycing.sources import BendingMagnet
    bl.DW = BendingMagnet(
        bl,
        name="DW",
        center=[0,0,0],
        B0=0.6,
        eE=7.0,
        eI=0.1,
        eEpsilonX=6.37,
        eEpsilonZ=1.29,
        eEspread=0.0,
        eSigmaX=91.0,
        eSigmaZ=30.0,
        xPrimeMax=1.0,  # in mrad, to match accepting slit
        zPrimeMax=0.5,  # in mrad, to match accepting slit
        distE="eV",
        eMin=43780.0,
        eMax=44220.0,
        nrays=500000,
        )

    #
    # Component index: 1 (Mask1)
    #    
    from xrt.backends.raycing.apertures import RectangularAperture
    bl.Mask1 = RectangularAperture(
        bl,
        name='Mask1',
        center=[0, 10000, 0],
        kind=('left', 'right', 'bottom', 'top'),
        opening=[-10.0/2, 10.0/2, -5.0/2, 5.0/2])                

    #
    # Component index: 2 (Bragg)
    #    
    from xrt.backends.raycing.oes import JohannCylinder
    from orangecontrib.xrt.util.bent_crys_2d import bentSi_Crys  # also installs the hook passing ray directions to bentSi_Crys
    from orangecontrib.xrt.util.stiffness_tensor import stiffness_tensor
    CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7, 63.9, 79.6], hkl=[5, -1, -1], sag=[0, 1, -1], unitfac=1e3)
    alpha = 0.0
    crystal = bentSi_Crys(
        hkl=[3, 1, 1],
        geom='Bragg reflected',
        t=10.0,
        aLrad=alpha,
        volumetricDiffraction=False,
        useTT=False,
        Rs=1e100,
        Rm=217527.5769,
        method='ML',
        ML_sign=1.0,
        ST=ST)
    thetaB = crystal.get_Bragg_angle(44000.0) - crystal.get_dtheta_symmetric_Bragg(44000.0)
    bl.Bragg = JohannCylinder(
        bl,
        name='Bragg',
        center=[0, 26500, 0],
        pitch=np.pi - alpha - thetaB,
        alpha=alpha,
        material=(crystal,),
        Rm=217527.5769,
        crossSection='circular',
        limPhysX=[-50, 50],
        limPhysY=[-500, 500],
        )

    #
    # Component index: 3 (Laue)
    #    
    from orangecontrib.xrt.util.bent_crys_2d import BentCrys2D, bentSi_Crys
    from orangecontrib.xrt.util.stiffness_tensor import stiffness_tensor
    CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7, 63.9, 79.6], hkl=[5, -1, -1], sag=[0, 1, -1], unitfac=1e3)
    alpha_hkl_nrm = np.pi / 2 - np.arccos(np.dot([3, 1, 1], [5, -1, -1]) / np.linalg.norm([3, 1, 1]) / np.linalg.norm([5, -1, -1]))
    alpha = -alpha_hkl_nrm
    crystal = bentSi_Crys(
        hkl=[3, 1, 1],
        geom='Laue reflected',
        t=0.6,
        aLrad=alpha,
        volumetricDiffraction=False,
        useTT=False,
        Rs=1686.0604,
        Rm=-28546.8644,
        method='ML',
        ML_sign=-1.0,
        PP_center_fraction=0.5,
        ST=ST)
    thetaB = crystal.get_Bragg_angle(44000.0) - crystal.get_dtheta_symmetric_Bragg(44000.0)
    bl.Laue = BentCrys2D(
        bl,
        name='Laue',
        center=[0, 26787.3039270, -50],
        pitch=np.pi/2 + alpha - thetaB,
        alpha=alpha,
        material=(crystal,),
        Rm=-28546.8644,
        Rs=1686.0604,
        Fcurv=1.0,
        limPhysX=[-50, 50],
        limPhysY=[-50, 50],
        )

    #
    # Component index: 4 (SCRexp)
    #    
    from xrt.backends.raycing.screens import Screen
    bl.SCRexp = Screen(
        bl,
        name="SCRexp",
        center=[0, 52000, -50],
        )

    #
    #
    #
    return bl


#
# run_process_code
#
def run_process(bl):

    import numpy as np
    t0 = time.time()

    beams_to_plot = dict()

    #
    # Component index: 0 (DW: BendingMagnet)
    #
    beam = bl.DW.shine()
    if bl.dump_beams_flag: dump_beam(bl, "DW", beam)

    #
    # Component index: 1 (Mask1: RectangularAperture)
    #
    beam_local = bl.Mask1.propagate(beam)
    if bl.dump_beams_flag: dump_beam(bl, "Mask1", beam_local)

    #
    # Component index: 2 (Bragg: JohannCylinder)
    #
    beam, _ = bl.Bragg.reflect(beam)
    if bl.dump_beams_flag: dump_beam(bl, "Bragg", beam)

    #
    # Component index: 3 (Laue: BentCrys2D)
    #
    beam, _ = bl.Laue.reflect(beam)
    if bl.dump_beams_flag: dump_beam(bl, "Laue", beam)

    #
    # Component index: 4 (SCRexp: Screen)
    #
    beam_local = bl.SCRexp.expose(beam)
    beams_to_plot['SCRexp'] = beam_local
    if bl.dump_beams_flag: dump_beam(bl, "SCRexp", beam)

    #
    #
    #
    dt = time.time() - t0
    print("Time needed to create source and trace system %.3f sec" % dt)
    if showIn3D: bl.prepare_flow()
    return beams_to_plot



def make_plot(bl, beamName, limitsH=None, limitsV=None, bins=1024, climits=None, cbins=256, ):
    xaxis = XYCAxis(label='x', limits=limitsH, unit='um', bins=bins, ppb=int(1024/bins))
    yaxis = XYCAxis(label='z', limits=limitsV, unit='um', bins=bins, ppb=int(1024/bins))
    caxis = XYCAxis(label='energy', limits=climits,
                    unit='eV', bins=cbins, fwhmFormatStr="%.2f", ppb=int(512/cbins))

    fname = os.path.join(bl.dump_scores_folder, "%s.png" % beamName)
    plot = XYCPlot(beam=beamName, xaxis=xaxis, yaxis=yaxis, caxis=caxis, saveName=fname)
    return plot
    
def dump_beam(bl, name, beam):
    n = np.random.randint(1, 100000)
    path = os.path.join(bl.dump_beams_folder, '%s_%06d' % (name, n))
    beamToSave = BeamProxy(beam)
    np.save(path, beamToSave)
    
def do_after_script(bl, plots):
    out_str1 = ""
    fwhm_h = np.array([plot.dx for plot in plots])
    fwhm_v = np.array([plot.dy for plot in plots])
    ecen = np.array([plot.cE/1e3 for plot in plots])
    flux = np.array([plot.flux for plot in plots])
    profile_h = np.array([[p.xaxis.binCenters, p.xaxis.total1D]
                          for p in plots])
    profile_v = np.array([[p.yaxis.binCenters, p.yaxis.total1D]
                          for p in plots])
    profile_e = np.array([[p.caxis.binCenters, p.caxis.total1D]
                          for p in plots])

    epeak = []
    for k, plot in enumerate(plots):
        fname = os.path.join(bl.dump_scores_folder,
                             "res_%s_profile_h.txt" % plot.beam)
        np.savetxt(fname, profile_h[k].T)
        fname = os.path.join(bl.dump_scores_folder,
                             f"res_%s_profile_v.txt" % plot.beam)
        np.savetxt(fname, profile_v[k].T)
        fname = os.path.join(bl.dump_scores_folder,
                             f"res_%s_profile_e.txt" % plot.beam)
        np.savetxt(fname, profile_e[k].T)
        e, y = np.loadtxt(fname, unpack=True)
        epeak.append(e[np.argmax(savgol_filter(y, 51, 3))])
    epeak = np.array(epeak)/1e3

    out_str2 = [
        "epeak  =  " + repr(epeak.tolist()),
        "ecen   =  " + repr(ecen.tolist()),
        "fwhm_h =  " + repr(fwhm_h.tolist()),
        "fwhm_v =  " + repr(fwhm_v.tolist()),
        "flux   =  " + repr(flux.tolist()),
    ]

    out_str2 = chr(10).join(out_str2)
    print(out_str2)
    out_str1 += chr(10) + out_str2 + chr(10)
    with open(os.path.join(bl.dump_scores_folder, "results.txt"), 'w') as f:
        f.write(out_str1)

#
# main
#
def main():
    #
    # xrt beamline
    #
    bl = build_beamline(name="APS_BLM_44keV")
    
    #
    # dumps
    #
    bl.dump_beams_flag = 0
    
    bl.dump_scores_folder = os.path.join(bl.name, "scores")
    os.makedirs(bl.dump_scores_folder, exist_ok=True)
    
    if bl.dump_beams_flag:
        bl.dump_beams_folder = os.path.join(bl.name, "beams")
        os.makedirs(bl.dump_beams_folder, exist_ok=True)
        
    #
    # plots
    #
    if showIn3D:
        bl.glow(scale=[4000, 5, 2000], centerAt='m1')
        return
    
    screens_to_plot = ['SCRexp']
    limitsH = [[-1000.0, 1000.0]]
    limitsV = [[-5000.0, 5000.0]]

    plots = []
    for i in range(len(screens_to_plot)):
        plots.append(make_plot(bl, screens_to_plot[i], limitsH=limitsH[i], limitsV=limitsV[i], bins=1024, climits=None, cbins=256,))
            
    #
    # declare run_process() that makes the tracing
    #
    xrt.backends.raycing.run.run_process = run_process
        
    #
    # run
    #
    run_ray_tracing(beamLine=bl,
                    plots=plots,
                    repeats=1,
                    threads=1,
                    processes=1,
                    afterScript=do_after_script,
                    afterScriptArgs=[bl, plots])


if __name__ in ["__main__", "builtins"]:
    main()
