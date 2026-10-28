
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
    # Component index: 0 (DW100)
    #    
    from xrt.backends.raycing.sources import Wiggler
    bl.DW100 = Wiggler(
        bl,
        name="DW100",
        center=[0,0,0],
        K=16.5,
        period=100.0,
        n=68,
        eE=3.0,
        eI=0.1,
        eEpsilonX=0.9042,
        eEpsilonZ=0.00784,
        eEspread=0.0,
        eSigmaX=137.0,
        eSigmaZ=4.9,
        xPrimeMax=1.1,  # in mrad, to match accepting slit
        zPrimeMax=0.1,  # in mrad, to match accepting slit
        distE="eV",
        eMin=39920.0,
        eMax=40080.0,
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
        opening=[-11000.0/2, 11000.0/2, -1000.0/2, 1000.0/2])                

    #
    # Component index: 2 (Filter1)
    #    
    from xrt.backends.raycing.oes import Plate
    bl.Filter1 = Plate(
        bl,
        name='Filter1',
        center=[0, 10000, 0],
        pitch=np.pi/2,
        material=Material('C', rho=3.52, kind='plate'),
        t=2.0,
        surface='diamond 2 mm')                             

    #
    # Component index: 3 (Filter2)
    #    
    from xrt.backends.raycing.oes import Plate
    bl.Filter2 = Plate(
        bl,
        name='Filter2',
        center=[0, 10004, 0],
        pitch=np.pi/2,
        material=Material(elements=['Si', 'C'], quantities=[1, 1], kind='plate', rho=3.21),
        t=4.8,
        surface='SiC 4.8 mm')                             

    #
    # Component index: 4 (LaueC1)
    #    
    from orangecontrib.xrt.util.bent_crys_2d import BentCrys2D, bentSi_Crys
    from orangecontrib.xrt.util.stiffness_tensor import stiffness_tensor
    CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7, 63.9, 79.6], hkl=[1, 0, 0], sag=[0, -1, 1], unitfac=1e3)
    alpha_hkl_nrm = np.pi / 2 - np.arccos(np.dot([1, 1, 1], [1, 0, 0]) / np.linalg.norm([1, 1, 1]) / np.linalg.norm([1, 0, 0]))
    alpha = alpha_hkl_nrm
    crystal = bentSi_Crys(
        hkl=[1, 1, 1],
        geom='Laue reflected',
        t=0.7,
        aLrad=alpha,
        volumetricDiffraction=False,
        useTT=False,
        Rs=-1575.6831,
        Rm=50670.4309,
        method='ML',
        ML_sign=1.0,
        PP_center_fraction=0.5,
        ST=ST)
    thetaB = crystal.get_Bragg_angle(40000.0) - crystal.get_dtheta_symmetric_Bragg(40000.0)
    bl.LaueC1 = BentCrys2D(
        bl,
        name='LaueC1',
        center=[0, 32300, 0],
        pitch=np.pi/2 + alpha + thetaB,
        alpha=alpha,
        material=(crystal,),
        Rm=50670.4309,
        Rs=-1575.6831,
        Fcurv=1.02,
        limPhysX=[-50, 50],
        limPhysY=[-50, 50],
        )

    #
    # Component index: 5 (LaueC2)
    #    
    from orangecontrib.xrt.util.bent_crys_2d import BentCrys2D, bentSi_Crys
    from orangecontrib.xrt.util.stiffness_tensor import stiffness_tensor
    CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7, 63.9, 79.6], hkl=[1, 0, 0], sag=[0, -1, 1], unitfac=1e3)
    alpha_hkl_nrm = np.pi / 2 - np.arccos(np.dot([1, 1, 1], [1, 0, 0]) / np.linalg.norm([1, 1, 1]) / np.linalg.norm([1, 0, 0]))
    alpha = -alpha_hkl_nrm
    crystal = bentSi_Crys(
        hkl=[1, 1, 1],
        geom='Laue reflected',
        t=0.7,
        aLrad=alpha,
        volumetricDiffraction=False,
        useTT=False,
        Rs=-1575.6831,
        Rm=50670.4309,
        method='ML',
        ML_sign=-1.0,
        PP_center_fraction=0.5,
        ST=ST)
    thetaB = crystal.get_Bragg_angle(40000.0) - crystal.get_dtheta_symmetric_Bragg(40000.0)
    bl.LaueC2 = BentCrys2D(
        bl,
        name='LaueC2',
        center=[0, 32803.8652733817, 50],
        pitch=np.pi/2 + alpha + thetaB,
        alpha=alpha,
        material=(crystal,),
        Rm=50670.4309,
        Rs=-1575.6831,
        Fcurv=1.02,
        limPhysX=[-50, 50],
        limPhysY=[-50, 50],
        )

    #
    # Component index: 6 (M1)
    #    
    from xrt.backends.raycing.oes import ToroidMirror
    bl.M1 = ToroidMirror(
        bl,
        name='M1',
        center=[0, 40100, 50],
        material=Material('Pt', rho=19.305),
        R=13253444.2302,
        r=1e+30,
        pitch=0.0017488571428571423,
        yaw=0.0,
        limPhysX=[-25, 25],
        limPhysY=[-650, 650],
        ) 

    #
    # Component index: 7 (SCRexp)
    #    
    from xrt.backends.raycing.screens import Screen
    bl.SCRexp = Screen(
        bl,
        name="SCRexp",
        center=[0, 56400, 107.01274285714284],
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
    # Component index: 0 (DW100: Wiggler)
    #
    beam = bl.DW100.shine()
    if bl.dump_beams_flag: dump_beam(bl, "DW100", beam)

    #
    # Component index: 1 (Mask1: RectangularAperture)
    #
    beam_local = bl.Mask1.propagate(beam)
    if bl.dump_beams_flag: dump_beam(bl, "Mask1", beam_local)

    #
    # Component index: 2 (Filter1: Plate)
    #
    beam, beam_local1, beam_local2 = bl.Filter1.double_refract(beam)
    if bl.dump_beams_flag: dump_beam(bl, "Filter1_1", beam_local1)
    if bl.dump_beams_flag: dump_beam(bl, "Filter1_2", beam_local2)

    #
    # Component index: 3 (Filter2: Plate)
    #
    beam, beam_local1, beam_local2 = bl.Filter2.double_refract(beam)
    if bl.dump_beams_flag: dump_beam(bl, "Filter2_1", beam_local1)
    if bl.dump_beams_flag: dump_beam(bl, "Filter2_2", beam_local2)

    #
    # Component index: 4 (LaueC1: BentCrys2D)
    #
    beam, _ = bl.LaueC1.reflect(beam)
    if bl.dump_beams_flag: dump_beam(bl, "LaueC1", beam)

    #
    # Component index: 5 (LaueC2: BentCrys2D)
    #
    beam, _ = bl.LaueC2.reflect(beam)
    if bl.dump_beams_flag: dump_beam(bl, "LaueC2", beam)

    #
    # Component index: 6 (M1: ToroidMirrorDistorted)
    #
    beam, _ = bl.M1.reflect(beam)
    if bl.dump_beams_flag: dump_beam(bl, "M1", beam)

    #
    # Component index: 7 (SCRexp: Screen)
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
    bl = build_beamline(name="BNL_SBDL_40keV")
    
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
    limitsH = [[-1500.0, 1500.0]]
    limitsV = [[-500.0, 500.0]]

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
