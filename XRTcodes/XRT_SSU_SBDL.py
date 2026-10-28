

import numpy as np
import matplotlib.pyplot as plt

import xrt.backends.raycing as raycing
import xrt.backends.raycing.sources as rs
import xrt.backends.raycing.apertures as ra
import xrt.backends.raycing.oes as roe
import xrt.backends.raycing.run as rr
import xrt.backends.raycing.materials as rm
import xrt.backends.raycing.materials_compounds as rmc
import xrt.backends.raycing.screens as rsc
import xrt.plotter as xrtp
import xrt.runner as xrtr

from geometry_correction import optimize_geometry_correction
from stiffness_tensor import stiffness_tensor

import BentCrys2D as BC2D

def angle_vec(v1, v2):
    v1_u = v1 / np.linalg.norm(v1)
    v2_u = v2 / np.linalg.norm(v2)
    dot_product = np.dot(v1_u, v2_u)
    dot_product = np.clip(dot_product, -1.0, 1.0)
    return np.arccos(dot_product)

def beam_intensity(beam, good_only=True):
    """
    Return total statistical intensity carried by an xrt beam.
    """
    intensity = np.real(
        np.asarray(beam.Jss) + np.asarray(beam.Jpp)
    )

    good = np.isfinite(intensity)

    if good_only and hasattr(beam, "state"):
        good &= np.asarray(beam.state) > 0

    return np.sum(intensity[good])

def fwhm_from_hist(axis):
    x = axis.binCenters
    y = axis.total1D
    good = np.isfinite(x) & np.isfinite(y)
    x = x[good]
    y = y[good]
    if len(y) == 0:
        return np.nan
    hm = 0.5 * np.max(y)
    above = y >= hm
    if np.sum(above) < 2:
        return np.nan
    i = np.where(above)[0]

    return x[i[-1]] - x[i[0]]

def print_crystal_transmission(name, beam_in, beam_out):
    """
    Print the actual intensity-weighted transmission of one crystal.
    """
    intensity_in = beam_intensity(beam_in)
    intensity_out = beam_intensity(beam_out)

    transmission = (
        intensity_out / intensity_in
        if intensity_in > 0.0
        else np.nan
    )

    print(
        f"{name}: "
        f"Iin={intensity_in:.8e}, "
        f"Iout={intensity_out:.8e}, "
        f"T={transmission:.8e}"
    )

    return transmission

# ==========================================================================
# beamline parameters

# source parameters

bener = 3. #ring energy [GeV]
# e-beam size and div RMS
esizeh = 19.3
esizev = 1.9
edivh = 12.9
edivv = 1.3
B0 = 6.6 # T

delta_E = 0.05

# beam acceptance
acch = 1.6 #~[mrad]
accv = 0.66 # [mrad]

# distances
d_slit1 = 10.0 * 1e3 # [mm] distance slits and filters
d_dblm_c1 = 31. * 1e3 # [mm] distance source - SBDLM C1
d_dblm_c2 = 32. * 1e3 # [mm] distance source - SBDLM C2
#d_exp = 44. * 1e3 # [mm] distance source - exp1
d_exp = 56. * 1e3 # [mm] distance source - exp2
fixedexit = 50. # [mm]

# DBLM crystals
c_thick = 1.0 #[mm]
crystal = 'Si' # crystal Si(111) on (100) substrate 
hkl = [1,1,1]
nrm = [1,0,0]
sag = [0,1,1]
h, k, l = hkl
s1, s2, s3 = sag
str_hkl = 'Si' + "".join(map(str, hkl))
str_nrm = "".join(map(str, nrm))
str_sag = "".join(map(str, sag))

aLcalc = angle_vec(hkl, nrm)
aLrad = np.pi/2. - aLcalc
print('aL',np.rad2deg(aLrad))

# standard length unit [mm] - unitfac=1e3
unitfac = 1e3
# mechanical properties of Laue crystals

c_fac = 0.5
CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7,63.9,79.6], hkl=nrm, sag=sag, unitfac=unitfac)
poL =  - ST[1,2] / ST[2,2] # ST'(ij)= ST(ij) / ST(33) - poL=ST'(23)
Sf = 'S_tensor.txt'
np.savetxt(Sf, ST)
print('pol',poL)

# energy, scan parameters
energies = np.linspace(35000, 200000, 34)

# simulation parameters
repeats = 1
nrays = 500000

title0 = 'SSU_SBDL_56m_' + str_hkl + '_nrm' + str_nrm + '_sag' + str_sag +  '_PP_'

# ==================================================================
def build_beamline():
    height = 0.0
   
    beamLine = raycing.BeamLine(azimuth=0.0, height=height, alignE='auto')
   
    # wiggler source
    beamLine.Source = rs.BendingMagnet(beamLine, name='SWLS', center=(0,0,0),
        nrays=nrays, B0=B0, eE=bener, eI=0.1,
        eSigmaX=esizeh, eSigmaZ=esizev, eEpsilonX=esizeh*edivh*1e-3, eEpsilonZ=esizev*edivv*1e-3,
        eMin=energies[0]*(1-delta_E),eMax=energies[0]*(1+delta_E),
        xPrimeMax=acch, zPrimeMax=accv)
    # ----------------------------------------------------------------------------------
    # aperture and filters @ 10m
    s1_xwidth = d_slit1 * acch * 1e-3
    s1_zwidth = d_slit1 * accv * 1e-3
    filterDiamond = rm.Material('C', rho=3.52, kind='plate')
    filterSi = rm.Material('Si', rho=2.33, kind='plate')

    beamLine.Mask1 = ra.RectangularAperture(beamLine, 'Mask1', [0., d_slit1, height],
        ('left', 'right', 'bottom', 'top'), [-s1_xwidth/2., s1_xwidth/2., -s1_zwidth/2., s1_zwidth/2.])
    beamLine.Screen1 = rsc.Screen(beamLine, 'SCR1', [0., d_slit1, height])
    beamLine.Filter1 = roe.Plate(beamLine, 'Filter1', (0., d_slit1, height), pitch=np.pi/2,
        limPhysX=(-s1_xwidth/2., s1_xwidth/2.), limPhysY=(-s1_zwidth/2., s1_zwidth/2.),
        material=(filterDiamond,), t=0.05)
    #beamLine.Filter2 = roe.Plate(beamLine, 'Filter2', (0., d_slit1, height), pitch=np.pi/2,
    #    limPhysX=(-s1_xwidth/2., s1_xwidth/2.), limPhysY=(-s1_zwidth/2., s1_zwidth/2.),
    #    material=(filterSi,), t=0.4)
    beamLine.Screen2 = rsc.Screen(beamLine, 'SCR2', [0., d_slit1, height])

    beamLine.LaueC1 = BC2D.BentCrys2D(
        name='LaueC1',
        bl=beamLine,
        alpha=aLrad,
        limPhysX=[-50.,50.],
        limPhysY=[-50.,50.]
    )
    beamLine.LaueC2 = BC2D.BentCrys2D(
        name='LaueC2',
        bl=beamLine,
        alpha=-aLrad,
        limPhysX=[-50.,50.],
        limPhysY=[-50.,50.]
    )

    beamLine.SCRexp = rsc.Screen(beamLine, 'SCRexp', [0., d_exp, fixedexit])
    # -----------------------------------------------------------------------------------

    return beamLine

def run_process(beamLine):
    """
    Double-Laue beamline with:
    1. one upstream trace;
    2. common Fcurv optimization using a reduced beam and ML;
    3. one final full PP or ML trace through C1-C2.
    """
    # ================================================================
    # Settings
    # ================================================================

    optimization_ray_count = 50000
    optimization_seed = 12345
    # Geometry optimization is faster and more stable with ML.
    optimization_method = "ML"
    initial_Fcurv = 1.15
    initial_step = 0.05
  
    # ================================================================
    # Helpers
    # ================================================================

    def copy_beam_subset(input_beam, indices):
        """
        Copy an xrt beam while retaining only selected ray indices.
        """
        output_beam = rs.Beam(copyFrom=input_beam)
        number_of_rays = len(input_beam.x)

        for attribute_name, attribute_value in vars(output_beam).items():
            if not isinstance(attribute_value, np.ndarray):
                continue
            if attribute_value.ndim == 0:
                continue
            if attribute_value.shape[0] != number_of_rays:
                continue
            setattr(output_beam,attribute_name,np.array(attribute_value[indices],copy=True,))

        return output_beam


    def make_reduced_valid_beam(input_beam,maximum_rays,seed):
        """
        Return a deterministic subset containing valid rays only.
        """
        if hasattr(input_beam, "state"):
            valid_indices = np.flatnonzero(np.asarray(input_beam.state) > 0)
        else:
            valid_indices = np.arange(len(input_beam.x))

        if valid_indices.size == 0:
            raise RuntimeError("No valid rays are available.")

        number_to_keep = min(int(maximum_rays),valid_indices.size,)

        if number_to_keep < valid_indices.size:
            rng = np.random.default_rng(seed)
            selected_indices = rng.choice(valid_indices,size=number_to_keep,replace=False)
            selected_indices = np.sort(selected_indices)
        else:
            selected_indices = valid_indices

        return copy_beam_subset(input_beam,selected_indices)

    def horizontal_fwhm(values, weights, bins=300):
        values = np.asarray(values, dtype=float)
        weights = np.asarray(weights, dtype=float)
        good = (
            np.isfinite(values)
            & np.isfinite(weights)
            & (weights > 0.0)
        )
        values = values[good]
        weights = weights[good]
        if values.size < 2:
            return np.nan
        hist, edges = np.histogram(values,bins=bins,weights=weights)
        if (hist.size < 2 or not np.any(hist > 0.0)):
            return np.nan
        centers = 0.5 * (edges[:-1] + edges[1:])
        peak_index = int(np.argmax(hist))
        half = 0.5 * hist[peak_index]
        # connected peak containing maximum
        il = peak_index
        while (il > 0 and hist[il - 1] >= half):
            il -= 1
        ir = peak_index
        while (ir < len(hist) - 1 and hist[ir + 1] >= half):
            ir += 1
        if il == ir:
            return float(edges[1] - edges[0])
        return float(centers[ir] - centers[il])

    def horizontal_caustic(beam,y_positions,bins=300,minimum_rays=100):
        x0 = np.asarray(beam.x,dtype=float)
        y0 = np.asarray(beam.y,dtype=float)
        a = np.asarray(beam.a,dtype=float)
        b = np.asarray(beam.b,dtype=float)
        weights = np.real(np.asarray(beam.Jss) + np.asarray(beam.Jpp))
        good = (
            (np.asarray(beam.state) > 0)
            & np.isfinite(x0)
            & np.isfinite(y0)
            & np.isfinite(a)
            & np.isfinite(b)
            & np.isfinite(weights)
            & (weights > 0.0)
            & (np.abs(b) > 1e-15)
        )
        x0 = x0[good]
        y0 = y0[good]
        a = a[good]
        b = b[good]
        weights = weights[good]
        if x0.size < minimum_rays:
            raise RuntimeError("Too few rays for horizontal caustic.")
        widths = np.full(len(y_positions),np.nan)
        for i, y_plane in enumerate(y_positions):
            tau = (y_plane - y0) / b
            x_plane = (x0 + a * tau)
            widths[i] = horizontal_fwhm(x_plane,weights,bins=bins)
        finite = np.isfinite(widths)
        if not np.any(finite):
            return np.nan, np.nan
        valid_indices = np.flatnonzero(finite)
        i_min = valid_indices[np.argmin(widths[finite])]
        return (float(y_positions[i_min]),float(widths[i_min]))
        
    # ================================================================
    # Trace all fixed upstream elements once
    # ================================================================

    beamSource = beamLine.Source.shine()
    beamMask1 = beamLine.Mask1.propagate(beamSource)
    beamSCR1 = beamLine.Screen1.expose(beamMask1)

    beamFilter1Global,beamFilter1Local1,beamFilter1Local2 = beamLine.Filter1.double_refract(beamMask1)
    beamSCR2 = beamLine.Screen2.expose(beamFilter1Global)
    # Fixed input immediately before C1.
    full_input_before_C1 = rs.Beam(copyFrom=beamFilter1Global)

    # ================================================================
    # Prepare reduced beam for Fcurv optimization
    # ================================================================

    reduced_input_before_C1 = make_reduced_valid_beam(
        input_beam=full_input_before_C1,
        maximum_rays=optimization_ray_count,
        seed=optimization_seed,
    )
    material_C1 = beamLine.LaueC1.material[0]
    material_C2 = beamLine.LaueC2.material[0]
    requested_method_C1 = material_C1.method
    requested_method_C2 = material_C2.method

    # ================================================================
    # Common-Fcurv trial
    # ================================================================
    def trace_Fcurv_trial(Fcurv):
        """
        Apply one common Fcurv to C1 and C2 and trace the fixed
        reduced input beam.
        """
        Fcurv = float(Fcurv)
        beamLine.LaueC1.Fcurv = Fcurv
        beamLine.LaueC2.Fcurv = Fcurv
        trial_input = rs.Beam(copyFrom=reduced_input_before_C1)
        beamC1Trial, _ = beamLine.LaueC1.reflect(beam=trial_input)
        beamC2Trial, _ = beamLine.LaueC2.reflect(beam=beamC1Trial)
        payload = {"Fcurv": Fcurv}

        return beamC2Trial, payload

    
    # ==============================================================
    # Fcurv optimization
    # First energy:
    #     broad optimization -> fine optimization
    # Following energies:
    #     previous optimum -> fine optimization only
    # ==============================================================

    try:
        material_C1.method = optimization_method
        material_C2.method = optimization_method
        # 1. BROAD SEARCH -- ONLY FIRST ENERGY
        if not hasattr(beamLine, "_previous_Fcurv"):
            print("\nBroad Fcurv optimization...")
            correction_broad = optimize_geometry_correction(
                trace_trial=trace_Fcurv_trial,
                plane_point=beamLine.SCRexp.center,
                plane_normal=(0.0, 1.0, 0.0),
                transverse_axis=(1.0, 0.0, 0.0),
                # broad starting search
                initial_Fcurv=1.15,
                initial_step=0.05,
                bounds=(0.5, 1.5),
                max_bracket_expansions=5,
                # only coarse refinement needed here
                refinement_iterations=2,
                grid_points=5,
                Fcurv_tolerance=1e-3,
                bin_width=0.002,
                use_intensity_weights=True,
                minimum_rays=500,
                verbose=False
            )
            fine_center = float(correction_broad.Fcurv)
            print(
                "\nBroad Fcurv result:",
                f"{fine_center:.9f}"
            )
        # FOLLOWING ENERGIES
        else:
            fine_center = float(beamLine._previous_Fcurv)
            print(
                "\nPrevious-energy Fcurv:",
                f"{fine_center:.9f}"
            )
        # 2. FINE SEARCH -- EVERY ENERGY
        fine_half_range = 0.015
        fine_lower = (fine_center - fine_half_range)
        fine_upper = (fine_center + fine_half_range)
        print(
            "Fine Fcurv search:",
            f"{fine_lower:.9f}",
            "to",
            f"{fine_upper:.9f}"
        )
        correction = optimize_geometry_correction(
            trace_trial=trace_Fcurv_trial,
            plane_point=beamLine.SCRexp.center,
            plane_normal=(0.0, 1.0, 0.0),
            transverse_axis=(1.0, 0.0, 0.0),
            initial_Fcurv=fine_center,
            # much smaller initial step
            initial_step=0.005,
            bounds=(
                fine_lower,
                fine_upper
            ),
            # no large bracket expansion required
            max_bracket_expansions=2,
            refinement_iterations=3,
            grid_points=5,
            Fcurv_tolerance=1e-4,
            bin_width=0.002,
            use_intensity_weights=True,
            minimum_rays=500,
            verbose=False
        )
        # Store optimum for next energy
        best_Fcurv = float(correction.Fcurv)
        beamLine._previous_Fcurv = (best_Fcurv)

    finally:
        material_C1.method = requested_method_C1
        material_C2.method = requested_method_C2

    best_Fcurv = float(correction.Fcurv)
    beamLine.LaueC1.Fcurv = best_Fcurv
    beamLine.LaueC2.Fcurv = best_Fcurv

    print("Optimized double-Laue correction:",f"Fcurv={best_Fcurv:.9f}",)

    beamC1Input = rs.Beam(copyFrom=full_input_before_C1)
    beamC1Global, beamC1Local = (beamLine.LaueC1.reflect(beam=beamC1Input))

    beamC2Input = rs.Beam(copyFrom=beamC1Global)
    beamC2Global, beamC2Local = (beamLine.LaueC2.reflect(beam=beamC2Input))

    beamSCRexp = beamLine.SCRexp.expose(beamC2Global)

    # ================================================================
    # Horizontal caustic
    # Broad scan: first energy only
    # Fine scan: every energy
    # ================================================================
    nominal_focus_y = float(beamLine.SCRexp.center[1])
    # BROAD horizontal caustic -- only first energy
    if not hasattr(beamLine,"_coarse_horizontal_focus"):
        print("\nCalculating broad horizontal caustic...")
        coarse_positions = np.arange(nominal_focus_y - 12000.0,nominal_focus_y + 3000.0 + 1.0,100.0)
        (coarse_focus,coarse_width) = horizontal_caustic(beam=beamC2Global,y_positions=coarse_positions,bins=200)
        beamLine._coarse_horizontal_focus = (coarse_focus)
        print(
            "Broad horizontal focus:",
            f"{coarse_focus:.3f} mm",
            "FWHM:",
            f"{coarse_width:.6e} mm"
        )
    # FINE horizontal caustic -- every energy
    if hasattr(beamLine,"_previous_horizontal_focus"):
        fine_center = float(beamLine._previous_horizontal_focus)
    else:
        fine_center = float(beamLine._coarse_horizontal_focus)
    fine_positions = np.arange(fine_center - 200.0,fine_center + 200.0 + 1.0,1.0)
    (horizontal_focus,horizontal_width) = horizontal_caustic(beam=beamC2Global,y_positions=fine_positions,bins=300)
    edge_margin = 10.0
    if (horizontal_focus <= fine_positions[0] + edge_margin or horizontal_focus >= fine_positions[-1] - edge_margin):
        print(
            "Horizontal focus near fine-scan boundary; "
            "expanding caustic."
        )
        expanded_positions = np.arange(fine_center - 1000.0,fine_center + 1000.0 + 2.0,2.0)
        (horizontal_focus,horizontal_width) = horizontal_caustic(beam=beamC2Global,y_positions=expanded_positions,bins=300)

    # store for next energy
    beamLine._previous_horizontal_focus = (horizontal_focus)
    print(
        "Horizontal focus:",
        f"{horizontal_focus:.6f} mm"
    )
    print(
        "Minimum horizontal FWHM:",
        f"{horizontal_width:.6e} mm"
    )
    print(
        "Horizontal focus offset from SCRexp:",
        f"{horizontal_focus - nominal_focus_y:.3f} mm"
    )

    # ================================================================
    # Return beams required by define_plots()
    # ================================================================

    return {
        "beamSource": beamSource,
        "beamMask1": beamMask1,
        "beamSCR1": beamSCR1,
        "beamFilter1Global": beamFilter1Global,
        "beamSCR2": beamSCR2,
        "beamC1Global": beamC1Global,
        "beamC1Local": beamC1Local,
        "beamC2Global": beamC2Global,
        "beamC2Local": beamC2Local,
        "beamSCRexp": beamSCRexp,
    }

rr.run_process = run_process


def define_plots(beamLine):
    plots = []
    # at mask @ 10m - only aperture
    plot_s = xrtp.XYCPlot('beamSource', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.2f'),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.2f'),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f'),
        fluxKind='total',
        title='source_div')
    plots.append(plot_s)
    # after filtering @ 10m
    plot_src = xrtp.XYCPlot('beamSCR2', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x$", r"mm", limits=None, fwhmFormatStr='%.3f'),
        yaxis=xrtp.XYCAxis(r"$z$", r"mm", limits=None, fwhmFormatStr='%.3f'),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f'),
        fluxKind='total',
        title='SCR2')
    plots.append(plot_src)

    plot_C1div = xrtp.XYCPlot('beamC1Global', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f'),
        fluxKind='total',
        title='LaueC1')
    plots.append(plot_C1div)

    # after at focus
    plot_bs = xrtp.XYCPlot('beamSCRexp', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x$", r"mm", limits=[-0.1,0.1], fwhmFormatStr='%.3f'),
        yaxis=xrtp.XYCAxis(r"$z$", r"mm", limits=None, fwhmFormatStr='%.3f'),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f'),
        fluxKind='total',
        title=title0+'SCRexp1-bs')
    plots.append(plot_bs)
    plot_div = xrtp.XYCPlot('beamSCRexp', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f'),
        fluxKind='total',
        title=title0+'SCRexp1-div')
    plots.append(plot_div)
   
    return plots

def update_energy_axes(plots, energy_eV, relative_bandwidth):
    """
    Update every plot color-energy axis for the current scan energy.
    XYCAxis internal energy values and limits are in eV, even when the
    displayed unit is keV. The axis.factor performs the display conversion.
    """
    energy_keV = energy_eV * 1e-3
    e_min = energy_keV * (1.0 - relative_bandwidth)
    e_max = energy_keV * (1.0 + relative_bandwidth)

    for plot in plots:
        caxis = plot.caxis
        # Only modify axes representing photon energy.
        if "energy" not in caxis.label.lower():
            continue
        caxis.limits = [e_min, e_max]
        # Show a narrow scale relative to the nominal energy.
        caxis.offset = energy_keV 

def update_scan_plot_axes(plots):
    """
    Reset or explicitly update spatial and angular plot limits
    before each scan point.

    Axis limits are expressed in the displayed plot units:
        position: mm
        divergence: mrad
    """

    # plots[3]: beam size at experimental screen
    plots[3].xaxis.limits = [-0.1,0.1]
    plots[3].yaxis.limits = None

    # plots[4]: divergence at experimental screen
    plots[4].xaxis.limits = None
    plots[4].yaxis.limits = None

def energy_scan(beamLine, plots, energies):
    flux0 = np.zeros_like(energies)
    flux = np.zeros_like(energies)
    flux_dens = np.zeros_like(energies)
    beam_x = np.zeros_like(energies)
    beam_z = np.zeros_like(energies)
    beam_dx = np.zeros_like(energies)
    beam_dz = np.zeros_like(energies)
    beam_dE = np.zeros_like(energies)
    fOut = open('res_' + title0 + 'energy.txt', 'w')
    fOut.write("E_eV\tflux0\tflux\tflux_density\tbeam_x\tbeam_z\tbeam_dx\tbeam_dz\tdE\tRs\tRm\n")
    for ie, e in enumerate(energies):
        print(f'energy {e:.1f} eV, {ie+1} of {len(energies)}')
        # thetaB and distnace between crystals
        # Bragg angle corrected for refraction (xrt)
        crystal_Bragg = rm.CrystalSi(hkl=hkl)
        thetaB = crystal_Bragg.get_Bragg_angle(e) - crystal_Bragg.get_dtheta_symmetric_Bragg(e)  # radians
        d_LaueC1_C2 = fixedexit / np.tan(2. * thetaB)

        # ==========================================
        # Focusing (Laue bending) - Rs, Rm
        pL = d_dblm_c1
        qL = d_exp - d_dblm_c1
        fs = 1 / (1 / pL + 1 / qL)
        Rs = 4 * np.sin(thetaB) * np.sin(aLrad) * fs
        Rm = Rs / c_fac / poL
        # crystal material parameters
        Cryst1 = BC2D.bentSi_Crys(hkl=hkl,geom='Laue Reflected',t=c_thick, aLrad=aLrad, 
                                    volumetricDiffraction=False,useTT=False,
                                    Rs=-Rs, Rm=Rm, method='PP',ML_sign=1.0,
                                    PP_center_fraction=0.5,ST=ST)
        Cryst2 = BC2D.bentSi_Crys(hkl=hkl,geom='Laue Reflected',t=c_thick, aLrad=-aLrad,
                                    volumetricDiffraction=False,useTT=False,
                                    Rs=-Rs, Rm=Rm, method='PP',ML_sign=-1.0,
                                    PP_center_fraction=0.5,ST=ST)

        pitch1 = np.pi/2 + aLrad + thetaB
        pitch2 = np.pi/2 - aLrad + thetaB
    
        Fcurv = 1.0

        beamLine.LaueC1.material = (Cryst1,)
        beamLine.LaueC1.Rs = -Rs
        beamLine.LaueC1.Rm = Rm
        beamLine.LaueC1.pitch = pitch1
        beamLine.LaueC1.Fcurv = Fcurv
        beamLine.LaueC1.center = [0, d_dblm_c1, 0]

        beamLine.LaueC2.material = (Cryst2,)
        beamLine.LaueC2.Rs = -Rs
        beamLine.LaueC2.Rm = Rm
        beamLine.LaueC2.pitch = pitch2
        beamLine.LaueC2.Fcurv = Fcurv
        beamLine.LaueC2.center = [0, d_dblm_c1 + d_LaueC1_C2, fixedexit]

        beamLine.alignE = e
        beamLine.Source.eMin = e * (1 - delta_E)
        beamLine.Source.eMax = e * (1 + delta_E)
        if hasattr(beamLine.Source, "reset"):
            beamLine.Source.reset()

        for plot in plots:
            plot.saveName = [plot.title + f'_{e/1000:.0f}keV.png']
            plot.caxis.limits = [e/1000*(1-delta_E), e/1000*(1+delta_E)]

        update_energy_axes(
            plots=plots,
            energy_eV=e,
            relative_bandwidth=delta_E,
        )

        update_scan_plot_axes(plots=plots)
                  
        for plot in plots[-2:]:
            plot.saveName = (
                f"{plot.title}_{e / 1000.0:.1f}keV.png"
            )
        
        yield

        psrc = plots[1]
        ps = plots[3]
        pd = plots[4]
        sx95 = fwhm_from_hist(ps.xaxis)
        sz95 = fwhm_from_hist(ps.yaxis)
        dx95 = fwhm_from_hist(pd.xaxis)
        dz95 = fwhm_from_hist(pd.yaxis)
        flux0[ie] = psrc.flux
        flux[ie] = ps.flux
        flux_dens[ie] = ps.flux / ps.dx / ps.dy
        beam_x[ie] = sx95
        beam_z[ie] = sz95
        beam_dx[ie] = dx95
        beam_dz[ie] = dz95
        beam_dE[ie] = ps.dE
        outStr = (
            f"{e:.6f}\t"
            f"{flux0[ie]:.8e}\t"
            f"{flux[ie]:.8e}\t"
            f"{flux_dens[ie]:.8e}\t"
            f"{beam_x[ie]:.8e}\t"
            f"{beam_z[ie]:.8e}\t"
            f"{beam_dx[ie]:.8e}\t"
            f"{beam_dz[ie]:.8e}\t"
            f"{beam_dE[ie]:.8e}\t"
            f"{Rs:.6f}\t"
            f"{Rm:.6f}\n"
        )
        print(outStr)
        fOut.write(outStr)
        fOut.flush()
    
    fOut.close()

    plt.plot(energies, flux)
    plt.show()

def main():
    beamLine = build_beamline()
    plots = define_plots(beamLine)
    xrtr.run_ray_tracing(plots=plots, repeats=repeats, beamLine=beamLine, 
                         generator=energy_scan, generatorArgs=[beamLine, plots, energies])
    
if __name__ == '__main__':
    main()