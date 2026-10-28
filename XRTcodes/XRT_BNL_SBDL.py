

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

def beam_width(values,weights=None,statistic="FWHM",bins=300,bin_width=None,qlo=0.16,qhi=0.84):
    """
    Calculate a 1D beam width from ray coordinates.
    Parameters
    ----------
    values : array-like
        Ray coordinates, for example x or z positions in mm.
    weights : array-like or None
        Statistical ray weights, normally Jss + Jpp.
        If None, all rays have equal weight.
    statistic : {"FWHM", "RMS", "95%"}
        Width definition:
        FWHM Full width at half maximum of a weighted histogram.
        RMS Standard deviation sigma of the weighted ray distribution.
        95% Central quantile width from qlo to qhi.
            With qlo=0.025 and qhi=0.975, this contains 95% of
            the statistical beam intensity.
    bins : int or sequence
        Number of histogram bins, used only for FWHM.
    qlo, qhi : float
        Lower and upper quantiles used for statistic="95%".

    Returns
    -------
    width : float
        Requested width. For RMS, this is one sigma.

    left : float
        Left boundary:
        - half-maximum crossing for FWHM
        - mean - sigma for RMS
        - lower quantile for 95%

    right : float
        Right boundary:
        - half-maximum crossing for FWHM
        - mean + sigma for RMS
        - upper quantile for 95%

    center : float
        Weighted mean position.

    Notes
    -----
    For RMS, the returned width is sigma, not 2*sigma.

    FWHM is histogram-dependent. Increasing `bins` improves spatial
    resolution but also increases statistical noise.
    """
    if bin_width is None:
        histogram, edges = np.histogram(values,bins=bins,weights=weights)
    else:
        xmin = np.min(values)
        xmax = np.max(values)
        left_edge = (np.floor(xmin / bin_width) * bin_width)
        right_edge = (np.ceil(xmax / bin_width) * bin_width)
        edges = np.arange(left_edge,right_edge + 1.001 * bin_width,bin_width)
        histogram, edges = np.histogram(values,bins=edges,weights=weights)

    values = np.asarray(values, dtype=float)
    if weights is None:
        weights = np.ones_like(values, dtype=float)
    else:
        weights = np.asarray(weights, dtype=float)
    if values.shape != weights.shape:
        raise ValueError("values and weights must have the same shape.")
    good = (np.isfinite(values) & np.isfinite(weights) & (weights > 0.0))
    values = values[good]
    weights = weights[good]
    if values.size < 2:
        return np.nan, np.nan, np.nan, np.nan
    weight_sum = np.sum(weights)
    if not np.isfinite(weight_sum) or weight_sum <= 0.0:
        return np.nan, np.nan, np.nan, np.nan
    center = np.sum(weights * values) / weight_sum
    statistic = statistic.strip().upper()

    # RMS width
    if statistic == "RMS":
        variance = (np.sum(weights * (values - center)**2) / weight_sum)
        variance = max(float(variance), 0.0)
        sigma = np.sqrt(variance)
        return (sigma,center - sigma,center + sigma,center)

    # Central quantile width
    if statistic in ("95%", "PERCENTILE", "QUANTILE"):
        if not (0.0 <= qlo < qhi <= 1.0):
            raise ValueError("Require 0 <= qlo < qhi <= 1.")
        order = np.argsort(values)
        x_sorted = values[order]
        w_sorted = weights[order]
        cumulative = np.cumsum(w_sorted)
        cumulative /= cumulative[-1]
        left = np.interp(qlo,cumulative,x_sorted)
        right = np.interp(qhi,cumulative,x_sorted)
        return (right - left,left,right,center)

    # Histogram FWHM
    if statistic == "FWHM":
        from scipy.ndimage import gaussian_filter1d
        # ----------------------------------------------------------
        # Fixed physical histogram resolution
        # ----------------------------------------------------------
        bin_width = 0.01   # mm for beam-position distributions
        xmin = np.min(values)
        xmax = np.max(values)
        left_edge = np.floor(xmin / bin_width) * bin_width
        right_edge = np.ceil(xmax / bin_width) * bin_width
        edges = np.arange(left_edge,right_edge + 1.001 * bin_width,bin_width)
        histogram, edges = np.histogram(values,bins=edges,weights=weights)
        histogram = np.asarray(histogram,dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        if (histogram.size < 2 or not np.any(histogram > 0.0)):
            return np.nan, np.nan, np.nan, center
        # ----------------------------------------------------------
        # Suppress single-bin statistical spikes
        # sigma=1 means only mild smoothing
        # ----------------------------------------------------------
        histogram_smooth = gaussian_filter1d(histogram,sigma=1.5)
        peak_index = int(np.argmax(histogram_smooth))
        half_maximum = (0.5 * histogram_smooth[peak_index])
        # ----------------------------------------------------------
        # Connected FWHM around physical main peak
        # ----------------------------------------------------------
        left_index = peak_index
        while (left_index > 0 and histogram_smooth[left_index - 1] >= half_maximum):
            left_index -= 1
        right_index = peak_index
        while (right_index < histogram_smooth.size - 1 and histogram_smooth[right_index + 1] >= half_maximum):
            right_index += 1
        # ----------------------------------------------------------
        # Interpolated left crossing
        # ----------------------------------------------------------
        if left_index == 0:
            left = edges[0]
        else:
            x1 = centers[left_index - 1]
            x2 = centers[left_index]
            y1 = histogram_smooth[left_index - 1]
            y2 = histogram_smooth[left_index]
            if y2 == y1:
                left = x2
            else:
                left = (x1 + (half_maximum - y1) * (x2 - x1) / (y2 - y1))
        # ----------------------------------------------------------
        # Interpolated right crossing
        # ----------------------------------------------------------
        if right_index == histogram_smooth.size - 1:
            right = edges[-1]
        else:
            x1 = centers[right_index]
            x2 = centers[right_index + 1]
            y1 = histogram_smooth[right_index]
            y2 = histogram_smooth[right_index + 1]
            if y2 == y1:
                right = x1
            else:
                right = (x1 + (half_maximum - y1) * (x2 - x1) / (y2 - y1))
        return (right - left,left,right,center)

    raise ValueError("Unknown statistic. Use 'FWHM', 'RMS' or '95%'.")

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
esizeh = 137.0
esizev = 4.9
edivh = 6.6
edivv = 1.6
#espread = 7.99e-4
Nperiod = 68
perlen = 100.0
keff = 16.5

delta_E = 0.002

# beam acceptance
acch = 1.1 #~[mrad]
accv = 0.1 # [mrad]

# distances
d_slit1 = 10.0 * 1e3 # [mm] distance slits and filters
d_dblm_c1 = 32.3 * 1e3 # [mm] distance source - SBDLM C1
d_dblm_c2 = 32.5 * 1e3 # [mm] distance source - SBDLM C2
d_m1 = 40.1 * 1e3 # [mm] distance source - first mirror
d_exp = 56.4 * 1e3 # [mm] distance source - exp
fixedexit = 50. # [mm]

# DBLM crystals
c_thick = 0.7 #[mm]
crystal = 'Si' # crystal Si(111) on (100) substrate 
hkl = [1,1,1]
nrm = [1,0,0]
sag = [0,-1,1]

h, k, l = hkl
str_hkl = 'Si' + "".join(map(str, hkl))
str_nrm = "".join(map(str, nrm))
str_sag = "".join(map(str, sag))

#Mirror M1
M1len = 1300.0
M1width = 50.0
M1mat = rm.Material('Pt', rho=19.305)

aLcalc = angle_vec(hkl, nrm)
aLrad = np.pi/2. - aLcalc
print('aL',np.rad2deg(aLrad))

# standard length unit [mm] - unitfac=1e3
unitfac = 1e3
# mechanical properties of Laue crystals

c_fac = 0.5
CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7,63.9,79.6], hkl=nrm, sag=sag, unitfac=unitfac)
poL = - ST[1,2] / ST[2,2]  # ST'(ij)= ST(ij) / ST(33) - poL=ST'(23)
Sf = 'S_tensor.txt'
np.savetxt(Sf, ST)
print('pol',poL)

# energy, scan parameters
energies = np.linspace(30000, 70000, 5)

# simulation parameters
repeats = 1
nrays = 500000

title0 = 'BNL_SBDL_highflux' + str_hkl + '_nrm' + str_nrm + '_sag' + str_sag +  '_ML_'

# ==================================================================
def build_beamline():
    height = 0.0
   
    beamLine = raycing.BeamLine(azimuth=0.0, height=height, alignE='auto')
   
    # dumping wiggler source
    beamLine.Source = rs.Wiggler(beamLine, name='DW100', center=(0,0,0),
        nrays=nrays, K=keff, eE=bener, eI=0.1, n=Nperiod, period=perlen, 
        eSigmaX=esizeh, eSigmaZ=esizev, eEpsilonX=esizeh*edivh*1e-3, eEpsilonZ=esizev*edivv*1e-3,
        eMin=energies[0]*(1-delta_E),eMax=energies[0]*(1+delta_E),
        xPrimeMax=acch, zPrimeMax=accv)
    # ----------------------------------------------------------------------------------
    # aperture and filters @ 10m
    s1_xwidth = d_slit1 * acch #* 1e-3
    s1_zwidth = d_slit1 * accv #* 1e-3

    beamLine.Mask1 = ra.RectangularAperture(beamLine, 'Mask1', [0., d_slit1, height],
        ('left', 'right', 'bottom', 'top'), [-s1_xwidth/2., s1_xwidth/2., -s1_zwidth/2., s1_zwidth/2.])
    beamLine.Screen1 = rsc.Screen(beamLine, 'SCR1', [0., d_slit1, height])

    filterDiamond = rm.Material('C',rho=3.52,kind='plate')
    filterSiC = rm.Material(
        elements=['Si', 'C'], 
        quantities=[1, 1], 
        kind='plate', 
        rho=3.21  # Dichte von SiC in g/cm³
    )
    beamLine.Filter1 = roe.Plate(
        beamLine,
        'Filter1',
        (0., d_slit1, height),
        pitch=np.pi/2,
        limPhysX=(-s1_xwidth/2., s1_xwidth/2.),
        limPhysY=(-s1_zwidth/2., s1_zwidth/2.),
        material=filterDiamond,
        t=2.
    )
    beamLine.Filter2 = roe.Plate(
        beamLine,
        'Filter2',
        (0., d_slit1+4, height),
        pitch=np.pi/2,
        limPhysX=(-s1_xwidth/2., s1_xwidth/2.),
        limPhysY=(-s1_zwidth/2., s1_zwidth/2.),
        material=filterSiC,
        t=4.8
    )
    
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

    beamLine.M1 = roe.ToroidMirror(
        beamLine, 'M1',r=1e30,
        surface=('Pt',), material=(M1mat,),
        limPhysX=(-M1width/2, M1width/2), limPhysY=(-M1len/2, M1len/2), limOptX=(-M1width/2, M1width/2),
        limOptY=(-M1len/2, M1len/2), alarmLevel=0.)

    beamLine.SCRexp = rsc.Screen(beamLine, 'SCRexp')
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

    def directional_caustic(beam,y_positions,direction="horizontal",bin_width=0.002,minimum_rays=100):
        y0 = np.asarray(beam.y,dtype=float)
        b = np.asarray(beam.b,dtype=float)
        if direction.lower() == "horizontal":
            pos0 = np.asarray(beam.x,dtype=float)
            slope = np.asarray(beam.a,dtype=float)
        elif direction.lower() == "vertical":
            pos0 = np.asarray(beam.z,dtype=float)
            slope = np.asarray(beam.c,dtype=float)
        else:
            raise ValueError("direction must be " "'horizontal' or 'vertical'")
        weights = np.real(np.asarray(beam.Jss) + np.asarray(beam.Jpp))
        good = (
            (np.asarray(beam.state) > 0)
            & np.isfinite(pos0)
            & np.isfinite(y0)
            & np.isfinite(slope)
            & np.isfinite(b)
            & np.isfinite(weights)
            & (weights > 0.0)
            & (np.abs(b) > 1e-15)
        )
        pos0 = pos0[good]
        y0 = y0[good]
        slope = slope[good]
        b = b[good]
        weights = weights[good]
        if pos0.size < minimum_rays:
            raise RuntimeError(f"Too few rays for {direction} caustic.")
        widths = np.full(len(y_positions),np.nan)
        for i, y_plane in enumerate(y_positions):
            tau = (y_plane - y0) / b
            position_plane = (pos0 + slope * tau)
            widths[i] = beam_fwhm_fixed_bin(position_plane,weights,bin_width=bin_width)
        finite = np.isfinite(widths)
        if not np.any(finite):
            return np.nan, np.nan
        valid_indices = np.flatnonzero(finite)
        i_min = valid_indices[np.argmin(widths[finite])]
        return (float(y_positions[i_min]),float(widths[i_min]))
        
    def beam_fwhm_fixed_bin(values,weights,bin_width=0.002):
        """
        Intensity-weighted FWHM using a fixed spatial bin width.
        The histogram expands automatically to contain the whole beam,
        but the resolution stays constant at every caustic position.
        For separated peaks, the FWHM spans all bins above half maximum.
        """
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
        xmin = np.min(values)
        xmax = np.max(values)
        if not (np.isfinite(xmin) and np.isfinite(xmax)):
            return np.nan
        if xmax <= xmin:
            return float(bin_width)
        # Fixed absolute histogram grid
        left = (np.floor(xmin / bin_width) * bin_width)
        right = (np.ceil(xmax / bin_width) * bin_width)
        edges = np.arange(left,right + 1.0001 * bin_width,bin_width)
        if edges.size < 3:
            return float(bin_width)
        hist, edges = np.histogram(values,bins=edges,weights=weights)
        if not np.any(hist > 0.0):
            return np.nan
        centers = 0.5 * (edges[:-1] + edges[1:])
        half = 0.5 * np.max(hist)
        above = np.flatnonzero(hist >= half)
        if above.size == 0:
            return np.nan
        if above.size == 1:
            return float(bin_width)
        return float(centers[above[-1]] - centers[above[0]])

    def vertical_focus_for_M1R(R):

        beamLine.M1.R = float(R)
        trial_input = rs.Beam(copyFrom=beamC2Global)
        beamM1Trial, _ = beamLine.M1.reflect(beam=trial_input)
        positions = np.arange(d_exp - 5000.0,d_exp + 2000.0 + 5.0,5.0)
        focus, width = directional_caustic(
            beam=beamM1Trial,
            y_positions=positions,
            direction="vertical",
            bin_width=0.002
        )
        return focus, width
    
    # ================================================================
    # Trace all fixed upstream elements once
    # ================================================================

    beamSource = beamLine.Source.shine()
    beamMask1 = beamLine.Mask1.propagate(beamSource)
    beamSCR1 = beamLine.Screen1.expose(beamMask1)
    beamFilter1global, beamFilter1local1, beamFilter1local2 = beamLine.Filter1.double_refract(beamSCR1)
    beamFilter2global, beamFilter2local1, beamFilter2local2 = beamLine.Filter2.double_refract(beamFilter1global)
    I0 = beam_intensity(beamSCR1)
    I_D = beam_intensity(beamFilter1global)
    I_SiC = beam_intensity(beamFilter2global)
    print("\nFILTER TRANSMISSION")
    print(f"before filters : {I0:.8e}")
    print(f"after diamond  : {I_D:.8e}   " f"T_D={I_D/I0:.6e}")
    print(f"after SiC      : {I_SiC:.8e}   " f"T_SiC={I_SiC/I_D:.6e}")
    print(f"total filters  : " f"T={I_SiC/I0:.6e}")

    # Fixed input immediately before C1.
    full_input_before_C1 = rs.Beam(copyFrom=beamFilter2global)

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

    finally:
        material_C1.method = requested_method_C1
        material_C2.method = requested_method_C2

    best_Fcurv = float(correction.Fcurv)

    def horizontal_focus_for_Fcurv(Fcurv):
        beamLine.LaueC1.Fcurv = float(Fcurv)
        beamLine.LaueC2.Fcurv = float(Fcurv)
        trial_input = rs.Beam(copyFrom=reduced_input_before_C1)
        beamC1Trial, _ = beamLine.LaueC1.reflect(beam=trial_input)
        beamC2Trial, _ = beamLine.LaueC2.reflect(beam=beamC1Trial)
        # We know focus should be somewhere around d_exp.
        # Broad enough for final correction.
        positions = np.arange(d_exp - 2500.0,d_exp + 2500.0 + 5.0,5.0)
        focus, width = directional_caustic(
            beam=beamC2Trial,
            y_positions=positions,
            direction="horizontal",
            bin_width=0.002
        )
        return focus, width

    F_candidates = best_Fcurv + np.array([
        -0.010,
        -0.005,
        0.000,
        0.005,
        0.010
    ])
    focus_results = []

    for F in F_candidates:

        focus, width = horizontal_focus_for_Fcurv(F)

        focus_results.append(
            (F, focus, width)
        )

        print(
            f"F={F:.6f}  "
            f"focus={focus:.1f} mm  "
            f"offset={focus-d_exp:.1f} mm  "
            f"FWHM={width:.4f} mm"
        )

    best = min(
        focus_results,
        key=lambda x: abs(x[1] - d_exp)
    )

    best_Fcurv = best[0]

    best_Fcurv = float(best[0])
    best_focus = float(best[1])
    best_width = float(best[2])

    beamLine.LaueC1.Fcurv = best_Fcurv
    beamLine.LaueC2.Fcurv = best_Fcurv

    # IMPORTANT: store the FINAL selected value for next energy
    beamLine._previous_Fcurv = best_Fcurv

    print(
        "\nFinal horizontal Fcurv:",
        f"{best_Fcurv:.9f}"
    )

    print(
        "Reduced-beam horizontal focus:",
        f"{best_focus:.3f} mm"
    )

    print(
        "Reduced-beam focus residual:",
        f"{best_focus - d_exp:.3f} mm"
    )

    print(
        "Horizontal FWHM at waist:",
        f"{best_width:.6e} mm"
    )

   
    print("Optimized double-Laue correction:",f"Fcurv={best_Fcurv:.9f}",)

    beamC1Input = rs.Beam(copyFrom=full_input_before_C1)
    beamC1Global, beamC1Local = (beamLine.LaueC1.reflect(beam=beamC1Input))

    beamC2Input = rs.Beam(copyFrom=beamC1Global)
    beamC2Global, beamC2Local = (beamLine.LaueC2.reflect(beam=beamC2Input))

    R0 = float(beamLine.M1.R)
    R_factors = np.arange(0.95,1.251,0.05)
    R_results = []
    for factor in R_factors:
        R = R0 * factor
        focus, width = vertical_focus_for_M1R(R)
        R_results.append((R, focus, width))
        print(
            f"R/R0={factor:.3f}  "
            f"R={R:.3f} mm  "
            f"focus={focus:.1f} mm  "
            f"offset={focus-d_exp:.1f} mm  "
            f"FWHM={width:.4f} mm"
        )
    R_focus, _, _ = min(R_results,key=lambda result:abs(result[1] - d_exp))

    beamLine.M1.R = float(R_focus)
    T1 = print_crystal_transmission("C1", beamC1Input, beamC1Global)
    T2 = print_crystal_transmission("C2", beamC2Input, beamC2Global)

    print(f"Double transmission = {T1*T2:.4f}")
    beamM1Global, beamM1Local = (beamLine.M1.reflect(beam=beamC2Global))
    beamSCRexp = beamLine.SCRexp.expose(beamM1Global)
    # ================================================================
    # Final beam widths directly from rays
    # ================================================================
    state = np.asarray(beamSCRexp.state)
    weights = np.real(np.asarray(beamSCRexp.Jss) + np.asarray(beamSCRexp.Jpp))
    good = ((state > 0) & np.isfinite(weights) & (weights > 0.0))
    x = np.asarray(beamSCRexp.x,dtype=float)[good]
    z = np.asarray(beamSCRexp.z,dtype=float)[good]
    a = np.asarray(beamSCRexp.a,dtype=float)[good]
    b = np.asarray(beamSCRexp.b,dtype=float)[good]
    c = np.asarray(beamSCRexp.c,dtype=float)[good]
    w = weights[good]
    # ------------------------------------------------
    # Position FWHM [mm]
    # ------------------------------------------------
    beamLine.final_beam_x, _, _, _ = beam_width(x,weights=w,statistic="FWHM",bin_width=0.005)
    beamLine.final_beam_z, _, _, _ = beam_width(z,weights=w,statistic="FWHM",bin_width=0.01)
    # ------------------------------------------------
    # Angular FWHM [mrad]
    # Use actual slopes relative to global y.
    # ------------------------------------------------
    xprime = (a / b) * 1e3
    zprime = (c / b) * 1e3
    beamLine.final_beam_dx, _, _, _ = beam_width(xprime,weights=w,statistic="FWHM",bin_width=0.005)
    beamLine.final_beam_dz, _, _, _ = beam_width(zprime,weights=w,statistic="FWHM",bin_width=0.002)
    print("\nFINAL BEAM WIDTHS FROM RAYS")
    print("x FWHM :", f"{beamLine.final_beam_x:.6e} mm")
    print("z FWHM :", f"{beamLine.final_beam_z:.6e} mm")
    print("x' FWHM:", f"{beamLine.final_beam_dx:.6e} mrad")
    print("z' FWHM:", f"{beamLine.final_beam_dz:.6e} mrad")
    
    # ================================================================
    # CAUSTICS
    # Horizontal: beam after C2 -> sagittal double-Laue focus
    # Vertical: beam after M1 -> meridional mirror focus
    # Broad scan: first energy only
    # Fine scan: every energy
    # ================================================================
    nominal_focus_y = float(beamLine.SCRexp.center[1])
    # ================================================================
    # HORIZONTAL CAUSTIC -- DOUBLE LAUE
    # ================================================================
    if not hasattr(beamLine,"_coarse_horizontal_focus"):
        print("\nCalculating broad horizontal caustic...")
        coarse_positions_x = np.arange(nominal_focus_y - 12000.0,nominal_focus_y + 3000.0 + 1.0,100.0)
        (coarse_horizontal_focus,coarse_horizontal_width) = directional_caustic(
            beam=beamC2Global,
            y_positions=coarse_positions_x,
            direction="horizontal",
            bin_width=0.002
        )
        beamLine._coarse_horizontal_focus = (coarse_horizontal_focus)
        print(
            "Broad horizontal focus:",
            f"{coarse_horizontal_focus:.3f} mm",
            "FWHM:",
            f"{coarse_horizontal_width:.6e} mm"
        )
    if hasattr(beamLine,"_previous_horizontal_focus"):
        fine_center_x = float(beamLine._previous_horizontal_focus)
    else:
        fine_center_x = float(beamLine._coarse_horizontal_focus)
    fine_positions_x = np.arange(fine_center_x - 200.0,fine_center_x + 200.0 + 1.0,1.0)
    (horizontal_focus,horizontal_width) = directional_caustic(
        beam=beamC2Global,
        y_positions=fine_positions_x,
        direction="horizontal",
        bin_width=0.001
    )
    edge_margin = 10.0
    if (horizontal_focus <= fine_positions_x[0] + edge_margin or horizontal_focus >= fine_positions_x[-1] - edge_margin):
        print("Horizontal focus near fine-scan boundary; ""expanding caustic.")
        expanded_positions_x = np.arange(fine_center_x - 1000.0,fine_center_x + 1000.0 + 2.0,2.0)
        (horizontal_focus,horizontal_width) = directional_caustic(
            beam=beamC2Global,
            y_positions=expanded_positions_x,
            direction="horizontal",
            bin_width=0.001
        )
    beamLine._previous_horizontal_focus = (horizontal_focus)

    # ================================================================
    # VERTICAL CAUSTIC -- AFTER M1
    # ================================================================
    if not hasattr(beamLine,"_coarse_vertical_focus"):
        print("\nCalculating broad vertical caustic...")
        coarse_positions_z = np.arange(d_exp - 8000.0,d_exp + 2000.0 + 1.0,50.0)
        (coarse_vertical_focus,coarse_vertical_width) = directional_caustic(
            beam=beamM1Global,
            y_positions=coarse_positions_z,
            direction="vertical",
            bin_width=0.002
        )
        beamLine._coarse_vertical_focus = (coarse_vertical_focus)
        print(
            "Broad vertical focus:",
            f"{coarse_vertical_focus:.3f} mm",
            "offset:",
            f"{coarse_vertical_focus - d_exp:.3f} mm",
            "FWHM:",
            f"{coarse_vertical_width:.6e} mm"
        )
    if hasattr(beamLine,"_previous_vertical_focus"):
        fine_center_z = float(beamLine._previous_vertical_focus)
    else:
        fine_center_z = float(beamLine._coarse_vertical_focus)
    fine_positions_z = np.arange(fine_center_z - 300.0,fine_center_z + 300.0 + 1.0,2.0)
    (vertical_focus,vertical_width) = directional_caustic(
        beam=beamM1Global,
        y_positions=fine_positions_z,
        direction="vertical",
        bin_width=0.001
    )
    edge_margin_z = 20.0
    if (vertical_focus <= fine_positions_z[0] + edge_margin_z or vertical_focus >= fine_positions_z[-1] - edge_margin_z):
        print("Vertical focus near fine-scan boundary; " "expanding caustic.")
        expanded_positions_z = np.arange(fine_center_z - 1500.0,fine_center_z + 1500.0 + 2.0,5.0)
        (vertical_focus,vertical_width) = directional_caustic(
            beam=beamM1Global,
            y_positions=expanded_positions_z,
            direction="vertical",
            bin_width=0.002
        )
    beamLine._previous_vertical_focus = (vertical_focus)

    print("\nCAUSTIC RESULTS")
    print(
        "Horizontal focus:",
        f"{horizontal_focus:.6f} mm"
    )
    print(
        "Horizontal focus offset:",
        f"{horizontal_focus - d_exp:.3f} mm"
    )
    print(
        "Minimum horizontal FWHM:",
        f"{horizontal_width:.6e} mm"
    )
    print(
        "Vertical focus:",
        f"{vertical_focus:.6f} mm"
    )
    print(
        "Vertical focus offset:",
        f"{vertical_focus - d_exp:.3f} mm"
    )
    print(
        "Minimum vertical FWHM:",
        f"{vertical_width:.6e} mm"
    )
    # ================================================================
    # Return beams required by define_plots()
    # ================================================================

    return {
        "beamSource": beamSource,
        "beamMask1": beamMask1,
        "beamSCR1": beamSCR1,
        "beamFilter1global": beamFilter1global,
        "beamFilter2global": beamFilter2global,
        "beamFilter1local1": beamFilter1local1,
        "beamFilter2local1": beamFilter2local1,
        "beamC1Global": beamC1Global,
        "beamC1Local": beamC1Local,
        "beamC2Global": beamC2Global,
        "beamC2Local": beamC2Local,
        "beamM1Global": beamM1Global,
        "beamM1Local": beamM1Local,
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

    plot_C2div = xrtp.XYCPlot('beamC2Global', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f'),
        fluxKind='total',
        title='LaueC2')
    plots.append(plot_C2div)

    plot_M1div = xrtp.XYCPlot('beamM1Global', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.3f'),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f'),
        fluxKind='total',
        title='M1')
    plots.append(plot_M1div)

    # after at focus
    plot_bs = xrtp.XYCPlot('beamSCRexp', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x$", r"mm", limits=None, fwhmFormatStr='%.3f'),
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

def update_scan_plot_axes(plots, M1pitch):
    """
    Reset or explicitly update spatial and angular plot limits
    before each scan point.

    Axis limits are expressed in the displayed plot units:
        position: mm
        divergence: mrad
    """

    # plots[0]: source divergence
    plots[0].xaxis.limits = None
    plots[0].yaxis.limits = None

    # plots[1]: divergence after Laue C2
    plots[1].xaxis.limits = None
    plots[1].yaxis.limits = None

    # plots[2]: divergence after M1
    plots[2].xaxis.limits = None

    # plots[3]: beam size at experimental screen
    plots[3].xaxis.limits = None
    plots[3].yaxis.limits = None

    # plots[4]: divergence at experimental screen
    plots[4].xaxis.limits = None

    # Nominal vertical angle after reflection by M1.
    zprime_center_mrad = 2.0 * M1pitch * 1e3

    # Use explicit moving limits for both post-M1 divergence plots.
    # Increase 0.20 if the distribution is wider.
    zprime_half_range_mrad = 0.20

    dynamic_zprime_limits = [
        zprime_center_mrad - zprime_half_range_mrad,
        zprime_center_mrad + zprime_half_range_mrad,
    ]

    plots[2].yaxis.limits = dynamic_zprime_limits
    plots[4].yaxis.limits = dynamic_zprime_limits

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
    fOut.write("E_eV\tflux0\tflux\tflux_density\tbeam_x\tbeam_z\tbeam_dx\tbeam_dz\tbeam_dE\tRs\tRm\n")
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
                                    Rs=-Rs, Rm=Rm, method='ML',ML_sign=1.0,
                                    PP_center_fraction=0.5, ST=ST)
        Cryst2 = BC2D.bentSi_Crys(hkl=hkl,geom='Laue Reflected',t=c_thick, aLrad=-aLrad,
                                    volumetricDiffraction=False,useTT=False,
                                    Rs=-Rs, Rm=Rm, method='ML',ML_sign=-1.0,
                                    PP_center_fraction=0.5, ST=ST)

        pitch1 = np.pi/2 + aLrad + thetaB
        pitch2 = np.pi/2 - aLrad + thetaB
    
        Fcurv = 1.0

        Earray = [70000, 60000, 50000, 40000, 30000]
        angle = [1.14e-3, 1.32e-3, 1.52e-3, 1.75e-3, 2.00e-3]
        coeff = np.polyfit(Earray, angle, deg=3)
        angE = np.poly1d(coeff)
        M1pitch = angE(e)
        print('pitch', M1pitch)
        # focusing
        M1R = 2.0 / (1.0 / d_m1 + 1 / (d_exp - d_m1)) / np.sin(M1pitch)
        # collimating
        #M1R = 2.0 * d_m1 / np.sin(M1pitch)
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

        beamLine.M1.R = M1R
        beamLine.M1.pitch = M1pitch
        beamLine.M1.center = [0, d_m1, fixedexit]

        yafterM1 = (d_exp - d_m1) * 2.0 * M1pitch
        print('yafterM1', yafterM1)

        beamLine.SCRexp.center = [0., d_exp, fixedexit + yafterM1]

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

        update_scan_plot_axes(
            plots=plots,
            M1pitch=M1pitch,
        )
                  
        for plot in plots[-2:]:
            plot.saveName = (
                f"{plot.title}_{e / 1000.0:.1f}keV.png"
            )
        
        yield
        # now all plots for this scan point are ready
        psrc = plots[0]
        ps = plots[3]
        pd = plots[4]
        sx_fwhm = beamLine.final_beam_x
        sz_fwhm = beamLine.final_beam_z
        dx_fwhm = beamLine.final_beam_dx
        dz_fwhm = beamLine.final_beam_dz
        flux0[ie] = psrc.flux
        flux[ie] = ps.flux
        flux_dens[ie] = ps.flux / ps.dx / ps.dy
        beam_x[ie] = sx_fwhm
        beam_z[ie] = sz_fwhm
        beam_dx[ie] = dx_fwhm
        beam_dz[ie] = dz_fwhm
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