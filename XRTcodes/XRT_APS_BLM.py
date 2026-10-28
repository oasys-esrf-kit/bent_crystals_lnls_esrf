

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

from dataclasses import dataclass
@dataclass
class CausticResult:
    longitudinal_positions: np.ndarray
    horizontal_width: np.ndarray
    vertical_width: np.ndarray
    horizontal_center: np.ndarray
    vertical_center: np.ndarray
    horizontal_left: np.ndarray
    horizontal_right: np.ndarray
    vertical_bottom: np.ndarray
    vertical_top: np.ndarray
    horizontal_focus: float
    vertical_focus: float
    minimum_horizontal_width: float
    minimum_vertical_width: float
    statistic: str

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

'''
def width_from_hist(axis, qlo=0.025, qhi=0.975):
    x = axis.binCenters
    w = axis.total1D
    good = np.isfinite(x) & np.isfinite(w) & (w > 0)
    x = x[good]
    w = w[good]

    if len(w) == 0:
        return np.nan, np.nan, np.nan

    cw = np.cumsum(w)
    cw /= cw[-1]
    xlo = np.interp(qlo, cw, x)
    xhi = np.interp(qhi, cw, x)
    return xhi - xlo, xlo, xhi
'''

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

def beam_width(values,weights=None,statistic="FWHM",bins=300,qlo=0.16,qhi=0.84):
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
        histogram, edges = np.histogram(values,bins=bins,weights=weights)
        histogram = np.asarray(histogram,dtype=float)
        centers = 0.5 * (edges[:-1] + edges[1:])
        if (histogram.size < 2 or not np.any(histogram > 0.0)):
            return np.nan, np.nan, np.nan, center
        peak_index = int(np.argmax(histogram))
        half_maximum = 0.5 * histogram[peak_index]
        # Search only the connected peak containing the global maximum.
        # This avoids spanning two separated peaks in a multimodal beam.
        left_index = peak_index
        while (left_index > 0 and histogram[left_index - 1] >= half_maximum):
            left_index -= 1
        right_index = peak_index
        while (right_index < histogram.size - 1 and histogram[right_index + 1] >= half_maximum):
            right_index += 1
        # Left half-maximum crossing.
        if left_index == 0:
            left = edges[0]
        else:
            x1 = centers[left_index - 1]
            x2 = centers[left_index]
            y1 = histogram[left_index - 1]
            y2 = histogram[left_index]
            if y2 == y1:
                left = x2
            else:
                left = x1 + (half_maximum - y1) * (x2 - x1) / (y2 - y1)
        # Right half-maximum crossing.
        if right_index == histogram.size - 1:
            right = edges[-1]
        else:
            x1 = centers[right_index]
            x2 = centers[right_index + 1]
            y1 = histogram[right_index]
            y2 = histogram[right_index + 1]

            if y2 == y1:
                right = x1
            else:
                right = x1 + (half_maximum - y1) * (x2 - x1) / (y2 - y1)
        return (right - left,left,right,center)

    raise ValueError("Unknown statistic. Use 'FWHM', 'RMS' or '95%'.")


def beam_caustic(beam,longitudinal_positions,statistic="FWHM",bins=300,qlo=0.025,qhi=0.975,good_only=True,minimum_rays=20):
    """
    Calculate horizontal and vertical beam caustics.
    The input beam is propagated analytically to a sequence of planes
    perpendicular to the global y direction.
    Parameters
    ----------
    beam : xrt.backends.raycing.sources.Beam
        Beam after the last optical element.
        Required arrays:
            beam.x, beam.y, beam.z
            beam.a, beam.b, beam.c
        Statistical weights are calculated from:
            beam.Jss + beam.Jpp
    longitudinal_positions : array-like
        Global y positions of the caustic planes, in mm.
    statistic : {"FWHM", "RMS", "95%"}
        Width definition used independently for x and z.
    bins : int
        Histogram bins used for FWHM.
    qlo, qhi : float
        Quantiles used for statistic="95%".
    good_only : bool
        If True, use only rays with beam.state > 0.
    minimum_rays : int
        Minimum number of valid rays required at each plane.
    Returns
    -------
    CausticResult
        Horizontal and vertical width curves and their minima.
    Notes
    -----
    The propagation is:
        tau = (y_plane - y_ray) / b_ray
        x_plane = x_ray + a_ray * tau
        z_plane = z_ray + c_ray * tau
    Therefore the scan planes are perpendicular to the global y axis.
    """
    y_scan = np.asarray(longitudinal_positions,dtype=float)
    if y_scan.ndim != 1 or y_scan.size < 2:
        raise ValueError(
            "longitudinal_positions must be a 1D array "
            "with at least two positions."
        )
    
    x0 = np.asarray(beam.x, dtype=float)
    y0 = np.asarray(beam.y, dtype=float)
    z0 = np.asarray(beam.z, dtype=float)
    a = np.asarray(beam.a, dtype=float)
    b = np.asarray(beam.b, dtype=float)
    c = np.asarray(beam.c, dtype=float)
    number_of_rays = x0.size

    for name, array in (
        ("beam.y", y0),
        ("beam.z", z0),
        ("beam.a", a),
        ("beam.b", b),
        ("beam.c", c),
    ):
        if array.size != number_of_rays:
            raise ValueError(f"{name} has a different number of rays.")

    if hasattr(beam, "Jss") and hasattr(beam, "Jpp"):
        weights = np.real(np.asarray(beam.Jss) + np.asarray(beam.Jpp))
    else:
        weights = np.ones(number_of_rays,dtype=float)

    valid = (
        np.isfinite(x0)
        & np.isfinite(y0)
        & np.isfinite(z0)
        & np.isfinite(a)
        & np.isfinite(b)
        & np.isfinite(c)
        & np.isfinite(weights)
        & (weights > 0.0)
        & (np.abs(b) > 1e-15)
    )
    if good_only and hasattr(beam, "state"):
        valid &= (np.asarray(beam.state) > 0)

    x0 = x0[valid]
    y0 = y0[valid]
    z0 = z0[valid]
    a = a[valid]
    b = b[valid]
    c = c[valid]
    weights = weights[valid]

    if x0.size < minimum_rays:
        raise ValueError(
            f"Only {x0.size} valid rays are available; "
            f"minimum_rays={minimum_rays}."
        )

    number_of_planes = y_scan.size
    horizontal_width = np.full(number_of_planes,np.nan)
    vertical_width = np.full(number_of_planes,np.nan)
    horizontal_center = np.full(number_of_planes,np.nan)
    vertical_center = np.full(number_of_planes,np.nan)
    horizontal_left = np.full(number_of_planes,np.nan)
    horizontal_right = np.full(number_of_planes,np.nan)
    vertical_bottom = np.full(number_of_planes,np.nan)
    vertical_top = np.full(number_of_planes,np.nan)

    for index, y_plane in enumerate(y_scan):
        tau = (y_plane - y0) / b
        x_plane = x0 + a * tau
        z_plane = z0 + c * tau
        (horizontal_width[index],horizontal_left[index],horizontal_right[index],horizontal_center[index],
            ) = beam_width(values=x_plane,weights=weights,statistic=statistic,bins=bins,qlo=qlo,qhi=qhi)
        (vertical_width[index],vertical_bottom[index],vertical_top[index],vertical_center[index],
            ) = beam_width(values=z_plane,weights=weights,statistic=statistic,bins=bins,qlo=qlo,qhi=qhi)
    finite_horizontal = np.isfinite(horizontal_width)
    finite_vertical = np.isfinite(vertical_width)

    if np.any(finite_horizontal):
        horizontal_indices = np.flatnonzero(finite_horizontal)
        local_index = np.argmin(horizontal_width[finite_horizontal])
        horizontal_focus_index = (horizontal_indices[local_index])
        horizontal_focus = y_scan[horizontal_focus_index]
        minimum_horizontal_width = (horizontal_width[horizontal_focus_index])
    else:
        horizontal_focus = np.nan
        minimum_horizontal_width = np.nan

    if np.any(finite_vertical):
        vertical_indices = np.flatnonzero(finite_vertical)
        local_index = np.argmin(vertical_width[finite_vertical])
        vertical_focus_index = (vertical_indices[local_index])
        vertical_focus = y_scan[vertical_focus_index]
        minimum_vertical_width = (vertical_width[vertical_focus_index])
    else:
        vertical_focus = np.nan
        minimum_vertical_width = np.nan

    return CausticResult(
        longitudinal_positions=y_scan,
        horizontal_width=horizontal_width,
        vertical_width=vertical_width,
        horizontal_center=horizontal_center,
        vertical_center=vertical_center,
        horizontal_left=horizontal_left,
        horizontal_right=horizontal_right,
        vertical_bottom=vertical_bottom,
        vertical_top=vertical_top,
        horizontal_focus=horizontal_focus,
        vertical_focus=vertical_focus,
        minimum_horizontal_width=(minimum_horizontal_width),
        minimum_vertical_width=(minimum_vertical_width),
        statistic=statistic.upper(),
    )


# ==========================================================================
# beamline parameters

# source parameters

bener = 7. #ring energy [GeV]
# e-beam size and div RMS
esizeh = 91.
esizev = 30.
edivh = 70.
edivv = 43.
# BM
B0 = 0.6 # [T]

delta_E = 0.005

# beam acceptance
acch = 1.0 #~[mrad]
accv = 0.5 # [mrad]

# distances
d_slit1 = 10.0 *1e3 # [mm] distance slits and filters
d_bragg = 26.5 *1e3 # [mm] distance source - DBLM C1
d_laue = 26.7 *1e3 # [mm] distance source - DBLM C2
d_exp = 52.0 * 1e3 # [mm] distance source - exp
fixedexit = -50. # [mm]

# crystals
c_thick_laue = 0.6 #[mm]
crystal = 'Si' # crystal Si(111) on (100) substrate 
hkl = [3,1,1]
nrm = [5,-1,-1]
sag = [0,1,-1]
h, k, l = hkl
str_hkl = 'Si' + "".join(map(str, hkl))
str_nrm = "".join(map(str, nrm))
str_sag = "".join(map(str, sag))

aB = 0.0 # symmetric Bragg
aLcalc = angle_vec(hkl, nrm)
aLrad = np.pi/2. - aLcalc
print('aL',np.rad2deg(aLrad))

# standard length unit [mm] - unitfac=1e3
unitfac = 1e3
# mechanical properties of Laue crystals
c_fac = 0.7
CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7,63.9,79.6], hkl=nrm, sag=sag, unitfac=unitfac)
poL = poisson[1,2]  # ST'(ij)= ST(ij) / ST(33) - poL=ST'(23)
Sf = 'S_tensor.txt'
np.savetxt(Sf, ST)

# energy, scan parameters
energies = np.linspace(30000, 60000, 16)

# simulation parameters
repeats = 1
nrays = 500000

title0 = 'APS_BLM_' + str_hkl + '_nrm' + str_nrm + '_sag' + str_sag +  '_ML_'
# ==================================================================
def build_beamline():
    height = 0.0
   
    beamLine = raycing.BeamLine(azimuth=0.0, height=height, alignE='auto')
   
    # wiggler source
    beamLine.source = rs.BendingMagnet(beamLine, name='DW', center=(0,0,0),
        nrays=nrays, B0=B0, eE=bener, eI=0.1,
        eSigmaX=esizeh, eSigmaZ=esizev, eEpsilonX=esizeh*edivh*1e-3, eEpsilonZ=esizev*edivv*1e-3,
        eMin=energies[0]*(1-delta_E),eMax=energies[0]*(1+delta_E),
        xPrimeMax=acch, zPrimeMax=accv)
    # ----------------------------------------------------------------------------------
    # aperture and filters @ 10m
    s1_xwidth = d_slit1 * acch * 1e-3
    s1_zwidth = d_slit1 * accv * 1e-3
    filterDiamond = rm.Material('C', rho=3.52, kind='plate')
    filterSiC = rmc.SiliconCarbide()
    #filterSiC = rm.Material(filterSiCm, rho=3.52, kind='plate')
    beamLine.mask1 = ra.RectangularAperture(beamLine, 'Mask1', [0., d_slit1, height],
        ('left', 'right', 'bottom', 'top'), [-s1_xwidth/2., s1_xwidth/2., -s1_zwidth/2., s1_zwidth/2.])
    beamLine.screen1 = rsc.Screen(beamLine, 'SCR1', [0., d_slit1, height])

    beamLine.Bragg = roe.JohannCylinder(
        name='Bragg',
        bl=beamLine,
        alpha=aB,
        limPhysX=[-50.,50.],
        limPhysY=[-500.,500.]
    )
    beamLine.Laue = BC2D.BentCrys2D(
        name='Laue',
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
    Trace source -> Bragg once.
    Use a fixed deterministic subset of the Bragg beam to optimize
    the Laue sagittal-curvature correction Fcurv.
    Finally trace the full Bragg beam through the Laue crystal.
    """
    # ==============================================================
    # Optimization settings
    # ==============================================================
    optimization_ray_count = 50000
    optimization_seed = 12345
    optimization_method = "ML"

    # ==============================================================
    # Source -> mask -> Bragg
    # ==============================================================
    beamSource = (beamLine.source.shine())
    beamMask1 = (beamLine.mask1.propagate(beamSource))
    beamSCR1 = (beamLine.screen1.expose(beamMask1))
    beamBraggGlobal, beamBraggLocal = (beamLine.Bragg.reflect(beam=rs.Beam(copyFrom=beamMask1)))
    # ==============================================================
    # Function to copy a beam subset
    # ==============================================================
    def copy_beam_subset(input_beam,indices):
        output_beam = rs.Beam(copyFrom=input_beam)
        number_of_rays = len(input_beam.x)
        for (attribute_name,attribute_value) in vars(output_beam).items():
            if not isinstance(attribute_value,np.ndarray):
                continue
            if (attribute_value.ndim == 0):
                continue
            if (attribute_value.shape[0] != number_of_rays):
                continue
            setattr(output_beam,attribute_name,np.array(attribute_value[indices],copy=True))

        return output_beam

    # ==============================================================
    # Fixed optimization subset
    # ==============================================================
    if hasattr(beamBraggGlobal,"state"):
        candidate_indices = (np.flatnonzero(beamBraggGlobal.state > 0))
    else:
        candidate_indices = (np.arange(len(beamBraggGlobal.x)))
    if candidate_indices.size == 0:
        raise RuntimeError("No rays survived the Bragg crystal.")

    number_for_optimization = min(optimization_ray_count,candidate_indices.size)
    rng = np.random.default_rng(optimization_seed)
    # deterministic permutation
    permutation = rng.permutation(candidate_indices)
    optimization_indices = np.sort(permutation[:number_for_optimization])
    beamBraggOptimization = (copy_beam_subset(beamBraggGlobal,optimization_indices))

    def trace_Fcurv_trial(Fcurv):
        beamLine.Laue.Fcurv = (float(Fcurv))
        trial_input = rs.Beam(copyFrom=beamBraggOptimization)
        (beamTrial, _) = beamLine.Laue.reflect(beam=trial_input)
        return (beamTrial,{"Fcurv":float(Fcurv)})

    # ==============================================================
    # Laue-only Fcurv trial
    # ==============================================================
    laue_material = (beamLine.Laue.material[0])
    final_laue_method = (laue_material.method)

    # ==============================================================
    # Fcurv optimization
    # First energy:
    #     broad optimization -> fine optimization
    # Following energies:
    #     previous optimum -> fine optimization only
    # ==============================================================

    try:
        laue_material.method = optimization_method
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
        laue_material.method = (final_laue_method)

    beamLine.Laue.Fcurv = (best_Fcurv)
    print(
        "\nOptimized geometry correction:",
        f"Fcurv={best_Fcurv:.9f}",
        f"(optimization method={optimization_method}, "
        f"final method={final_laue_method})"
    )

    # ==============================================================
    # Final FULL Laue trace
    # ==============================================================
    (beamLaueGlobal,beamLaueLocal) = beamLine.Laue.reflect(
        beam=rs.Beam(copyFrom=beamBraggGlobal))
    beamSCRexp = (beamLine.SCRexp.expose(beamLaueGlobal))

    # ==============================================================
    # CAUSTIC
    # Broad search: only once
    # Fine search: every energy
    # ==============================================================
    nominal_focus_y = float(beamLine.SCRexp.center[1])
    # 1. BROAD CAUSTIC -- ONLY ONCE
    if not hasattr(beamLine, "_coarse_caustic_focus"):
        print("\nCalculating broad caustic...")
        coarse_positions = np.arange(nominal_focus_y - 12000.0,nominal_focus_y + 3000.0 + 1.0,100.0)
        coarse_caustic = beam_caustic(
            beam=beamLaueGlobal,
            longitudinal_positions=coarse_positions,
            statistic="FWHM",
            bins=200,
            good_only=True
        )
        beamLine._coarse_caustic_focus = float(coarse_caustic.horizontal_focus)
        print(
            "Broad horizontal focus:",
            f"{beamLine._coarse_caustic_focus:.3f} mm"
        )
    # 2. FINE CAUSTIC -- EVERY ENERGY
    if hasattr(beamLine, "_previous_fine_focus"):
        fine_center = (beamLine._previous_fine_focus)
    else:
        fine_center = (beamLine._coarse_caustic_focus)
    fine_positions = np.arange(fine_center - 200.0,fine_center + 200.0 + 1.0,1.0)
    fine_caustic = beam_caustic(
        beam=beamLaueGlobal,
        longitudinal_positions=fine_positions,
        statistic="FWHM",
        bins=300,
        good_only=True
    )
    horizontal_focus = float(fine_caustic.horizontal_focus)
    horizontal_width = float(fine_caustic.minimum_horizontal_width)
    vertical_focus = float(fine_caustic.vertical_focus)
    vertical_width = float(fine_caustic.minimum_vertical_width)
    # Store horizontal focus for next energy
    beamLine._previous_fine_focus = (horizontal_focus)
    print(
        "Horizontal focus:",
        f"{horizontal_focus:.6f} mm"
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
        "Minimum vertical FWHM:",
        f"{vertical_width:.6e} mm"
    )
    # ==============================================================
    # Return beams for XRT plots
    # ==============================================================

    return {
        "beamSource":beamSource,
        "beamMask1":beamMask1,
        "beamSCR1":beamSCR1,
        "beamBraggGlobal":beamBraggGlobal,
        "beamBraggLocal":beamBraggLocal,
        "beamLaueGlobal":beamLaueGlobal,
        "beamLaueLocal":beamLaueLocal,
        "beamSCRexp":beamSCRexp,
    }

rr.run_process = run_process

def define_plots(beamLine):
    plots = []
    bins = 200
    # at mask @ 10m - only aperture
    plot_scr = xrtp.XYCPlot('beamSCR1', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x$", r"mm", limits=None, fwhmFormatStr='%.2f',bins=bins),
        yaxis=xrtp.XYCAxis(r"$z$", r"mm", limits=None, fwhmFormatStr='%.2f',bins=bins),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.5f',bins=bins),
        fluxKind='total',
        title='SCR1')
    plots.append(plot_scr)
    # after Bragg
    plot_Bragg = xrtp.XYCPlot('beamBraggGlobal', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.3f',bins=bins),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.3f',bins=bins),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.5f',bins=bins),
        fluxKind='total',
        title='Bragg_div')
    plots.append(plot_Bragg)

     # after Laue
    plot_Laue = xrtp.XYCPlot('beamLaueGlobal', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.3f',bins=bins),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.3f',bins=bins),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f',bins=bins),
        fluxKind='total',
        title='Laue_div')
    plots.append(plot_Laue)

    # at focus
    plot_bs = xrtp.XYCPlot('beamSCRexp', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x$", r"mm", limits=None, fwhmFormatStr='%.3f',bins=bins),
        yaxis=xrtp.XYCAxis(r"$z$", r"mm", limits=[-5,5], fwhmFormatStr='%.3f',bins=bins),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f',bins=bins),
        fluxKind='total',
        title=title0+'SCRexp-bs')
    plots.append(plot_bs)
    plot_div = xrtp.XYCPlot('beamSCRexp', aspect='auto',
        xaxis=xrtp.XYCAxis(r"$x'$", r"mrad", limits=None, fwhmFormatStr='%.3f',bins=bins),
        yaxis=xrtp.XYCAxis(r"$z'$", r"mrad", limits=None, fwhmFormatStr='%.3f',bins=bins),
        caxis=xrtp.XYCAxis('energy', 'keV', fwhmFormatStr='%.3f',bins=bins),
        fluxKind='total',
        title=title0+'SCRexp-div')
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
    fOut.write("E_eV\tflux0\tflux\tflux_density\tbeam_x\tbeam_z\tbeam_dx\tbeam_dz\tdE\tRmB\tRsL\tRmL\n")
    for ie, e in enumerate(energies):
        print(f'energy {e:.1f} eV, {ie+1} of {len(energies)}')
        # thetaB and distnace between crystals
        # Bragg angle corrected for refraction (xrt)
        crystal_Bragg = rm.CrystalSi(hkl=hkl)
        thetaB = crystal_Bragg.get_Bragg_angle(e) - crystal_Bragg.get_dtheta_symmetric_Bragg(e)  # radians
        d_bragg_laue = abs(fixedexit) / np.tan(2. * thetaB)
        # Horizontal Focusing (Laue bending) - Rs, Rm
        pL = d_bragg + d_bragg_laue
        FL2 = d_exp - (d_bragg + d_bragg_laue) # Laue q
        qL = FL2
        fs = 1 / (1 / pL + 1 / qL)
        RsL = 2 * np.sin(thetaB) * np.sin(aLrad) * fs
        RmL = RsL / c_fac / poL
        FB2 = d_exp - d_bragg # Bragg q
        FB2v = - np.cos(np.abs(aLrad) + thetaB) / (-2. / RmL + np.cos(np.abs(aLrad) - thetaB) / FL2) + FB2 - FL2
        RmB = 2.0 / ((np.sin(aB + thetaB)/d_bragg) - (np.sin(aB - thetaB)/FB2v))

        # crystal material parameters
        CrystB = BC2D.bentSi_Crys(hkl=hkl, geom='Bragg reflected', t=10., aLrad=aB,
            volumetricDiffraction=False, useTT=False,
            Rs=1e100, Rm=RmB, method='ML')
        CrystL = BC2D.bentSi_Crys(hkl=hkl,geom='Laue Reflected',t=c_thick_laue, aLrad=np.pi/2. - aLrad, 
            volumetricDiffraction=False, useTT=False,
            Rs=RsL, Rm=-RmL, method='ML',ML_sign=-1.0,
            PP_center_fraction=0.5)

        pitch_Bragg = np.pi - aB - thetaB
        pitch_Laue = np.pi/2. - aLrad - thetaB
   
        Fcurv = 1.0
            
        beamLine.Bragg.material = (CrystB,)
        beamLine.Bragg.Rm = RmB
        beamLine.Bragg.pitch = pitch_Bragg
        beamLine.Bragg.center = [0, d_bragg, 0]

        beamLine.Laue.material = (CrystL,)
        beamLine.Laue.Rs = RsL
        beamLine.Laue.Rm = - RmL
        beamLine.Laue.pitch = pitch_Laue
        beamLine.Laue.Fcurv = Fcurv
        beamLine.Laue.center = [0, d_bragg + d_bragg_laue, fixedexit]

        beamLine.SCRexp.center = [0, d_bragg + d_bragg_laue + FL2, fixedexit]
        beamLine.alignE = e
        beamLine.source.eMin = e * (1 - delta_E)
        beamLine.source.eMax = e * (1 + delta_E)
        if hasattr(beamLine.source, "reset"):
            beamLine.source.reset()

        update_energy_axes(
            plots=plots,
            energy_eV=e,
            relative_bandwidth=delta_E,
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
        lowcut = 0.025
        highcut = 0.975
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
            f"{RmB:.6f}\t"
            f"{RsL:.6f}\t"
            f"{RmL:.6f}\n"
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