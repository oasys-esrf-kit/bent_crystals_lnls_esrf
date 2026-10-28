"""
Generic sagittal-curvature correction for XRT.
The optimizer changes Fcurv and minimizes the intensity-weighted
FWHM on a specified target plane.
Important:
    - FWHM uses a FIXED histogram range for every Fcurv trial.
    - If the beam is split into two peaks and both peaks are above
      half maximum, the width spans both peaks.
    - No parabolic assumption is made for FWHM.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any, Callable, Dict, Optional, Sequence, Tuple
import numpy as np
# =====================================================================
# Results
# =====================================================================
@dataclass
class GeometryTrial:
    Fcurv: float
    sigma: float       # contains FWHM; name retained for compatibility
    centroid: float
    n_rays: int
    final_beam: Any
    payload: Any = None

@dataclass
class GeometryCorrectionResult:
    Fcurv: float
    sigma: float       # contains FWHM
    centroid: float
    n_rays: int
    final_beam: Any
    payload: Any
    history: Sequence[GeometryTrial]

# =====================================================================
# Utilities
# =====================================================================
def _unit(vector, name):
    v = np.asarray(vector, dtype=float)
    if v.shape != (3,):
        raise ValueError(f"{name} must contain exactly three components.")
    norm = np.linalg.norm(v)
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError(f"{name} must be finite and non-zero.")
    return v / norm

# =====================================================================
# Infer propagation direction
# =====================================================================
def infer_beam_direction(beam,minimum_rays=200,use_intensity_weights=True):
    directions = np.column_stack(
        (
            np.asarray(beam.a, dtype=float),
            np.asarray(beam.b, dtype=float),
            np.asarray(beam.c, dtype=float),
        )
    )
    good = ((np.asarray(beam.state) > 0) & np.all(np.isfinite(directions),axis=1))
    if np.count_nonzero(good) < minimum_rays:
        raise RuntimeError("Too few rays to infer beam direction.")
    directions = directions[good]
    central_direction = None
    if (use_intensity_weights and hasattr(beam, "Jss") and hasattr(beam, "Jpp")):
        weights = np.real(
            np.asarray(beam.Jss, dtype=float)[good]
            + np.asarray(beam.Jpp, dtype=float)[good]
        )
        wgood = (np.isfinite(weights) & (weights > 0.0))
        if np.count_nonzero(wgood) >= minimum_rays:
            central_direction = np.average(directions[wgood],axis=0,weights=weights[wgood])
    if central_direction is None:
        central_direction = np.median(directions,axis=0)
    return _unit(central_direction,"inferred plane normal")

# =====================================================================
# FWHM
# =====================================================================
def _weighted_fwhm(values,weights,bin_width=0.002,):
    """
    Intensity-weighted FWHM with a FIXED spatial bin width.
    The histogram range automatically expands to contain the complete
    beam, so defocused rays are never discarded.
    For split beams the FWHM spans all histogram bins above half
    maximum, thereby penalizing a double focus.
    """
    values = np.asarray(values,dtype=float)
    weights = np.asarray(weights,dtype=float)
    good = (np.isfinite(values) & np.isfinite(weights) & (weights > 0.0))
    values = values[good]
    weights = weights[good]
    if values.size < 2:
        return np.nan
    # --------------------------------------------------------------
    # Fixed histogram grid
    # --------------------------------------------------------------
    xmin = np.min(values)
    xmax = np.max(values)
    if not (np.isfinite(xmin) and np.isfinite(xmax)):
        return np.nan
    if xmax <= xmin:
        return float(bin_width)

    # Align every histogram to the same absolute grid.
    # Example for bin_width = 0.002 mm:
    # ..., -0.004, -0.002, 0, 0.002, 0.004, ...
    left = (np.floor(xmin / bin_width) * bin_width)
    right = (np.ceil(xmax / bin_width) * bin_width)
    # guarantee at least two bins
    if right <= left:
        right = left + 2.0 * bin_width
    edges = np.arange(left,right + 1.0001 * bin_width,bin_width)
    if edges.size < 3:
        edges = np.array([left,left + bin_width,left + 2.0 * bin_width])
    histogram, edges = np.histogram(values,bins=edges,weights=weights)
    histogram = np.asarray(histogram,dtype=float)
    if (histogram.size < 2 or not np.any(histogram > 0.0)):
        return np.nan
    centers = 0.5 * (edges[:-1] + edges[1:])
    # --------------------------------------------------------------
    # FWHM
    # --------------------------------------------------------------
    half_maximum = (0.5 * np.max(histogram))
    above = (histogram >= half_maximum)
    indices = np.flatnonzero(above)
    if indices.size == 0:
        return np.nan
    if indices.size == 1:
        return float(bin_width)
    # Deliberately span separated peaks.
    width = (centers[indices[-1]] - centers[indices[0]])
    return float(width)

# =====================================================================
# Beam FWHM on arbitrary plane
# =====================================================================
def beam_fwhm_on_plane(
    beam,
    plane_point,
    plane_normal,
    transverse_axis=(1.0, 0.0, 0.0),
    use_intensity_weights=True,
    minimum_rays=100,
    require_forward_intersection=True,
    bin_width=0.002,
):
    p = np.asarray(plane_point,dtype=float)
    if p.shape != (3,):
        raise ValueError("plane_point must contain three components.")
    n = _unit(plane_normal,"plane_normal")
    e = _unit(transverse_axis,"transverse_axis")
    # transverse direction must lie in screen plane
    e = e - np.dot(e, n) * n
    en = np.linalg.norm(e)
    if (not np.isfinite(en)or en <= 1e-14):
        raise ValueError("transverse_axis is parallel to plane_normal.")
    e /= en
    # ---------------------------------------------------------------
    # ray positions and directions
    # ---------------------------------------------------------------
    r0 = np.column_stack((
        np.asarray(beam.x, dtype=float),
        np.asarray(beam.y, dtype=float),
        np.asarray(beam.z, dtype=float),))
    v = np.column_stack((
        np.asarray(beam.a, dtype=float),
        np.asarray(beam.b, dtype=float),
        np.asarray(beam.c, dtype=float),))
    state_good = (np.asarray(beam.state) > 0)
    finite = (np.all(np.isfinite(r0),axis=1) & np.all(np.isfinite(v),axis=1))
    denominator = v @ n
    valid = (state_good & finite & np.isfinite(denominator) & (np.abs(denominator) > 1e-14))
    # ---------------------------------------------------------------
    # intersect rays with target plane
    # ---------------------------------------------------------------
    t = np.full(len(state_good),np.nan,dtype=float)
    t[valid] = ((p - r0[valid]) @ n) / denominator[valid]
    if require_forward_intersection:
        valid &= (t >= 0.0)
    if (np.count_nonzero(valid) < minimum_rays):
        raise RuntimeError("Too few rays reach the target plane.")
    intersection = (r0[valid] + t[valid, None] * v[valid])
    coordinate = (intersection - p) @ e
    finite_coordinate = (np.isfinite(coordinate))
    coordinate = coordinate[finite_coordinate]
    # ---------------------------------------------------------------
    # intensity weights
    # ---------------------------------------------------------------
    if (use_intensity_weights and hasattr(beam, "Jss") and hasattr(beam, "Jpp")):
        weights = np.real(np.asarray(beam.Jss)[valid][finite_coordinate] + np.asarray(beam.Jpp)[valid][finite_coordinate])
        good_weights = (np.isfinite(weights) & (weights > 0.0))
        coordinate = coordinate[good_weights]
        weights = weights[good_weights]
    else:
        weights = np.ones_like(coordinate,dtype=float)
    if coordinate.size < minimum_rays:
        raise RuntimeError("Too few weighted rays remain.")
    # ---------------------------------------------------------------
    # centroid
    # ---------------------------------------------------------------
    weight_sum = np.sum(weights)
    centroid = float(np.sum(coordinate * weights) / weight_sum)
    # ---------------------------------------------------------------
    # FWHM
    # ---------------------------------------------------------------
    width = _weighted_fwhm(coordinate,weights,bin_width=bin_width)
    if not np.isfinite(width):
        raise RuntimeError("Could not determine FWHM.")
    return (width,centroid,int(coordinate.size))

# =====================================================================
# Fcurv optimizer
# =====================================================================
def optimize_geometry_correction(
    trace_trial: Callable[[float], Tuple[Any, Any]],
    plane_point,
    plane_normal=None,
    transverse_axis=(1.0, 0.0, 0.0),
    initial_Fcurv=1.0,
    initial_step=0.05,
    bounds=(0.5, 1.5),
    max_bracket_expansions=5,
    refinement_iterations=4,
    grid_points=5,
    Fcurv_tolerance=1e-4,
    bin_width=0.002,
    use_intensity_weights=True,
    minimum_rays=200,
    verbose=False,
):
    """
    Minimize transverse FWHM at the fixed target plane.
    Search:
        1. automatically bracket the Fcurv minimum
        2. repeatedly refine that bracket with a small regular grid
    No parabolic approximation is used.
    """
    fmin, fmax = map(float,bounds)
    if not (0.0 < fmin < fmax):
        raise ValueError("Invalid Fcurv bounds.")
    if initial_step <= 0.0:
        raise ValueError("initial_step must be > 0.")
    if grid_points < 5:
        raise ValueError("grid_points must be >= 5.")
    history = []
    cache: Dict[float,GeometryTrial] = {}
    # ---------------------------------------------------------------
    # target-plane normal
    # ---------------------------------------------------------------
    if plane_normal is None:
        resolved_plane_normal = None
    else:
        resolved_plane_normal = _unit(plane_normal,"plane_normal")
    # ---------------------------------------------------------------
    # evaluate one Fcurv
    # ---------------------------------------------------------------
    def evaluate(Fcurv):
        nonlocal resolved_plane_normal
        F = float(np.clip(Fcurv,fmin,fmax))
        key = round(F,12)
        if key in cache:
            return cache[key]
        final_beam, payload = (trace_trial(F))
        if resolved_plane_normal is None:
            resolved_plane_normal = (
                infer_beam_direction(
                    final_beam,
                    minimum_rays=minimum_rays,
                    use_intensity_weights=use_intensity_weights))
            if verbose:
                print("Inferred target-plane normal:",resolved_plane_normal)
        width, centroid, n_rays = (
            beam_fwhm_on_plane(
                final_beam,
                plane_point=plane_point,
                plane_normal=resolved_plane_normal,
                transverse_axis=transverse_axis,
                use_intensity_weights=use_intensity_weights,
                minimum_rays=minimum_rays,
                bin_width=bin_width))
        trial = GeometryTrial(
            Fcurv=F,
            sigma=width,
            centroid=centroid,
            n_rays=n_rays,
            final_beam=final_beam,
            payload=payload
        )
        cache[key] = trial
        history.append(trial)
        if verbose:
            print(
                f"Geometry trial: "
                f"Fcurv={F:.9f}, "
                f"FWHM={width:.9g} mm, "
                f"centroid={centroid:.9g} mm, "
                f"rays={n_rays}"
            )
        return trial

    # ===============================================================
    # STEP 1 -- bracket minimum
    # ===============================================================
    center = float(np.clip(initial_Fcurv,fmin,fmax))
    step = float(initial_step)
    left = evaluate(center - step)
    middle = evaluate(center)
    right = evaluate(center + step)
    for _ in range(max_bracket_expansions):
        if (middle.sigma <= left.sigma and middle.sigma <= right.sigma):
            break
        if left.sigma < middle.sigma:
            right = middle
            middle = left
            step *= 1.6
            candidate = (middle.Fcurv - step)
            if candidate <= fmin:
                left = evaluate(fmin)
                break
            left = evaluate(candidate)
        else:
            left = middle
            middle = right
            step *= 1.6
            candidate = (middle.Fcurv + step)
            if candidate >= fmax:
                right = evaluate(fmax)
                break
            right = evaluate(candidate)
    bracket = sorted((left,middle,right),key=lambda item:item.Fcurv)
    left, middle, right = bracket
    # ===============================================================
    # STEP 2 -- regular grid refinement
    # ===============================================================
    for iteration in range(refinement_iterations):
        lo = left.Fcurv
        hi = right.Fcurv
        if hi <= lo:
            break
        Fgrid = np.linspace(lo,hi,grid_points)
        trials = [evaluate(F) for F in Fgrid]
        ibest = min(range(len(trials)),key=lambda i:trials[i].sigma)
        if ibest == 0:
            left = trials[0]
            middle = trials[0]
            right = trials[1]
        elif ibest == len(trials) - 1:
            left = trials[-2]
            middle = trials[-1]
            right = trials[-1]
        else:
            left = trials[ibest - 1]
            middle = trials[ibest]
            right = trials[ibest + 1]
        interval = (right.Fcurv - left.Fcurv)
        if verbose:
            print(
                f"Geometry refinement "
                f"{iteration}: "
                f"Fcurv={middle.Fcurv:.9f}, "
                f"FWHM={middle.sigma:.9g} mm, "
                f"interval={interval:.9g}"
            )

        if interval <= Fcurv_tolerance:
            break
    # ===============================================================
    # Best actually evaluated solution
    # ===============================================================
    best = min(history,key=lambda item:item.sigma)
    if verbose:
        print(
            "Final geometry correction: "
            f"Fcurv={best.Fcurv:.9f}, "
            f"FWHM={best.sigma:.9g} mm, "
            f"centroid={best.centroid:.9g} mm, "
            f"rays={best.n_rays}"
        )
    return GeometryCorrectionResult(
        Fcurv=best.Fcurv,
        sigma=best.sigma,
        centroid=best.centroid,
        n_rays=best.n_rays,
        final_beam=best.final_beam,
        payload=best.payload,
        history=tuple(history)
    )