"""Figures for the slides: flux and Fcurv vs energy, old / fixed / test runs."""
import glob, os, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

X = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "XRTcodes")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
COL = {"old": "#2a78d6", "fixed": "#eb6834", "test": "#1baf7a"}
LAB = {"old": "old scripts (before the fixes)", "fixed": "scripts after the fixes",
       "test": "test: after the fixes + entrance-surface depth + Laue rocking"}
MRK = {"old": "o", "fixed": "s", "test": "^"}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b",
                     "xtick.color": "#52514e", "ytick.color": "#52514e"})

def res(d):
    f = glob.glob(os.path.join(X, d, "res_*_energy.txt"))
    return np.loadtxt(f[0], skiprows=1, ndmin=2) if f else None

def fcurv(d):
    t = open(os.path.join(X, d, "run.log")).read()
    return np.array([float(v) for v in re.findall(r"Fcurv=([0-9.]+)", t)])

BL = [("SSU_SBDL", "SSU_SBDL (double Laue, PP)"), ("BNL_SBDL", "BNL_SBDL (double Laue, ML)"),
      ("SSU_BLM", "SSU_BLM (Bragg + Laue, ML)"), ("APS_BLM", "APS_BLM (Bragg + Laue, ML)")]
SUFFIX = {"old": "", "fixed": "_fixed", "test": "_test_fc0_rock"}

for quantity in ("flux", "fcurv"):
    fig, axs = plt.subplots(2, 2, figsize=(9, 5.2) if quantity == "fcurv" else (7, 6.4),
                            sharey=(quantity == "fcurv"))
    for ax, (name, title) in zip(axs.flat, BL):
        for k in ("old", "fixed", "test"):
            d = "run_XRT_" + name + SUFFIX[k]
            r = res(d)
            if r is None:
                continue
            y = r[:, 2] if quantity == "flux" else fcurv(d)
            ax.plot(r[:, 0] / 1e3, y, color=COL[k], lw=2, marker=MRK[k], ms=5,
                    markeredgecolor="white", markeredgewidth=1, label=LAB[k])
        ax.set_title(title, fontsize=10, loc="left", color="#0b0b0b")
        ax.grid(axis="y", color="#e6e5e0", lw=0.8)
        if quantity == "flux":
            ax.set_yscale("log")
            ax.set_ylabel("flux at sample (ph/s)")
        else:
            ax.axhline(1.0, color="#9a9994", lw=1, ls="--")
            ax.set_ylabel("Fcurv")
        ax.set_xlabel("E (keV)")
    h, l = [], []
    for ax in axs.flat:
        for hh, ll in zip(*ax.get_legend_handles_labels()):
            if ll not in l:
                h.append(hh); l.append(ll)
    if quantity == "flux":
        fig.legend(h, l, loc="upper left", ncol=1, frameon=False, fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.87))
    else:
        fig.legend(h, l, loc="upper center", ncol=3, frameon=False, fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(OUT, quantity + ".pdf"))
    fig.savefig(os.path.join(OUT, quantity + ".png"), dpi=110)

# Laue rocking curve at 80 keV, SSU_BLM crystal: the last plot of
# XRTcodes/XOPPY/comparison_SSU_BLM.py, from the files it writes
D = os.path.join(X, "XOPPY", "SSU_BLM_radii")
pp = np.loadtxt(os.path.join(D, "diff_pat_PP.dat"), skiprows=5)
ml = np.loadtxt(os.path.join(D, "diff_pat_ML.dat"), skiprows=5)
internal = np.loadtxt(os.path.join(D, "comparison_internal.dat"))
fig, ax = plt.subplots(figsize=(8, 4.2))
curves = [(pp[:, 0], pp[:, -1], "XOPPY diff_pat (PP)", "#2a78d6", "-"),
          (ml[:, 0], ml[:, -1], "XOPPY diff_pat (ML)", "#eb6834", "-"),
          (internal[:, 0], internal[:, 1], "internal PP", "#1baf7a", "--"),
          (internal[:, 0], internal[:, 3], "internal ML", "#eda100", "--")]
for x, y, label, color, ls in curves:
    ax.plot(x, y, color=color, lw=2, ls=ls, label="%s: %.1f $\\mu$rad" % (label, np.trapezoid(y, x)))
ax.set_xlim(-100, 100)
ax.set_xlabel(r"$\theta - \theta_B$ (incident) ($\mu$rad)")
ax.set_ylabel("reflectivity (s-polarization)")
ax.grid(axis="y", color="#e6e5e0", lw=0.8)
ax.legend(frameon=False, loc="upper left", fontsize=9, title="integrated reflectivity", title_fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "laue_xoppy.pdf"))
fig.savefig(os.path.join(OUT, "laue_xoppy.png"), dpi=110)

# SSU_SBDL (scripts after the fixes) vs SSU_BLM (test run): same source, same 56 m screen
sbdl = res("run_XRT_SSU_SBDL_fixed")
blm = res("run_XRT_SSU_BLM_test_fc0_rock")
SERIES = [(sbdl, "SSU_SBDL: double Laue, PP (after the fixes)", "#2a78d6", "o"),
          (blm, "SSU_BLM: Bragg + Laue, ML (test run)", "#1baf7a", "^")]
PANELS = [(lambda r: r[:, 2], "flux at sample (ph/s)", True),
          (lambda r: r[:, 3], "flux density (ph/s/mm$^2$)", True),
          (lambda r: r[:, 8] / r[:, 0] * 1e5, r"bandwidth $\Delta E/E$ (%)", True),
          (lambda r: r[:, 5], "vertical FWHM (mm)", True)]
fig, axs = plt.subplots(2, 2, figsize=(7, 6.4))
for ax, (f, ylabel, logy) in zip(axs.flat, PANELS):
    for r, label, color, marker in SERIES:
        ax.plot(r[:, 0] / 1e3, f(r), color=color, lw=2, marker=marker, ms=5,
                markeredgecolor="white", markeredgewidth=1, label=label)
    if logy:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.set_xlabel("E (keV)")
    ax.grid(axis="y", color="#e6e5e0", lw=0.8)
h, l = axs.flat[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper left", ncol=1, frameon=False, fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.savefig(os.path.join(OUT, "ssu_sbdl_vs_blm.pdf"))
fig.savefig(os.path.join(OUT, "ssu_sbdl_vs_blm.png"), dpi=110)

# XOPPY integrated reflectivities of the SSU crystals (XRTcodes/XOPPY/integrated_SSU.py)
xi = np.loadtxt(os.path.join(X, "XOPPY", "SSU_integrated", "integrated_SSU.dat"))
E = xi[:, 0] / 1e3
B, L311, L111 = xi[:, 1], xi[:, 2], xi[:, 3]
fig, axs = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True)
ax = axs[0]
for y, label, color, marker in ((B, "Bragg Si(311), ML, 10 mm (SSU_BLM)", "#2a78d6", "o"),
                                (L311, "Laue Si(311), PP, 1 mm (SSU_BLM)", "#eb6834", "s"),
                                (L111, "Laue Si(111), PP, 1 mm (SSU_SBDL)", "#1baf7a", "^")):
    ax.plot(E, y, color=color, lw=2, marker=marker, ms=5, markeredgecolor="white", markeredgewidth=1, label=label)
ax.set_yscale("log")
ax.set_ylabel(r"integrated reflectivity ($\mu$rad)")
ax.set_ylim(2, 400)
ax.grid(axis="y", color="#e6e5e0", lw=0.8)
ax.legend(frameon=False, fontsize=9, loc="lower right")
ax = axs[1]
ratio_tracing = np.array([blm[np.argmin(abs(blm[:, 0] - e)), 2] / sbdl[np.argmin(abs(sbdl[:, 0] - e)), 2]
                          for e in xi[:, 0]])
ax.plot(E, L311 * B / L111**2, color="#2a78d6", lw=2, marker="o", ms=5, markeredgecolor="white",
        markeredgewidth=1, label=r"XOPPY: $(L_{311}\,B_{311}) / (L_{111}\,L_{111})$")
ax.plot(E, ratio_tracing, color="#eb6834", lw=2, marker="s", ms=5, markeredgecolor="white",
        markeredgewidth=1, label="ray tracing: flux SSU_BLM (test run) / SSU_SBDL")
ax.axhline(1.0, color="#9a9994", lw=1, ls="--")
ax.set_yscale("log")
ax.set_ylabel("SSU_BLM / SSU_SBDL")
ax.set_xlabel("E (keV)")
ax.grid(axis="y", color="#e6e5e0", lw=0.8)
ax.legend(frameon=False, fontsize=9, loc="upper left")
fig.tight_layout()
fig.savefig(os.path.join(OUT, "integrated_ssu.pdf"))
fig.savefig(os.path.join(OUT, "integrated_ssu.png"), dpi=110)
