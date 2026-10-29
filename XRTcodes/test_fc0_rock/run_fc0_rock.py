"""TEST wrapper (not part of the scripts): run XRT_<name>.py with
 - bentSi_Crys depth reference at the entrance surface (PP_center_fraction = 0, test_fc0_rock/BentCrys2D.py
   that honours the parameter; otherwise identical to XRTcodes/BentCrys2D.py)
 - before each energy's run_process: rock the Laue pitch (offset -60..+160 urad, 100k rays,
   Fcurv = 1) and keep the offset that maximizes the intensity after the Laue crystal."""
import sys, numpy as np
import os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)  # BentCrys2D.py that honours PP_center_fraction
sys.path.insert(1, os.path.dirname(HERE))  # XRTcodes: the XRT_*.py scripts
import BentCrys2D
_init = BentCrys2D.bentSi_Crys.__init__
def init(self, *a, **k):
    k["PP_center_fraction"] = 0.0; _init(self, *a, **k)
BentCrys2D.bentSi_Crys.__init__ = init
import xrt.backends.raycing.sources as rs
import xrt.backends.raycing.run as rr
M = __import__(sys.argv[1])
original_run_process = M.run_process
def run_process(bl):
    nr = bl.source.nrays; bl.source.nrays = 100000
    b0 = bl.source.shine(); bl.mask1.propagate(b0); bl.source.nrays = nr
    bB, _ = bl.Bragg.reflect(rs.Beam(copyFrom=b0))
    p0 = bl.Laue.pitch; F0 = bl.Laue.Fcurv; bl.Laue.Fcurv = 1.0
    offs = np.arange(-60, 161, 5) * 1e-6; I = []
    for d in offs:
        bl.Laue.pitch = p0 + d
        bL, _ = bl.Laue.reflect(rs.Beam(copyFrom=bB)); g = bL.state == 1
        I.append((bL.Jss + bL.Jpp)[g].sum())
    I = np.array(I); k = int(np.argmax(I))
    bl.Laue.pitch = p0 + offs[k]; bl.Laue.Fcurv = F0
    print("ROCK Laue pitch offset %+.0f urad: intensity %.2f (at 0 offset %.2f)" % (offs[k]*1e6, I[k], I[np.argmin(np.abs(offs))]))
    return original_run_process(bl)
M.run_process = run_process
rr.run_process = run_process
M.main()
