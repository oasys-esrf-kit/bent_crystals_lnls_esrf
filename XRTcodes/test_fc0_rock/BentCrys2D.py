

import inspect
import numpy as np
import xrt.backends.raycing.materials as rm
try:  # xrt 1.x (with the patched oes_base that passes ray directions to get_amplitude)
    from xrt.backends.raycing.oes_base import OE
except ImportError:  # xrt >= 2.0
    from xrt.backends.raycing.oes import OE
from xrt.backends.raycing.physconsts import CH
from xrt.backends import raycing

class BentCrys2D(OE):
    # Bend crystal in Bragg or Laue with Rs, Rm

    def __init__(self, *args, **kwargs):
        """
        *Rm* and *Rs*: float
            Meridional and sagittal radii.
        """
        kwargs = self.pop_kwargs(**kwargs)
        OE.__init__(self, *args, **kwargs)

    def pop_kwargs(self, **kwargs):
        self.Rm = kwargs.pop('Rm', 51000.)  # R meridional
        self.Rs = kwargs.pop('Rs', None)  # R sagittal
        self.Fcurv = kwargs.pop('Fcurv',None) # curvature correction factor for asymmetric case
        if self.Rs is None:
            self.Rs = 1e16
        return kwargs

    def local_z(self, x, y):
        """Determines the surface of OE at (*x*, *y*) position."""
        z = 0.5 * x ** 2 / self.Rs + 0.5* y ** 2 / self.Rm
        return z
    
    def local_n_depth(self, x, y, z):
        # doesnt work for bent crystal using PP or ML method
        
        n0 = self.local_n_saddle(x, y, self.Rm, self.Rs, self.alpha)

        return n0

    def local_n_saddle(self, x, y, Rm, Rs, alpha):
        # Fcurv is just a curvature factor to obtain the focus at theoretical position
        # to facilitate simulation and analysis
        # Physical Rs in reality keeps without applying that factor
        # the rest is identically to the BentLaue2D function of XRT
        a = -x / (Rs * self.Fcurv)
        b = -y / Rm
        c = 1.0

        norm = np.sqrt(a**2 + b**2 + 1)
        a /= norm
        b /= norm
        c /= norm

        sinpitch = -b
        cospitch = np.sqrt(1 - b**2)

        sinroll = -a
        cosroll = np.sqrt(1 - a**2)

        #base Bragg normal 
        aB = np.zeros_like(a)
        bB = np.ones_like(a)
        cB = np.zeros_like(a)

        if alpha:
            bB, cB = raycing.rotate_x(bB, cB, self.cosalpha, -self.sinalpha)
        
        aB, cB = raycing.rotate_y(aB, cB, cosroll, -sinroll)
        bB, cB = raycing.rotate_x(bB, cB, cospitch, sinpitch)

        normB = (bB**2 + cB**2 + aB**2)**0.5
        aB /=normB
        bB /=normB
        cB /=normB
      
        return [aB, bB, cB, a, b, c]
        
    def local_n(self, x, y):
        return self.local_n_saddle(x, y, self.Rm, self.Rs, self.alpha)

# Crystal
class bentSi_Crys(rm.CrystalSi):
    def __init__(self, *args, **kwargs):
        kwargs = self.__pop_kwargs(**kwargs)
        # Silicon density [g/cm^3]
        kwargs.setdefault("rho", 2.3296)
        rm.CrystalSi.__init__(self, *args, **kwargs)

    def __pop_kwargs(self, **kwargs):
        self.Rs = kwargs.pop("Rs")
        self.Rm = kwargs.pop("Rm")
        self.aLrad = kwargs.pop("aLrad")
        self.method = kwargs.pop("method")
        self.ML_sign = kwargs.pop("ML_sign", 1.0)
        # PP deviation-parameter centering.
        # 0.0: original entrance-surface convention
        # 0.5: geometrical Bragg condition placed at the crystal mid-plane
        self.PP_center_fraction = kwargs.pop("PP_center_fraction",0.5,)
        self.ST = kwargs.pop("ST",None)

        return kwargs

    def _rays_from_context(self, xd, yd, beamInDotNormal):
        # Rebuild the arguments that the patched xrt 1.x oes_base passed to get_amplitude,
        # from the context recorded by the _reflect_local hook (see end of file).
        context = getattr(self, "_ray_context", None)
        if context is None or "a_in" not in context:
            raise RuntimeError("bentSi_Crys: no ray context recorded by the _reflect_local hook")
        lb = context["lb"]
        goodN = lb.state == 1
        if not (np.array_equal(context["x"], xd) and np.array_equal(context["y"], yd) and
                np.array_equal(lb.x[goodN], xd)):
            raise RuntimeError("bentSi_Crys: ray context does not match the rays to diffract")
        ad, bd, cd = context["a_in"], context["b_in"], context["c_in"]
        oeNormal = context["oeNormal"]
        inDotNormal = np.clip(ad * oeNormal[0] + bd * oeNormal[1] + cd * oeNormal[2], -1, 1)
        if not np.allclose(inDotNormal, beamInDotNormal, rtol=0, atol=1e-9):
            raise RuntimeError("bentSi_Crys: recorded directions/normals inconsistent with xrt")
        outBN = [lb.a[goodN], lb.b[goodN], lb.c[goodN]]
        return ad, bd, cd, oeNormal, outBN

    def get_amplitude(self, E, beamInDotSurfaceNormal, beamOutDotSurfaceNormal, beamInDotNormal, 
                             xd=None, yd=None, ad=None, bd=None, cd=None,
                             oeNormal=None,outBN=None):

        def build_deviation_reference(polar_factor,beta,csmall,ML_refinement):
            """
            Build consistent PP and ML deviation coordinates.
            The common reference is the deviation at the crystal mid-plane:
                eta_center = y_center
            PP:
                continuous entrance and exit coordinates over the true thickness t.
            ML:
                discrete lamella-center coordinates, symmetrically positioned around
                the same center.
            """
            # ------------------------------------------------------------
            # Common deviation reference at the crystal mid-plane
            # ------------------------------------------------------------
            center = deviation_common / polar_factor
            # ------------------------------------------------------------
            # PP: continuous interval across the true crystal thickness
            # ------------------------------------------------------------
            #
            # eta_out - eta_in = -beta*t
            #
            pp_delta_total = -beta * t
            fc = self.PP_center_fraction
            pp_in = center - fc * pp_delta_total
            pp_out = pp_in + pp_delta_total
            # ------------------------------------------------------------
            # ML: determine discrete lamella spacing
            # ------------------------------------------------------------
            if is_bragg:
                step_abs = 2.0
            else:
                step_abs = np.pi / 2.0

            # ============================================================
            # ML lamella discretization
            # ============================================================
            # Target normalized deviation step.
            ml_step_target = np.copysign(step_abs,csmall)
            # Target lamella thickness implied by the ML step criterion.
            cry_a_target = np.abs(ml_step_target / csmall)
            cry_t_target = np.median(cry_a_target * ext_l_s / polar_factor)
            # Number of lamellae required to cover the real crystal.
            nslab = max(1,int(np.ceil(t / cry_t_target)))
            # additional convergence refinement
            nslab *= ML_refinement
            cry_t = t / nslab
            ml_step_abs = np.median(np.abs(csmall * cry_t * polar_factor / ext_l_s))
            ml_step = np.copysign(ml_step_abs,np.median(csmall))
            # IMPORTANT:
            # force total ML thickness to equal the real crystal thickness.
            cry_t = (t / nslab)
            # Actual normalized deviation step corresponding to
            # this exact physical lamella thickness.
            ml_step_abs = np.median(np.abs(csmall * cry_t * polar_factor / ext_l_s))
            # Preserve bending/deviation sign.
            ml_step = np.copysign(ml_step_abs,np.median(csmall))
            # ------------------------------------------------------------
            # ML lamella centers
            # ------------------------------------------------------------
            ml_first_center = (center - fc * (nslab - 1) * ml_step)
            ml_last_center = (ml_first_center + (nslab - 1) * ml_step)
            # Geometrical boundaries surrounding the discrete lamella centers.
            ml_boundary_in = (ml_first_center - 0.5 * ml_step)
            ml_boundary_out = (ml_boundary_in + nslab * ml_step)

            return {
                "center": center,
                "pp_in": pp_in,
                "pp_out": pp_out,
                "pp_delta_total": pp_delta_total,
                "ml_step": ml_step,
                "ml_first_center": ml_first_center,
                "ml_last_center": ml_last_center,
                "ml_boundary_in": ml_boundary_in,
                "ml_boundary_out": ml_boundary_out,
                "cry_t": cry_t,
                "nslab": nslab,
            }
        
        def calc_csmall(polar_factor):
            return (
                dalpha_over_dt
                * bfac
                * (0.5 * ext_l_s / polar_factor)
                / (
                    polar_factor
                    * abs_PsiH
                    * sqrt_abs_bfac
                )
            ) / polar_factor**2
        
        if ad is None:  # xrt >= 2.0 does not pass ray directions and normals
            ad, bd, cd, oeNormal, outBN = self._rays_from_context(xd, yd, beamInDotNormal)

        is_bragg = "Bragg" in self.geom
        #is_laue  = "Laue"  in self.geom
        method = self.method

        # basic parameters for reflectivity calculation
        mu = self.get_absorption_coefficient(E)  # cm^-1
        mu /= 10 # [mm^-1]
        waveLength = CH / E # the word "lambda" is reserved
        lam_mm = waveLength * 1e-7

        # get Psi-factors
        dhkl = self.d #crystal d-spacing [A]
        F0, Fhkl, Fhkl_, Psi0, PsiH, PsiH_ = self.get_F_chi(E, 0.5 / dhkl)
        thetaB = np.arcsin(waveLength / dhkl / 2.) #Bragg angle in rad
        aLrad = self.aLrad # asymmetric angle
        # extinction depth [mm] - s,p - depending on E
        ext_l_s, ext_l_p = rm.CrystalSi.get_extinction_lengths(self,E)
        ext_l_s *= np.sin(thetaB) 
        ext_l_p *= np.sin(thetaB) 
        
        #polarisation factor [s,p]
        pfac_p = np.abs(np.cos(2.0 * thetaB))
        pfac_p0 = float(np.median(pfac_p)) 
        pfac = [1.0,pfac_p0]
        
        # Radii and thickness of crystals
        Rs = self.Rs # Rs [mm]
        Rm = self.Rm # Rm [mm]
        t = self.t # crystal thickness [mm]

        # calculate compliance matrix, G-vector and gamma from incoming and outgoing vectors
        # get compliance matrix
        # S-tensor.txt is created by the main code:
        # CT, ST, ET, GT, poisson = stiffness_tensor(C0=[165.7,63.9,79.6], hkl=nrm, sag=sag, unitfac=unitfac)
        # poL = poisson[1,2]  # ST'(ij)= ST(ij) / ST(33) - poL=ST'(23)
        # Sf = 'S_tensor.txt'
        # np.savetxt(Sf, ST)
        #ST = np.loadtxt('S_tensor.txt')
        ST=self.ST
        if ST is None:
            ST = np.loadtxt('S_tensor.txt')
        S11 = ST[0,0] 
        S12 = ST[0,1] 
        S13 = ST[0,2] 
        S14 = ST[0,3] 
        S21 = ST[1,0] 
        S22 = ST[1,1] 
        S23 = ST[1,2] 
        S24 = ST[1,3] 
        S31 = ST[2,0] 
        S32 = ST[2,1] 
        S41 = ST[3,0] 
        S42 = ST[3,1] 

        # Input, output, Bragg-plane directions cosine 
        Vbh = np.vstack((oeNormal[0], oeNormal[1], oeNormal[2]))
        Vbh /= np.sqrt(np.sum(Vbh * Vbh, axis=0))

        Vin = np.vstack((ad,bd,cd))
        Vin /= np.sqrt(np.sum(Vin * Vin, axis=0))

        # calculate Bragg-plane normal for any case
        beamInDotH_raw = np.einsum("ij,ij->j",Vbh,Vin)
        median_beamInDotH = np.median(beamInDotH_raw[np.isfinite(beamInDotH_raw)])
        if median_beamInDotH > 0.0:
            reciprocal_sign = -1.0
        else:
            reciprocal_sign = 1.0
        Vbh_diff = reciprocal_sign * Vbh
        beamInDotH = np.einsum("ij,ij->j",Vbh_diff,Vin)

        Vout = np.vstack(((outBN[0]), (outBN[1]), (outBN[2])))
        Vout /= np.sqrt(np.sum(Vout * Vout, axis=0))

        lamb_over_d = waveLength / dhkl

        # surface Normals = oeNormal[3:6]
        if len(oeNormal) == 6:
            Vsurf = np.vstack((oeNormal[-3], oeNormal[-2], oeNormal[-1]))
            Vsurf /= np.sqrt(np.sum(Vsurf * Vsurf, axis=0))
        else: # alpha=0 symmetric case
            Vsurf = np.vstack((oeNormal[0], oeNormal[1], oeNormal[2]))
            Vsurf /= np.sqrt(np.sum(Vsurf * Vsurf, axis=0))
            Vbh = Vsurf

        gamma0 = np.einsum('ij,ij->j', Vsurf, Vin)
        gammaH = np.einsum('ij,ij->j', Vsurf, Vout)
        bfac = gamma0 / gammaH   

        abs_bfac = np.abs(bfac)
        sqrt_abs_bfac = np.sqrt(abs_bfac)
        abs_PsiH = np.abs(PsiH)
        abs_PsiH_2 = np.abs(PsiH_)**2
        abs_gammaH = np.abs(gammaH)
        abs_gamma0 = np.abs(gamma0)
        sqrt_abs_gamma = np.sqrt(np.abs(gamma0 * gammaH))
        sqrt_abs_PsiH_ = np.sqrt(np.abs(PsiH * PsiH_))

        # alphaZ is positive always using the beamInDotH definition
        alphaZ = lamb_over_d**2 + 2.0 * lamb_over_d * beamInDotH

        # Common susceptibility denominator used by both PP and ML.
        chi_h_common = np.sqrt(PsiH * PsiH_)

        # Common normalized deviation before applying polarization.
        deviation_common = (((1.0 - bfac) * Psi0 + bfac * alphaZ) / (2.0 * sqrt_abs_bfac * chi_h_common)).real

        #calculate bending moments from Rs, Rm and ST
        compliance_fac = 1 / (S12 * S21 - S11 * S22)
        M1 = (S12 / Rm - S22 / Rs) * compliance_fac
        M2 = (S21 / Rs - S11 / Rm) * compliance_fac  

        #for ML
        gamma0_ml = beamInDotSurfaceNormal
        A1 = (-Vout[1] * Vin[1] * (Vout[2] - Vin[2]) + Vout[1] * Vin[2] * (Vout[1] - Vin[1]) + Vout[2] * Vin[1] * (Vout[1] - Vin[1]))
        A2 = Vout[2] * Vin[2] * (Vout[2] - Vin[2])
        A3 = Vout[2] * Vin[2] * (Vout[1] - Vin[1])
        VS = (A1 * (S21*M1 + S22*M2) + A2 * (S31*M1 + S32*M2) + A3 * (S41*M1 + S42*M2))
        # dimensionless curvature c
        dalpha_over_dt = -2.0 * VS / gamma0_ml # ml_c factor as calculated in crystal3
        dalpha_over_dt *= self.ML_sign

        # for PP
        # calculating G-vector
        G1 = Vin[2] * Vout[2] * Vbh_diff[2]
        G2 = Vin[2] * Vout[2] * Vbh_diff[1]
        G3 = (Vin[1] * Vout[2] * Vbh_diff[1] + Vin[2] * Vout[1] * Vbh_diff[1] - Vin[1] * Vout[1] * Vbh_diff[2])
        # continue G- vector
        wHd = 1.0 / (dhkl * 1e-7)
        G = wHd * (G1 * (M2 * S23 + M1 * S13) + G2 * (M2 * S24 + M1 * S14) + G3 * (M2 * S22 + M1 * S12))                
        beta_s = G * lam_mm / (pfac[0] * abs_PsiH * sqrt_abs_gamma)
        beta_p = G * lam_mm / (pfac[1] * abs_PsiH * sqrt_abs_gamma)

        csmall_s = calc_csmall(pfac[0])
        csmall_p = calc_csmall(pfac[1])

        reference_s = build_deviation_reference(polar_factor=pfac[0],beta=beta_s,csmall=csmall_s,ML_refinement=1)
        reference_p = build_deviation_reference(polar_factor=pfac[1],beta=beta_p,csmall=csmall_p,ML_refinement=1)
        references = (reference_s,reference_p)

        if method == 'PP': # use Penning-Polder method PP
            # ONLY WITH LAUE CRYSTAL!
        
            curS = []
            curP = []
            
            for pfaci, pfacn in enumerate(pfac):
                # pendelloesung length in intensity
                pend_len = lam_mm * abs_gamma0 / (pfacn * abs_PsiH) / sqrt_abs_bfac
                # beta as defined by Xianbo
                beta = G * lam_mm / (pfacn * abs_PsiH * sqrt_abs_gamma)
                betac = np.pi / (2 * pend_len) # critical strain gradient

                reference = references[pfaci]
                if pfaci == 0:
                    beta = beta_s
                else:
                    beta = beta_p
                #eta_0 = reference["center"]
                eta_in = reference["pp_in"]
                eta_out = reference["pp_out"]
                root_in = np.sqrt(eta_in**2 + bfac + 0j)
                root_out = np.sqrt(eta_out**2 + bfac + 0j)

                cep = eta_out + root_out
                cem = eta_out - root_out
                cip = eta_in + root_in
                cim = eta_in - root_in
                
                limit_Psi0_Img=1e-10
                if aLrad < np.pi:
                    limit_Psi0_Img = 1e-10
                else:
                    limit_Psi0_Img = 1e-12

                eps = []
                for psih, psih_, psi0 in zip(PsiH, PsiH_, Psi0):
                    if abs(psi0.imag) > limit_Psi0_Img:
                        eps.append(((psih * psih_) ** 0.5).imag / psi0.imag)
                    else:
                        eps.append(1.0)
        
                eloss = []
                for betac_th, beta_th in zip(betac, beta):
                    if abs(beta_th) > 1e-6:
                        eloss.append(np.exp(-2 * np.pi * betac_th / abs(beta_th)))
                    else:
                        eloss.append(0.0)

                for beta_th, cep_th, cem_th, cip_th, cim_th, eloss_th, eps_th, bfac_th, gamma0_th, mu_th in zip(
                        beta, cep, cem, cip, cim, eloss, eps, bfac, gamma0, mu):
                    
                    if np.real(cep_th * cip_th) > 0.0:
                        IRp_raw = np.real(cep_th**2 * bfac_th / (cep_th**2 + bfac_th) / (cip_th**2 + bfac_th))

                        if abs(beta_th) > 1e-6:
                            B22p = np.real((bfac_th - 1.0) * (cep_th - cip_th) / (2.0 * beta_th * t))
                            D22p = np.real(bfac_th * eps_th * np.log(cep_th / cip_th) / (beta_th * t))
                        else:
                            B22p = 0.0
                            D22p = np.real(-2.0 * cip_th / (1.0 + cip_th * cip_th))

                        absorption_p = np.real(np.exp(-mu_th * t * (1.0 + B22p + D22p) / abs(gamma0_th)))
                    else:
                        IRp_raw = 0.0
                        B22p = 0.0
                        D22p = 0.0
                        absorption_p = 1.0

                    if np.real(cem_th * cim_th) > 0.0:
                        IRm_raw = np.real(cem_th**2 * bfac_th / (cem_th**2 + bfac_th) / (cim_th**2 + bfac_th))

                        if abs(beta_th) > 1e-6:
                            B22m = np.real((bfac_th - 1.0) * (cem_th - cim_th) / (2.0 * beta_th * t))
                            D22m = np.real(bfac_th * eps_th * np.log(cem_th / cim_th) / (beta_th * t))
                        else:
                            B22m = 0.0
                            D22m = np.real(-2.0 * cim_th / (1.0 + cim_th * cim_th))

                        absorption_m = np.real(np.exp(-mu_th * t * (1.0 + B22m + D22m) / abs(gamma0_th)))
                    else:
                        IRm_raw = 0.0
                        B22m = 0.0
                        D22m = 0.0
                        absorption_m = 1.0

                    IRp = IRp_raw * absorption_p 
                    IRm = IRm_raw * absorption_m 
                    Reflectivity = (np.real(IRp + IRm) * (1.0 - eloss_th))

                    if pfaci == 0:
                        curS.append(Reflectivity)  
                    else: 
                        curP.append(Reflectivity)
              
            curS = np.asarray(curS, dtype=float)
            curP = np.asarray(curP, dtype=float)

            # curS and curP are intensity reflectivities.
            # xrt get_amplitude() expects field amplitudes.
            curS = np.nan_to_num(curS, nan=0.0, posinf=0.0, neginf=0.0)
            curP = np.nan_to_num(curP, nan=0.0, posinf=0.0, neginf=0.0)
            curS = np.clip(curS, 0.0, 1.0)
            curP = np.clip(curP, 0.0, 1.0)
            curS = np.sqrt(curS)
            curP = np.sqrt(curP)
        else: # use multilamellar method ML
        
            def ml_pol(polarFactor,reference):
             
                # the ML formulas (crystal3.f90) use susceptibilities with negative imaginary part;
                # xrt returns their complex conjugates
                qsmall = (bfac * np.conj(PsiH) * np.conj(PsiH_) * polarFactor**2)
                g = ((1.0 - bfac) * np.conj(Psi0).imag / (2.0 * polarFactor * sqrt_abs_bfac * sqrt_abs_PsiH_))
                deltay = reference["ml_step"]
                cry_t = reference["cry_t"]
                nslab = reference["nslab"]
                skh = cry_t / abs_gammaH
                asmall = (np.pi * cry_t / lam_mm / abs_gamma0)
                # Common center reference.
                y0 = reference["center"]
                # First evaluated lamella center.
                y = reference["ml_first_center"].copy()
                Rtot = np.zeros_like(y0,dtype=float)
                prodt = np.ones_like(y0,dtype=float)
                tsmall_old = np.ones_like(y0,dtype=float)
                ztmp = polarFactor * abs_PsiH * sqrt_abs_bfac
                
                for i in range(1, nslab + 1):
                    z = (y + 1j * g) * ztmp
                    qz = qsmall + z*z
                    usmall = np.sqrt(qz)
                    qbig = np.abs(qz)
                    zabsq = np.abs(z)**2
                    zu = np.conjugate(z) * usmall
                    w = usmall.imag
                    v = usmall.real

                    RC1 = np.exp(-2.0 * asmall * z.imag * (bfac + 1.0) / (bfac - 1.0))

                    s2 = np.sin(asmall * v)**2
                    aw = np.clip(asmall * w, -100.0, 100.0)
                    sh2 = np.sinh(aw)**2
                    d1 = qbig + (qbig + zabsq) * np.abs(sh2)
                    d2 = (qbig - zabsq) * s2
                    d3 = np.real(-zu) * np.sinh(2.0 * aw)
                    d4 = np.imag(zu) * np.sin(2.0 * asmall * v)

                    if is_bragg:
                        denom = np.real(d1 - d2 + d3 + d4)
                        ufr = ((s2 + np.abs(sh2)) * abs_bfac * abs_PsiH_2)
                        rsmall = np.real(polarFactor**2 * ufr / denom)
                        rsmall = np.nan_to_num(rsmall, nan=0.0, posinf=0.0, neginf=0.0)
                        rsmall = np.clip(rsmall, 0.0, 1.0)
                        tsmall = qbig * RC1 / denom
                        prodt *= tsmall_old
                        # absorption for one lamella traversal
                        att = np.exp(-mu * skh * (i-1))
                        att = np.nan_to_num(att, nan=0.0, posinf=0.0, neginf=0.0)
                        att = np.clip(att, 0.0, 1.0)
                        # ray reaches this lamella, reflects, and exits
                        Rtot += rsmall * prodt * att
                        tsmall_old = tsmall
                        y += deltay
                        
                    else: # Laue
                        denom = np.real(d1 - d2 - d3 - d4)
                        ufr = (s2 + np.abs(sh2)) * abs_bfac * np.abs(PsiH * PsiH_)
                        rsmall = polarFactor**2 * ufr / qbig
                        rsmall = np.nan_to_num(rsmall, nan=0.0, posinf=0.0, neginf=0.0)
                        rsmall = np.clip(np.real(rsmall), 0.0, 1.0)
                        tsmall = denom * RC1 / qbig
                        prodt *= tsmall_old
                        att = np.exp(-mu * skh * (nslab-i))
                        att = np.nan_to_num(att, nan=0.0, posinf=0.0, neginf=0.0)
                        att = np.clip(att, 0.0, 1.0)
                        Rtot += rsmall * prodt * att
                        tsmall_old = tsmall
                        y += deltay

                Rtot = np.nan_to_num(Rtot, nan=0.0, posinf=0.0, neginf=0.0)
                return np.real(np.clip(Rtot, 0.0, 1.0))

            curS = ml_pol(pfac[0],reference_s)
            curP = ml_pol(pfac[1],reference_p)
            curS = np.sqrt(curS)
            curP = np.sqrt(curP)
            #print("curS mean =", np.mean(curS))
            #print("curP mean =", np.mean(curP))
        return np.array(curS), np.array(curP)


def _install_ray_context_hook():
    """xrt >= 2.0 calls Crystal.get_amplitude only with (E, beamInDotSurfaceNormal,
    beamOutDotSurfaceNormal, beamInDotNormal, x, y), but bentSi_Crys also needs the incoming
    and outgoing ray directions and the Bragg/surface normals (the patched xrt 1.x oes_base
    passed them). Wrap OE._reflect_local so that, when the material is a bentSi_Crys, the
    local_n call on the diffracting rays records them (directions are then already in the
    local frame and not yet reflected). Applies to any OE class (e.g. JohannCylinder).
    Hooks installed by different copies of this module chain (each handles its own class)."""
    owner = next(c for c in OE.__mro__ if "_reflect_local" in vars(c))
    original = owner._reflect_local
    if getattr(original, "_bent_crys_class", None) is bentSi_Crys:  # already installed
        return
    signature = inspect.signature(original)

    def _reflect_local(self, *args, **kwargs):
        bound = signature.bind(self, *args, **kwargs)
        material = bound.arguments.get("material")
        materials = material if raycing.is_sequence(material) else [material]
        crystals = [m for m in materials if isinstance(m, bentSi_Crys)]
        if not crystals:
            return original(self, *args, **kwargs)

        lb = bound.arguments["lb"]
        local_n = bound.arguments.get("local_n") or self.local_n
        context = {"lb": lb}

        def recording_local_n(x, y):
            normals = local_n(x, y)
            goodN = lb.state == 1
            if np.size(x) == goodN.sum() and np.array_equal(x, lb.x[goodN]) and \
                    np.array_equal(y, lb.y[goodN]):
                context.update(x=np.array(x), y=np.array(y),
                               a_in=lb.a[goodN].copy(), b_in=lb.b[goodN].copy(),
                               c_in=lb.c[goodN].copy(),
                               oeNormal=np.asarray([n * np.ones_like(x) for n in normals]))
            return normals

        bound.arguments["local_n"] = recording_local_n
        for crystal in crystals: crystal._ray_context = context
        try:
            return original(*bound.args, **bound.kwargs)
        finally:
            for crystal in crystals: crystal._ray_context = None

    _reflect_local._bent_crys_class = bentSi_Crys
    owner._reflect_local = _reflect_local

_install_ray_context_hook()
