# calulates: stiffness tensor, E-modulus, Shear modulus, Poisson ratio of crystals
       
def stiffness_tensor(descriptor='Si', hkl=[1,0,0], sag=None, C0=[165.7,63.9,79.6], ansys=0, unitfac=1.):
    # descriptor: crystal material
    # Stiffness coefficient matrix (6,6) values
        # for Si: C11=165.7 GPa, C12=63.9 GPa, C44=79.6 GPa
    # hkl: deflection orientation of crystal (surface normal)
    # sag: sagittal axis used perpendicular to surface normal
    # ansys: if <>0 use ansys coordinate system
    # unit factor: if 1. units are SI, eg. [m] - use other factor for converting into cm or mm...

    import numpy as np
    from scipy.linalg import expm, norm
    from scipy.spatial.transform import Rotation as R

    hh, kk, ll = hkl

    # conventions
    
    # !standard definition! - as used in XRT (crystal defined in LOCAL ccordinate system)
    # !standard definition! - as used in SHADOW / XCRYSTAL
    # e1: direction vector norm(al) to diffraction plane ; z-axis
    # e2: in-plane along sagittal direction ; y-axis
    # e3: in-plane along meridional direction ; x-axis

    # !ANSYS definition!
    # e1: direction vector norm(al) to diffraction plane ; z-axis
    # e2: in-plane along sagittal direction ; x-axis
    # e3: in-plane along meridional direction ; y-axis
    # ansys=1 -> use ANSYS definition

    #calculate direction vectors perpendicular to the hkl (norm) vector
    if sag is None:
        hkl = np.array(hkl,dtype=float)
        e0_1 = np.where(hkl != 0)[0]
        e0_0 = np.where(hkl == 0)[0]
        if np.size(e0_0) > 0: 
            e1 = hkl
            e2 = np.array([0, 0, 0],dtype=float)
            e3 = np.array([0, 0, 0],dtype=float)
            if np.size(e0_1) == 1:
                e2[e0_0[0]] = hkl[e0_1[0]]
                e3[e0_0[1]] = hkl[e0_1[0]]
            else:
                e2[e0_1[0]] = hkl[e0_1[0]]
                e2[e0_1[1]] = -hkl[e0_1[1]]
                e3[e0_0[0]] = 1
        else:
            e1 = np.array([hh,kk,ll],dtype=float)
            e2 = np.array([0,ll,-kk],dtype=float)
            e3 = np.array([(kk*kk+ll*ll),-hh*kk,-hh*ll],dtype=float)
    else:
        #surface normal is hkl, sagittal axis sag is e2, e3 is cross product
        s1, s2, s3 = sag
        e1 = np.array([hh,kk,ll],dtype=float)
        e2 = np.array([s1,s2,s3],dtype=float)
        e3 = np.cross(hkl,sag)

    e1 = e1 / norm(e1)
    e2 = e2 / norm(e2)
    e3 = e3 / norm(e3)

    # Stiffness matrix Si(100)
    C11=C0[0]; C12=C0[1]; C44=C0[2]
    C100 = np.array([[C11,C12,C12,0,0,0],[C12,C11,C12,0,0,0],[C12,C12,C11,0,0,0],[0,0,0,C44,0,0],[0,0,0,0,C44,0],[0,0,0,0,0,C44]])

    mu=C11 - C12 - 2 * C44

    Ce=np.zeros((6,6))
    Ce[0,0] = mu; Ce[1,1] = mu; Ce[2,2] = mu
    
    if ansys == 0:
        M0 = np.array([e1*e1,e2*e2,e3*e3,e2*e3,e1*e3,e1*e2]) # common definition
    else:
        e1, e2, e3 = e3, e1, e2
        M0 = np.array([e1*e1,e2*e2,e3*e3,e1*e2,e2*e3,e3*e1]) # ANSYS definition

    # C is compliance matrix
    # units: GPa = 10E9 N/m²
    M = np.dot(M0,M0.transpose())
    C = C100 + mu * M - Ce
    C = np.round(C,decimals=14)

    # S is the stiffness matrix
    # units: GPa-1 = 10E9 m²/N
    S = np.linalg.inv(C)
    S = np.round(S,decimals=14)
    S[S == -0] = 0

    # E is the Young modulus
    # units: GPa = 10E9 N/m²
    E = [1/S[0,0],1/S[1,1],1/S[2,2]]
    E = np.array(E)


    # G is the shear modulus
    # units: GPa = 10E9 N/m²
    G = [1/S[3,3],1/S[4,4],1/S[5,5]]
    G = np.array(G)

    # unitless - between 0...1
    poisson = -S[:3,:3] / np.array([[S[0,0],S[0,0],S[0,0]],[S[1,1],S[1,1],S[1,1]],[S[2,2],S[2,2],S[2,2]]])

    # convert into other units 
    C = C / unitfac
    S = S * unitfac
    E = E / unitfac
    G = G / unitfac

    return C, S, E, G, poisson
