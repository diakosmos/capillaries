from math import *

def mfp(P,T,sigma=1.3e-25):
    # P in mTorr, T in K, sigma in cm^2; answer in cm
    r_eff = 0.5 * sqrt(sigma / pi)
    r_eff_A = r_eff / 1e10
    print r_eff_A
    ell = 23.0 * (1./P) * (T/300.) * (1/r_eff_A**2)
    return ell