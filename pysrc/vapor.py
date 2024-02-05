from math import *
import mks

# material constants:
# Dnu = 0.239 cm^2/sec, 67thCRC, p. F-45
Dnu = 0.239 * 1e-4  # diffusion coefficient for water vapor in air
Leh2o = 2.257e6     # specific enthalpy of vaporization of water, J/kg
kgair = 0.024       # thermal cond. of air, in W/m/K, www.engineeringtoolbox.com
mh2o = 18.0153 * mks.u * mks.NA # molar mass of water

def hPa2torr(phpa):
    # hPa in torr: P(torr) = 0.75006 * P(hPa)
    return phpa * 0.75006

def Pa2torr(ppa):
    return ppa * 0.0075006

def torr2hPa(prt):
    return prt / 0.75006

def torr2Pa(prt):
    return prt / 0.0075006

def mfp(P,T,sigma=1.3e-19):
    # mfp = 1 / (n sigma)
    # Note one mole per 22.4L at STP; L = m^3 / 1000
    # Note that default sigma=1.3e-19 comes from r_eff = 0.5 Angstrom,
    # which is approx value for He, and sigma = pi*(2*r)^2.
    # Effective cross section for air or N2 is higher.
    n = mks.NA / (22.4 * 0.001) * (P / mks.atm) * (273.15 / T)
    mfp = 1.0 / n / sigma
    return mfp

def prpsat(tek):
    """ Returns saturation water vapor pressure in Pa, prp,
    (Goff-Gratch equation) as a function of temp in Kelvin, tek.
    Taken from wikipedia entry for Goff-Gratch_equation, 2007/8/28"""
    # Formula in K and hPa. Routine converts hPa to Pa at end.
    c1 = -7.90298
    c2 = 5.02808
    c3 = -1.3816e-7
    c4 = 8.1328e-3
    e3 = 11.344
    e4 = -3.49149
    psthpa = 1013.25 # saturation pressure at steam point temperature, in hPa
    tekstp = 373.15  # steam point temperature, in K
    lgphpa = c1 * (tekstp / tek - 1) + c2 * log10(tekstp/tek)
    lgphpa += c3*(10**(e3*(1-tek/tekstp))-1)
    lgphpa += c4*(10**(e4*(tekstp/tek - 1))-1) + log10(psthpa)
    phpa = 10**lgphpa
    prp = 100.0 * phpa
    return prp

def RRy(q, srften=0.072):
    """Rayleigh stability limit for electrostatic explosion,
    correspoinding to energy of surface = coulomb energy.
    See p.20 of book on ESI/MS ed. by Richard Cole.
    Note charge q is charge number, i.e. q=1,2,3..."""
    qq = q*mks.ee
    R3 = qq**2 / 64 / pi**2 / mks.epsilon0 / srften
    R = pow(R3,1.0/3.0)
    return R

def Dnustar(Dnu):  # not coded up yet - for now let Dvstar = Dv.
    """ see Pruppacher and Klett 1997, p.506 & 509, for
    effective diffusion constant when Kn becomes of order unity"""
    return Dnu

#def dmdt(a, tekinf, tekdrp, prpinf):
#    """ Returns dm/dt of a water droplet in kg/s, given a (radius of particle
#    in meters), tekinf (ambient temperature in K), tekdrp (temp of
#    droplet in K), and prpinf (ambient pressure OF H2O!!!! in Pa) """
#    prpsatdrp = prpsat(tekdrp)
#    temp = (prpinf / tekinf - prpsatdrp / tekdrp)
#    temp *= 4 * pi * a * Dnustar(Dnu) * mh2o / mks.R
#    return temp

def tauinv(a, tekinf, tekdrp, prpinf):
    # 1/tau, i.e. evaporation/condensation rate.
    prpsatdrp = prpsat(tekdrp)
    temp = ((prpinf / tekinf) - (prpsatdrp / tekdrp))
    temp *= 3 * Dnustar(Dnu) * mh2o / (mks.rhoh2o * a * mks.R)
    return temp

#def tdrop(a, tekinf, dmdtd):
#    """ Temperature of evaporating droplet, warmed by
#    ambient air and cooled by heat of vaporization.
#    From Moyle, Smidansky & Lamb (2006?). kgair is
#    effective thermal conductivity, which needs to be
#    corrected for high Kn effects, but I'll do that later."""
#    return (tekinf + (Leh2o / (4 * pi * kgair)) * dmdtd)

def tdrop(a, tekinf, tauinvd):
    """ Temperature of evaporating droplet, warmed by
    ambient air and cooled by heat of vaporization.
    From Moyle, Smidansky & Lamb (2006?). kgair is
    effective thermal conductivity, which needs to be
    corrected for high Kn effects, but I'll do that later."""
    return (tekinf + (mks.rhoh2o * a**2 * Leh2o /(3 * kgair)) * tauinvd)

def mass(a):
    return (4.0 * pi / 3.0) * pow(a,3) * mks.rhoh2o

def radius(m):
    return pow(3.0*m /(4.0 * pi * mks.rhoh2o),1.0/3.0)
