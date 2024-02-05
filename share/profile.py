from myimports import *

from fanning import fanning

def profile(rar,D,Re,turb=True,u0=1.):
    # returns velocity (u) profile of laminar or turbulent flow in pipe,
    # where rar is array of R-values, D=pipe diameter, Re = Reynolds no,
    # turb = turbulence state (True or False, i.e. turbulent or laminar),
    # and u0 is the mean velocity, defined as Q = (pi/4)*D**2 * rho * u0,
    # assuming rho is independent of R. u is array of velocities output.
    u = zeros(size(rar))
    fF = fanning(Re,turb)
    if (turb == True):
        xhat = (1. - fanning(Re,False)/fanning(Re,True))**0.25
    else:
        xhat = 0.0
    factor = 0.125 * Re * fanning(Re,turb)
    for i in range(size(rar)):
        if (abs(2*rar[i]) >= xhat):
            u[i] = factor * (1. - 4. * rar[i]**2)
        else:
            u[i] = factor * (1. - xhat**2)
    return u
                
