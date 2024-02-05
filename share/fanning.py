from math import *
from myfuncs import *

# Returns Fanning friction factor
# using Prandtl's law, or my fit to
# Prandtl over 2000 < Re < 4000:
def fanning(Re,turb=True,method=0,acc=1e-6):
    if (turb==False):
        return 16./Re
    else:
        if (method==0):
            def prandtl(fF):
                return 4.0 * log10(Re * sqrt(fF)) - 0.4 - 1./sqrt(fF)
            fF = rtbis(prandtl,1.0e-8,1.0,acc)
            return fF
        else:
            return 0.00765 + 9.58/Re
