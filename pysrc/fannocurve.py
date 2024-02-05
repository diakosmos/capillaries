#
# Purpose: calculate Fanno curve on T-s diagram: Given M, find T/T* and s-s*. 
#

from math import *
from air import *   # constants for perfect-gas air: Rs, cp, cv, etc

# returns ratio of T to T*:
def fannocurve_T(M):
    return (g+1)/(2+(g-1)*M**2)

# returns difference s-s*:
def fannocurve_s(M):
    return Rs * (0.5 * (g+1)/(g-1) * log(fannocurve_T(M)) + log(M))

