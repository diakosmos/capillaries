from math import *
from numpy import *

# returns integral of (4 f_F / D) dx, as given e.g. by
# eq. 3.97, p. 113 in Anderson's book on compressible flow.
def fanno(M):
    temp = M**2
    temp /= 1 + 0.2 * M**2
    temp = log(temp)
    temp *= 2.4 / 2.8
    temp += 1 / (1.4 * M**2)
    temp *= -1
    return temp


