from math import *


def deltas(T,P):
    Rair = 287.0
    g = 1.4
    cp = Rair * g / (g-1)
    return cp*log(T/273.15) - Rair*log(P/100.e3)
