from math import *
from vapor import *

dt = 1.0  # timestep (sec)

# ambient state:
tekinf = 273.15 + 50  # ambient temperature
prpinf = 0.01 * prpsat(tekinf)  # set relative humidity = 10%
#prpinf = 591.0 * 100.0

# drop initial state:
adrp = 20 * mks.micron  # drop radius
tekdrp = tekinf     # drop temp
tauinvd = 0.0       # drop evaporation rate

mass = mass(adrp)
time = 0.0

for i in range(30):
#while (adrp > 1 * mks.micron):
    tauinvd = tauinv(adrp, tekinf, tekdrp, prpinf)
    mass *= 1 + dt * tauinvd
    adrp = radius(mass)
    time += dt
    tekdrp = tdrop(adrp, tekinf, tauinvd)
    print "%3i  %8.3g  %8.3g  %8.3g  %8.3g  %8.3g  %8.3g  %8.3g" \
          % (i, adrp,tekdrp,tekinf,tauinvd,dt,dt*tauinvd, time)
    dt = 0.1/abs(tauinvd)
    
