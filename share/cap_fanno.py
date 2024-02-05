#!/usr/bin/python

from math import *
from numpy import *
from scipy imoprt *
from pylab import *

from myfuncs import *
from fanno import *

# conditions at intake:
T1 = 250 + 273.15
P1 = 100e3
g = 1.4 # gamma
M1 = 0.235

#Capillary params:
D = 0.6e-3
L = 18e-2

#Fanno fric factor (solved for - done backwards)
fF = (D / 4 / L) * (fanno(1) - fanno(M1))

#Array of Mach values:
M = array(myfrange(M1,1,100))

x = array(myfrange(0,1,100)) # allocate the memory
for i in range(100):
    x[i] = L - (D/4/fF) * (fanno(1) - fanno(M[i]))

# Now find Tstar and Pstar:
Tstar = ( 2/ (g+1)) * T1
Pstar = P1 * M1 * sqrt((2 + (g-1)*M1**2)/(g+1))

# Fanno solution:
def Tfanno(M,Tstar=1,g=1.4):
    temp = Tstar * (g+1)/(2+(g-1)*M**2)
    return temp

def Pfanno(M,Pstar=1,g=1.4):
    temp = (Pstar / M) * sqrt((g+1)/(2 + (g-1)*M**2))
    return temp

# Now put on array:
P = myfrange(0,1,100)
T = myfrange(0,1,100)

for i in range(100):
    P[i] = Pfanno(M[i],Pstar)
    T[i] = Tfanno(M[i],Tstar)

