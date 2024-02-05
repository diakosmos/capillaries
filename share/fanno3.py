#!/usr/bin/python
from myfuncs import *
from math import *
from numpy import *
from scipy import *
import matplotlib
#matplotlib.use('wx')  # otherwise crashes. Think you need to switch to "ps" or "pdf" for vector graphics for publication work.
from pylab import *


from fanno import *

# conditions at intake:
#T1 = 250 + 273.15
T1 = 63. + 273.15
P1 = 100e3

# capillary parameters:
D = 600.e-6
L = 18.e-2

# gas parameters (don't change)
g = 1.4 # gamma
Rair = 287.0

nn = 10000 # no. of nodes

# Global variables
Tstar = 0.0
Pstar = 0.0

# Global arrays
x = array(myfrange(0,L,nn))
M = x * 0
T = x * 0
T0 = x * 0
P = x * 0
P0 = x * 0
rho = x * 0
a = x * 0
u = x * 0
fricint = x * 0
turbar = x * 0
Renar = x * 0

def fannofunc(M):
    temp =  - 1. / (1.4 * M**2) - (2.4 / 2.8) * log(M**2 / (1 + 0.2 * M**2))
    return temp

def turb(Ren,turbstate):
    if (Ren > 4000):
        return True
    if (Ren < 2000):
        return False
    return turbstate

def fFanno(Re,turbstate):
    if turbstate == False:
        fFanno = 16./Re
    else:
        fFanno = 0.00765 + 9.58/Re
        # this is a linear fit in 1/Re over range 2000-4000.
        # It would be better to fit to 1/sqrt(Re). Or, use
        # Prantdl's law, but that would be slower (iterative soln).
    return fFanno

# Reynolds number
def fRe(Q,T,D):
    mu = 2.53e-5 * sqrt(T / 473.15)    
    return (4/pi) * Q / D / mu

# The idea is to start with a guess for M1 (Mach at inlet)
# and then solve for M2 (Mach at outlet), which should be 1....
# (so then iterate to find actual M1).

# vars: M, T, P, rho



def flowsolve(QstpLm,D=1.0e-3,L=18e-2):
    A = (pi/4) * D**2
    dx = L / (nn-1)
    x = array(myfrange(0,L,nn))
#    myjmax = 40
#   myxacc = 1e-5
    Q = QstpLm * 2.125e-5
    M1 = Q / P1 / sqrt(g/(Rair * T1)) / A
    turbar[0] = 1 # entry always turbulent
    
    Tstar = T1 * (2 + (g-1)*M1**2)/(g+1)
    Pstar = P1 * M1 * sqrt((2 + (g-1)*M1**2)/(g+1))
    
    M[0] = M1
    T[0] = T1
#    P[0] = P1
    fricint[0] = 0
    Renar[0] = fRe(Q,T[0],D)
    fF = fFanno(fRe(Q,T1,D),turbar[0])

    for i in range(1,nn):
        fricint[i] = fricint[i-1] + 4 * fF * dx / D 
        def fannorootfunc(M):
            return fanno(M) - fanno(M1) - fricint[i]
        M[i] = rtbis(fannorootfunc,M1,1.0,1e-5)
        T[i] = Tstar * (g+1)/(2 + (g-1)*M[i]**2)
        Renar[i] = fRe(Q,T[i],D)
        turbar[i] = turb(Renar[i],turbar[i-1])
#        if (turbar == True):
#            print "turbulent."
        fF = fFanno(Renar[i],turbar[i])

#    P = (Pstar / M) * sqrt((g+1)/(2+(g-1)*M**2))

def flowsolve2(D,L=18e-2):
    def lfunc(q):
        flowsolve(q,D,L)
        return (-M[nn-1] + 0.99)
    ljmax = 40
    lxacc = 1e-5
    lx1 = 0
    lx2 = 100 # assume flow is between 0 and 100 L/min
    lfmid = lfunc(lx2)
    lf = lfunc(lx1)
    lrtbis = lx2
    ldx = lx1 - lx2
    for lj in range(1,ljmax):
        ldx = ldx * 0.5
        lxmid = lrtbis + ldx
        lfmid = lfunc(lxmid)
        if (lfmid <= 0):
            lrtbis = lxmid
        if ((abs(ldx) < lxacc) or (lfmid == 0)):
            return lrtbis
            break
        
        
Qstp = flowsolve2(D,L)
flowsolve(Qstp,D,L)
T0 = T * (1 + ((g-1)/2)*M**2)   # should be constant!
M1 = M[0]
Tstar = T1 * (2 + (g-1)*M1**2)/(g+1)
Pstar = P1 * M1 * sqrt((2 + (g-1)*M1**2)/(g+1))
P0star = Pstar / .528
P = (Pstar / M) * sqrt((g+1)/(2+(g-1)*M**2)) 
P0 = P0star * (1./M) * ((2+(g-1)*M**2)/(2+(g-1)))**((g+1)/(2*(g-1)))
rho = P / Rair / T
a = sqrt(g * Rair * T)
u = a * M

