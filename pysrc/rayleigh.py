#!/usr/bin/python

from math import *
from numpy import *
from scipy import *
from pylab import *

from myfuncs import *
from fanno import *

# conditions at intake:
T1 = 250 + 273.15
P1 = 100e3
g = 1.4 # gamma
Rair = 287.0

#Capillary params:
#D = 0.6e-3

nn = 100 # no. of nodes


def turb(Ren,turbstate):
    if (Ren > 4000):
        return True
    if (Ren < 2000):
        return False
    return turbstate

# friction factor:
def fFanno(Re,turbstate):
    if turbstate == False:
        fFanno = 16./Re
    else:
        fFanno = 0.00765 + 9.58/Re
    return fFanno

def fRe(Q,T,D):
    mu = 2.53e-5 * sqrt(T / 473.15)    
    return (4/pi) * Q / D / mu

# The idea is to start with a guess for M1 (Mach at inlet)
# and then solve for M2 (Mach at outlet), which should be 1....
# (so then iterate to find actual M1).

# vars: M, T, P, rho

x = array(myfrange(0,1,nn))
M = x * 0
T = x * 0
P = x * 0
fricint = x * 0
turbar = x * 0

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
    fF = fFanno(fRe(Q,T1,D),turbar[0])

    for i in range(1,nn):
        fricint[i] = fricint[i-1] + 4 * fF * dx / D 
        def fannorootfunc(M):
            return fanno(M) - fanno(M1) - fricint[i]
        M[i] = rtbis(fannorootfunc,M1,1.0,1e-5)
        T[i] = Tstar * (g+1)/(2 + (g-1)*M[i]**2)
        Ren = fRe(Q,T[i],D)
        turbar[i] = turb(Ren,turbar[i-1])
#        if (turbar == True):
#            print "turbulent."
        fF = fFanno(Ren,turbar[i])

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
        
        
