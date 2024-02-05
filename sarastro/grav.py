#!/usr/bin/python

from math import *
from numpy import *
from scipy import *
from pylab import *

from myfuncs import *


# conditions at intake:
T1 = 1.0
g = 1.4 # gamma
#Rair = 287.0

#Capillary params:
#D = 0.6e-3

nn = 100 # no. of nodes

def 


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

def energysolve(delphinorm):
    # solve energy equation to find Mach number given normalized delta phi
    def func(M):
        exp1 = (2 * g - 2) / (3 - g)
        exp2 = 4.0 / (3- g)
        temp = M**exp1 + (g-1)*M**exp2 + delphinorm - 0.5*(g+1)
        return temp
    M = r
        
        
