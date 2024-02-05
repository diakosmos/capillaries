from math import *
from myfuncs import *
from numpy import *
from pylab import *
#
# 
#
# input:
#   array of R vals, Rar, and max radius, R0
#   array of z values, zar
#   array of initial concentrations, chiar, and 2D array for
#     resultant concentrations, chiar2d(nR,nz)
#   2-D array of diffusion constants verus R and z, Kar2d

# K/u in (cm^2 / cm), i.e. cm
#cfl: (dR)**2/dz > K/u ---> dz < dR**2 / K/u
#
# K ethanol: 0.137 cm^2 /s
# u = 400m/s -> 800m/s in core;
# K/u = 1.7125e-8
#
nR = 10
R0 = 0.3e-3
L = 18e-2
K0 = 1.7125e-8 * 10 
dzcfl = (R0 / nR)**2 / K0
nz = int(1.5 * L/dzcfl)*2 # 1.5 is a safety factor

#nR = 10
#nz = 100

#R0 = 1.0
#K0 = 0.137 * 1e-4  # fiduciual diffusivity of ethanol in MKS units
#chi0 = 1.0 / ( pi * R0**2)  # initial concentration normalized so integral dA = 1.0
chi0 = 1.0


Rar = array(myfrange(0.0,R0 * (2*nR - 1.0)/(2.0*nR),nR))
zar = array(myfrange(0.0,L,nz))
totalions = zar * 0.0
Kar2d = ones([nR,nz]) * K0
chiar = ones(nR) * chi0
chiar2d = ones([nR,nz]) * chi0

# parabolic profile for Kar2d:
Kar1d = ones(nR) * K0
for i in range(nR):
    Kar1d[i] = K0 * 2 / (1 - Rar[i]**2 + 0.1)
for j in range(nz):
    Kar2d[:,j] = Kar1d

def radialdiffuse(Rar,R0,zar,Kar2d,chiar,chiar2d):
    nR = len(Rar)
    nz = len(zar)
    # check to make sure Rar[0] = 0.0:
    if (Rar[0] != 0.0):
        return
    # now create some arrays
    # First, array of dual node locations:
    Rdualar = Rar * 0.0
    # I don't think this is quite right but let's do this to start:
    for i in range(nR-1):
        Rdualar[i] = sqrt(0.5*(Rar[i]**2 + Rar[i+1]**2))
    Rdualar[nR-1] = R0
    # Next, array of Hodge map between the two:
    Hodge0 = Rar * 0.0
    for i in range(1,nR):
        Hodge0[i] = pi * (Rdualar[i]**2 - Rdualar[i-1]**2)
    Hodge0[0]  = pi * (Rdualar[0]**2)
    # Hodge1 map (trivial):
    Hodge1 = Rdualar * pi
    # Some other arrays:
    gradchiar = Rar * 0.0
    Xar = Rar * 0.0
    nextXar = Rar * 0.0
    Kar = Rar * 0.0
    far = Rar * 0.0
    Far = Rar * 0.0
    # Gradient operator:
    def grad(arout,arin):
        for i in range(nR-1):
            arout[i] = (arin[i+1]-arin[i])/(Rar[i+1]-Rar[i])
        arout[nR-1] = -arin[nR-1] / (R0 - Rar[nR-1])
    # timestep operator:
    def advance(Xout,Xin,F,dz):
        for i in range(1,nR):
            Xout[i] = Xin[i] + dz * (F[i-1]-F[i])
        Xout[0] = Xin[0] - dz * F[0]
    # Now initialize first column of chiar2d:
    chiar2d[:,0] = chiar
    #
    # Main loop:
    #
    for j in range(nz-1):
        chiar = chiar2d[:,j]
        #        print ("Here is chiar in step ", j, " : ", chiar)
        Kar = Kar2d[:,j]
        Xar = Hodge0 * chiar

        grad(gradchiar,chiar)
        far = - Kar * gradchiar
        Far = Hodge1 * far
        dz = zar[j+1]-zar[j]
        advance(nextXar,Xar,Far,dz)
        chiar2d[:,j+1] = (1.0 / Hodge0) * nextXar
        
    for j in range(nz):
        totalions[j] = sum(Hodge0 * chiar2d[:,j])

                                    
radialdiffuse(Rar,R0,zar,Kar2d,chiar,chiar2d)
