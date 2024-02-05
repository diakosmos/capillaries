from math import *
from numpy import *
#
# 
#
# input:
#   array of R vals, Rar, and max radius, R0
#   array of z values, zar
#   array of initial concentrations, chiar, and 2D array for
#     resultant concentrations, chiar2d(nR,nz)
#   2-D array of diffusion constants verus R and z, Kar2d

R0 = 0.3e-6
K0 = 0.137 * 1e-4  # fiduciual diffusivity of ethanol in MKS units
chi0 = 1.0 / ( pi * R0**2)  # initial concentration normalized so integral dA = 1.0
nR = 20
nz = 100
Rar = array(myfrange(0.0,R0 * (2*nR - 1.0)/nR))
zar = array(myfrange(0.0,18e-2,nz))
Kar2d = ones([nR,nz]) * K0
chiar = ones(nR) * chi0
chiar2d = ones([nR,nz]) * chi0

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
    for i in range(nR):
        Rdualar[i] = sqrt(0.5*(Rar[i]**2 + Rar[i+1]**2))
    Rdualar[nR] = R0
    # Next, array of Hodge map between the two:
    Hodge0 = Rar * 0.0
    for i in range(1,nR):
        Hodge0[i] = pi * (Rdualar[i]**2 - Rdualar[i-1]**2)
    Hodge0[0]  = pi * (Rdualar[0]**2)
    # Hodge1 map (trivial):
    Hodge1 = Rdualar * pi
    # Some other arrays:
    chiar = Rar * 0.0
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
        arout[nR] = -arin[nR] / (R0 - Rar[nR])
    # timestep operator:
    def advance(Xout,Xin,F,dz):
        for i in range(1,nR):
            Xout[i] = Xin[i] + dz * (F[i]-F[i-1])
        Xout[0] = Xin[0] + dz * F[0]
    # Now initialize first column of chiar2d:
    chiar2d[:,0] = chiar
    #
    # Main loop:
    #
    for j in range(nz-1):
        chiar = chiar2d[:,j]
        Kar = Kar2d[:,j]
        Xar = Hodge0 * chiar 
        grad(gradchiar,chiar)
        far = - kar * gradchiar
        Far = Hodge1 * far
        dz = zar[j+1]-zar[j]
        advance(nextXar,Xar,Far,dz)
        chiar2d[:,j+1] = (1.0 / Hodge0) * nextXar


radialdiffuse(Rar,R0,zar,Kar2d,chiar,chiar2d)
