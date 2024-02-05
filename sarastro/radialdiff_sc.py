from myimports import *
import fanno3

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
nR = 8
R0 = 0.3e-3
L = 18e-2
Kethanol = 1.37e-5
K = Kethanol

chi0 = 1.0

Rar = array(myfrange(0.0,R0 * (2*nR - 1.0)/(2.0*nR),nR))
zar = fanno3.x
nz = size(zar)
u1d = fanno3.u
parabolicprofile = 2. * (R0**2-Rar**2) / R0**2
totalions = zar * 0.0
Kar2d = ones([nR,nz])
chiar = ones(nR) * chi0
chiar2d = ones([nR,nz]) * chi0
u2d = multiply.outer(parabolicprofile,u1d)
Kar2d = K / u2d
maxKeff = K / u2d.min()
maxdz = max(diff(zar))
dzcfl = (R0 / nR)**2 / maxKeff
if (dzcfl < maxdz):
    print "Reduce number of radial nodes or increase number of axial nodes for stability."

# 2D effective diffusivity array:
for i in range(nR):
    for j in range(nz):
        Kar2d[i,j] = K / u2d[i,j]

# Electric potential array (just allocate for now):
V = ones([nR,nz])

def radialdiffuse(Rar,R0,zar,Kar2d,chiar,chiar2d,V):
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
