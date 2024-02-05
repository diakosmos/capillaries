#
# units: cgs
#
################################################################################
# parameters:

N = 100                                 # no. of zones
(xmin, xmax) = (-2.0, 2.0)              #
deltax = (xmax - xmin)/(N+1)

Amin = 0.1                             # min area (cm^2) of nozzle
xL = 1.0                                # lengthscale for nozzle area vs. x

################################################################################
# imports:

from myfuncs import *

################################################################################
# imports:

################################################################################
# allocations:

# set or allocate memory for x, A(x)
x = myfrange(xmin,xmax,N)
A = [0.0] * N

# allocate memory for fundamental quantities: ( _p: predicted, _c: corrected)
rho = rho_p = rho_c = [0.0] * N
tmp = tmp_p = tmp_c = [0.0] * N
v_x = v_x_p = v_x_c = [0.0] * N

# allocate memory for derived quantities:
pre = [0.0] * N

################################################################################
# set:

for i in range(N):
    A[i] = Amin + (x[i] / xL)**2






