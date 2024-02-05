#
# Constants for idealized (perfect) air
#

import constants

mu = 28.965 # effective molecular wt of air
Rs = constants.Rgas / mu
g = 1.4 # gamma
cp = Rs * g/(g-1)
cv = Rs / (g-1)
