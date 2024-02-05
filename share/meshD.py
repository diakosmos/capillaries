from fanno2 import *

nD = 100

Dmin = 100.0e-6
Dmax = 1000.0e-6

Dar = exp(array(myfrange(log(Dmin),log(Dmax),nD)))

Qar = zeros(nD)

for Di in range(nD):
    Qar[Di] = flowsolve2(Dar[Di],18.e-2)
            
