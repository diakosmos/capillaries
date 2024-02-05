from fanno4 import *

nT = 100

TCmin = 0.
TCmax = 350.

Tmin = TCmin + 273.15
Tmax = TCmax + 273.15

Tar = exp(array(myfrange(log(Tmin),log(Tmax),nT)))
TCar = Tar - 273.15

Qar = zeros(nT)

for Ti in range(nT):
    Qar[Ti] = flowsolve2(600.e-6,18.e-2,Tar[Ti])
            
