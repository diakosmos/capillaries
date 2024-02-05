from fanno2 import *

nD = 100
nL = 100

Dmin = 100.0e-6
Dmax = 1000.0e-6

Lmin = 1.0e-2
Lmax = 30.0e-2

Dar = exp(array(myfrange(log(Dmin),log(Dmax),nD)))
Lar = exp(array(myfrange(log(Lmin),log(Lmax),nL)))

Qar = zeros([nD,nL])

for Di in range(nD):
    for Li in range(nL):
        Qar[Li,Di] = flowsolve2(Dar[Di],Lar[Li])
            
