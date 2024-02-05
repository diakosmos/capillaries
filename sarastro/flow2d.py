from fanno3b import *
fannoflow(D,L,T1)
Qstp = flowsolve2(D,L,T1)
flowsolve(Qstp,D,L,T1)
T0 = T1 * (1 + ((g-1)/2)*M**2)   # should be constant!
M1 = M[0]
Tstar = T1 * (2 + (g-1)*M1**2)/(g+1)
Pstar = P1 * M1 * sqrt((2 + (g-1)*M1**2)/(g+1))
P0star = Pstar / .528
P = (Pstar / M) * sqrt((g+1)/(2+(g-1)*M**2)) 
P0 = P0star * (1./M) * ((2+(g-1)*M**2)/(2+(g-1)))**((g+1)/(2*(g-1)))
rho = P / Rair / T
a = sqrt(g * Rair * T)
u = a * M
