module fanno_mod

  use constants_mod
  use globals_mod

contains

  function fannointegral (M) result (fannint)
    ! returns integral of (4 f_F / D) dx, as given e.g. by
    ! eq. 3.97, p. 113 in Anderson's book on compressible flow.
    implicit none
    double precision, intent (in) :: M
    double precision :: fannint
    fannint = -(2.4 / 2.8) * log(M**2 / (1. + 0.2 * M**2)) - 1. / (1.4 * M**2)
    return
  end function fannointegral

  function flowtry (Qstp, D, L, T1, P1, x, n) result (Mout)
    implicit none
    ! integrates Fanno curve from inlet to outlet; result is M at outlet.
    ! If M<<1, then guess for Q was too low. If M > 1.0....
    double precision, intent (in) :: Qstp, D, L, T1, P1, fF
    double precision, dimension (n) :: x
    integer, intent (in) :: n
    double precision :: Mout

    double precision :: A, dx, Q, M1, Tstar, Pstar
    double precision :: pi
    pi = 4*atan2(1.0,1.0)

    A = (pi/4) * D**2
    dx = L / (n-1)
    Q = Qstp * 2.125e-5_8
    M1 = Q / P1 / sqrt(g/(Rair * T1)) / A
    turbar(1) = 1
    Tstar = T1 * (2 + (g-1)*M1**2)/(g+1)
    Pstar = P1 * M1 * sqrt((2.0 + (g-1)*M1**2)/(g+1))

    M(1) = M1
    T(1) = T1
    fricint(1) = 0.0
    Renar(1) = fRe(Q,T(1),D)
    fF = fFanno(fRe(Q,T1,D),turbar(1))
    do i=2,n
       fricint(i) = fricint(i-1) + 4 * fF * dx / D
    enddo
  end function flowtry

  function fRe(Q,T,D) result (Re)
    double precision :: Q,T,D,Re
    mu = 2.53e-5_8 * sqrt(T / 473.15_8)
    Re = (4 / pi) * Q / D / mu
    return
  end function fRe

  function fanning(Re,turb,method,acc) return (fF)
    double precision :: Re, acc
    integer :: turb, method
    double precision :: fF
    if (turb .eq. 0) then
       fF = 16. / Re
    else
       if (method .eq. 0) then
          fF = rtbis(prandtl,1.0e-8,1.0,acc)
       else
          fF = 0.00765 + 9.58/Re
       endif
    endif
    return
  end function fanning

  function prandtl(fF)
    double precision :: fF, prandtl
    prandtl = 4.0 * log10(Re * sqrt(fF)) - 0.4 - 1./sqrt(fF)
    return
  end function prandtl


end module fanno_mod
