program main

!  use rtbis_mod
!  use funcs_mod
  use constants_mod
  use fanno_mod
  use globals_mod

  implicit none
  integer, parameter :: n = 1000

  ! global arrays
  !double precision, dimension (1:n) :: x, M, T, P, P0, rho, a, u, fricint, turbar, Renar
  double precision :: xx
  double precision :: pi
  pi = 4*atan2(1.0,1.0)
  allocate (x (n), M (n), T (n), P (n), P0 (n), rho (n), a (n), u (n), fricint (n), turbar (n), Renar (n))

  xx = fannointegral(0.5_8)
  print *, xx
  
  x(1) = 4.0_8
  

end program main
