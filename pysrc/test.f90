program profile

  real :: x
  real :: ree, accx
  integer :: flag

  ree = 1.0
  accx = 2.0
  flag = 1

  x = fanning(ree,flag,accx)
  

contains

function fanning (re,turb,acc) result (fF)
  real, intent (in) :: re, acc
  integer, intent (in) :: turb
  real, intent (out) :: fF

  fF = 4.3

end function fanning

end program profile
