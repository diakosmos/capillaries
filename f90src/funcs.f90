module funcs_mod

contains

  function myfunc(x) result (f)
    implicit none
    real, intent (in) :: x
    real :: f
    f = cos(x)
    print *, x, "   ", f
    return
  end function myfunc

end module funcs_mod
