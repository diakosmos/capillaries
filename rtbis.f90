module rtbis_mod

contains

  function rtbis(func,x1,x2,xacc) result (root)
    implicit none
    double precision, intent (in) :: x1, x2, xacc
    double precision :: root
    interface
       function func(x)
         double precision, intent (in) :: x
         double precision :: func
       end function func
    end interface

    integer, parameter :: maxit = 80
    integer :: j
    double precision :: dx, f, fmid, xmid

    fmid = func(x2)
    f = func(x1)
    if (f*fmid .ge. 0.0) print *, 'rtbis: root must be bracketed'
    if (f .le. 0.0) then
       root=x1
       dx=x2-x1
    else
       root=x2
       dx=x1-x2
    endif
    do j=1,maxit
       dx = dx * 0.5
       xmid = root+dx
       fmid=func(xmid)
       if (fmid .le. 0.0) root = xmid
       if (abs(dx) .lt. xacc .or. fmid .eq. 0.0) return
    enddo
    print *, 'too many iterations in rtbis'
  end function rtbis

end module rtbis_mod
