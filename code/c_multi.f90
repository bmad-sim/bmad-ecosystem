!+
! Function c_multi (n, m)
!
! Subroutine to compute multipole factors: 
!          c_multi(n, m) =  +/- ("n choose m")/n! 
!
! Input:
!   n,m -- Integer: For n choose m
!
! Output:
!   c_multi  -- Real(rdef): Multipole factor
!-

#include "CESR_platform.inc"

function c_multi (n, m) result (c_out)

  use bmad

  implicit none

  integer, intent(in) :: n, m
  integer in, im
  
  real(rdef) c_out, factorial_n
  real(rdef), save :: c(n_pole_maxx, n_pole_maxx)

  logical, save :: init_needed = .true.

! The magnitude of c(n, m) is number of combinations normalized by n!

  if (init_needed) then

    c(0, 0) = 1

    do in = 1, n_pole_maxx
      c(in, 0) = 1
      c(in, in) = 1
      do im = 1, in-1
        c(in, im) = c(in-1, im-1) + c(in-1, im)
      enddo
    enddo

    factorial_n = 1

    do in = 0, n_pole_maxx
      if (in > 0) factorial_n = in * factorial_n
      do im = 0, in
        c(in, im) = c(in, im) / factorial_n
        if (mod(im, 4) == 0) c(in, im) = -c(in, im)
        if (mod(im, 4) == 3) c(in, im) = -c(in, im)
      enddo
    enddo

    init_needed = .false.

  endif

!

  c_out = c(n, m)

end function
