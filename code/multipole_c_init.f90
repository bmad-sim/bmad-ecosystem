!+
! Subroutine multipole_c_init (c, maxx)
!
! Subroutine to compute multipole factors: 
!          c(n, m) =  +/- ("n choose m")/n! 
!
! Input:
!   maxx -- Integer: Size of c matrix
!
! Output:
!   c(0:maxx, 0:maxx) -- Real: Multipole factors.  
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"



subroutine multipole_c_init (c, maxx)

  implicit none

  integer n, m, maxx

  real c(0:maxx, 0:maxx), factorial_n

! The magnitude of c(n, m) is number of combinations normalized by n!

  c(0, 0) = 1

  do n = 1, maxx
    c(n, 0) = 1
    c(n, n) = 1
    do m = 1, n-1
      c(n, m) = c(n-1, m-1) + c(n-1, m)
    enddo
  enddo

  factorial_n = 1
  do n = 0, maxx
    if (n > 0) factorial_n = n * factorial_n
    do m = 0, n
      c(n, m) = c(n, m) / factorial_n
      if (mod(m, 4) == 0) c(n, m) = -c(n, m)
      if (mod(m, 4) == 3) c(n, m) = -c(n, m)
    enddo
  enddo

  return
  end
