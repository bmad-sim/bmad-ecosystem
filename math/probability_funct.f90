!+
! Function probability_funct (x)
!
! Routine compute the normal probability function PROBABILITY_FUNCT.
!
! Defining relation:
!     probability_funct(x) = integral_{0, x} dt exp(-t^2/2) / sqrt(2*pi)
!
! Properties:
!     probability_funct(0)      = 0
!     probability_funct(\infty) = 0.5
!     probability_funct(-x)     = -probability_funct(x)
!
! In terms of the error function:
!     probability_funct(x) = erf(x/sqrt(2)) / 2
!-

!$Id$
!$Log$
!Revision 1.4  2003/07/09 01:29:31  dcs
!new bmad
!
!Revision 1.3  2002/02/23 20:34:50  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.2  2001/09/27 17:47:10  rwh24
!UNIX compatibility updates
!


function probability_funct(x)

  use nr


  use precision_def

  implicit none

  real(rp) probability_funct, x

  probability_funct = erf(x/sqrt(2.0)) / 2

end function
