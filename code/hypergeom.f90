!+
! Function HYPERGEOM (HGCX, ARG)
!
! Function to calculate a particular hypergeometric function
! -- Created by Daniel Fromowitz, January 1999.
!
! Input:
!     HGCX -- Integer: Determines which hypergeometric series is used
!       If HGCX = 1, use series coefficients of F (3/4, 5/4; 1; )
!                                              2 1
!
!       If HGCX = 2, use series coefficients of F (5/4, 7/4; 2; ) * 3/4
!                                              2 1
!
!     ARG  -- Real: Last (fourth) argument of the hypergeometric function
!
! Output:
!     HYPERGEOM -- Real: The hypergeometric function
!-

  function hypergeom (hgcx, arg)
  implicit none

  integer hgcx, i
  real arg, arg_power, hypergeom, next_term
  include 'hypergeom.inc'

  if (hgcx == 1) then
    hypergeom = 1.0
  else
    hypergeom = 3/4.e0
  endif

  i = 1
  arg_power = arg
   10 next_term = hgc(hgcx,i) * arg_power
  if (next_term/hypergeom > 1.d-5) then
    hypergeom = hypergeom + next_term
    i = i + 1
    arg_power = arg_power * arg
    if (i <= i_max) goto 10
  endif

  return
  end
