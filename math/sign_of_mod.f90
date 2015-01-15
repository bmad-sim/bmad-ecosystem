module sign_of_mod

use precision_def

implicit none

private sign_of_real, sign_of_int

!+
! Function sign_of (num) result (num_sign)
!
! Routine to return the sign of a number.
! Similar to the sign intrinsic function except here there is only one argument.
!
! Input:
!   num  -- integer or real(rp): Input number
!
! Output:
!   num_sign -- integer: +1 if num is non-negative, -1 if num is negative.
!-

interface sign_of
  module procedure sign_of_real
  module procedure sign_of_int
end interface

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function sign_of_real (num) result (num_sign)
!
! See documentation on sign_of function.
!-

function sign_of_real (num) result (num_sign)

real(rp) num
integer num_sign

!

if (num < 0) then
  num_sign = -1
else
  num_sign = 1
endif

end function sign_of_real

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function sign_of_int (num) result (num_sign)
!
! See documentation on sign_of function.
!-

function sign_of_int (num) result (num_sign)

integer num
integer num_sign

!

if (num < 0) then
  num_sign = -1
else
  num_sign = 1
endif

end function sign_of_int

end module
