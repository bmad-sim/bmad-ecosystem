module sign_of_mod

use utilities_mod

implicit none

private sign_of_real, sign_of_int

!+
! Function sign_of (num, zero_is_zero) result (num_sign)
!
! Routine to return the sign of a number.
! Note: Fortran instrinsic sign function is similar to sign_of with zero_is_zero = False.
!
! Input:
!   num           -- integer or real(rp): Input number
!   zero_is_zero  -- logical, optional: If True (default), num = 0 gives num_sign = 0.
!                     If False, num = 0 gives num_sign = 1.
!
! Output:
!   num_sign      -- integer or real(rp): +1 if num is positive, -1 if num is negative,
!                     and 0 or +1 if num is zero depending upon setting of zer_is_zero.
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
! Function sign_of_real (num, zero_is_zero) result (num_sign)
!
! Overloaded routine. See documentation on sign_of function.
!-

elemental function sign_of_real (num, zero_is_zero) result (num_sign)

real(rp), intent(in) :: num
integer num_sign
logical, optional, intent(in) :: zero_is_zero

!

if (num == 0) then
  if (logic_option(.true., zero_is_zero)) then
    num_sign = 0
  else
    num_sign = 1
  endif
elseif (num < 0) then
  num_sign = -1
else
  num_sign = 1
endif

end function sign_of_real

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Function sign_of_int (num, zero_is_zero) result (num_sign)
!
! Overloaded routine. See documentation on sign_of function.
!-

elemental function sign_of_int (num, zero_is_zero) result (num_sign)

integer, intent(in) :: num
integer num_sign
logical, optional, intent(in) :: zero_is_zero

!

if (num == 0) then
  if (logic_option(.true., zero_is_zero)) then
    num_sign = 0
  else
    num_sign = 1
  endif
elseif (num < 0) then
  num_sign = -1
else
  num_sign = 1
endif

end function sign_of_int

end module
