!+
! sim_utils_struct
!-

module sim_utils_struct

use precision_def
use physical_constants

implicit none

! A all_pointer_struct is just a pointer to either a real, integer, or logical variable.
! This is used to construct arrays of pointers.

type all_pointer_struct
  real(rp), pointer :: r => null()
  integer, pointer :: i => null()
  logical, pointer :: l => null()
end type 

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function is_true (param) result (this_true)
!
! Routine to translate from a real number to a boolian True or False.
! Translation: 0 = False, nonzero = True
!
! The typical use of this routine is for parameters in ele_struct%value(:) which
! is a real array. Some of the elements in the %value array are used to specify
! boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).
! 
! Input:
!   param    -- real(rp): Real number to be translated
!
! Output:
!   this_true -- logical: Set False if param is zero. True otherwise.
!-

function is_true (param) result (this_true)

real(rp) param
logical this_true

!

this_true = (param /= 0)

end function is_true

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function is_false (param) result (this_false)
!
! Routine to translate from a real number to a boolian True or False.
! Translation: 0 = False, nonzero = True
!
! The typical use of this routine is for parameters in ele_struct%value(:) which
! is a real array. Some of the elements in the %value array are used to specify
! boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).
! 
! Input:
!   param    -- real(rp): Real number to be translated
!
! Output:
!   this_false -- logical: Set True if param is zero. False otherwise.
!-

function is_false (param) result (this_false)

real(rp) param
logical this_false

!

this_false = (param == 0)

end function is_false

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Function value_of_all_ptr (a_ptr) result (value)
!
! Routine to return the value pointed to by an all_pointer_struct.
!
! Input:
!   a_ptr -- all_pointer_struct: Pointer to a variable
!
! Output:
!   value -- real(rp): Value. Set to real_garbage$ if number of pointer components of a_ptr that
!             are associated is not 1.
!-

function value_of_all_ptr (a_ptr) result (value)

type (all_pointer_struct) a_ptr
real(rp) value
integer n

!


n = 0

if (associated(a_ptr%r)) then
  n = n + 1
  value = a_ptr%r
endif

if (associated(a_ptr%i)) then
  n = n + 1
  value = a_ptr%r
endif

if (associated(a_ptr%l)) then
  n = n + 1
  if (a_ptr%l) then
    value = true$
  else
    value = false$
  endif
endif

if (n /= 1) then
  value = real_garbage$
endif

end function value_of_all_ptr

end module
