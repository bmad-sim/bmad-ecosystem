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

complex(rp), parameter :: i_imaginary = (0.0d0, 1.0d0)
  
! real_garbage$ and int_garbage$ can be used, for example, to identify
! variable that have not been set.

integer, parameter :: int_garbage$ = -987654
real(rp), parameter :: real_garbage$ = -987654.3

! lf$ (the line feed or LF character) can be used to encode a multiline string.
! EG: string = 'First Line' // lf$ // 'Second Line'

character(1), parameter :: lf$ = achar(10)

! 

integer, parameter :: x_axis$ = 1, y_axis$ = 2, z_axis$ = 3

! True and false

real(rp), parameter :: true$ = 1, false$ = 0
integer, parameter :: true_int$ = 1, false_int$ = 0
integer, parameter :: yes$ = 1, no$ = 0, maybe$ = 2

! Color escape sequences

character(*), parameter :: rl_prompt_start_ignore = achar(1)   ! For use with GNU readline routine.
character(*), parameter :: rl_prompt_end_ignore = achar(2)     ! For use with GNU readline routine.

character(*), parameter :: black_color = achar(27) // '[30m' 
character(*), parameter :: red_color = achar(27) // '[31m' 
character(*), parameter :: green_color = achar(27) // '[32m' 
character(*), parameter :: yellow_color = achar(27) // '[33m' 
character(*), parameter :: blue_color = achar(27) // '[34m' 
character(*), parameter :: magenta_color = achar(27) // '[35m' 
character(*), parameter :: cyan_color = achar(27) // '[36m' 
character(*), parameter :: gray_color = achar(27) // '[37m' 


character(*), parameter :: dark_gray_color = achar(27) // '[90m' 
character(*), parameter :: peach_color = achar(27) // '[91m' 
character(*), parameter :: light_green_color = achar(27) // '[92m' 
character(*), parameter :: light_yellow_color = achar(27) // '[93m' 
character(*), parameter :: light_blue_color = achar(27) // '[94m' 
character(*), parameter :: pink_color = achar(27) // '[95m' 
character(*), parameter :: aqua_color = achar(27) // '[96m' 
character(*), parameter :: white_color = achar(27) // '[97m' 

character(*), parameter :: blink_color = achar(27) // '[5m' 
character(*), parameter :: bold_color = achar(27) // '[1m' 

character(*), parameter :: reset_color = achar(27) // '[0m' 

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
