module utilities_mod

use precision_def

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function integer_option (integer_default, opt_integer)
!
! Function to retrun True or False dependending upon the state of an 
! optional integer.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   integer_default -- Integer: Default if opt_integer is not present.
!   opt_integer     -- Integer, optional: Input integer.
!
! Output:
!   integer_option -- Integer: 
!                       = opt_integer       if present.
!                       = integer_default   otherwise.
!-

function integer_option (integer_default, opt_integer) result (integer_out)

  implicit none

  integer, intent(in) :: integer_default
  integer, intent(in), optional :: opt_integer
  integer integer_out

!

  if (present(opt_integer)) then
    integer_out = opt_integer
  else
    integer_out = integer_default
  endif

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function logic_option (logic_default, opt_logic)
!
! Function to retrun True or False dependending upon the state of an 
! optional logical.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   logic_default -- Logical: Default if opt_logic is not present.
!   opt_logic     -- Logical, optional: Input logical.
!
! Output:
!   logic_option -- Logical: 
!                       = opt_logic       if present.
!                       = logic_default   otherwise.
!-

function logic_option (logic_default, opt_logic) result (logic_out)

  implicit none

  logical, intent(in) :: logic_default
  logical, intent(in), optional :: opt_logic
  logical logic_out

!

  if (present(opt_logic)) then
    logic_out = opt_logic
  else
    logic_out = logic_default
  endif

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function real_option (real_default, opt_real)
!
! Function to retrun True or False dependending upon the state of an 
! optional real.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   real_default -- Real(rp): Default if opt_real is not present.
!   opt_real     -- Real(rp), optional: Input real.
!
! Output:
!   real_option -- Real(rp): 
!                       = opt_real       if present.
!                       = real_default   otherwise.
!-

function real_option (real_default, opt_real) result (real_out)

  implicit none

  real(rp), intent(in) :: real_default
  real(rp), intent(in), optional :: opt_real
  real(rp) real_out

!

  if (present(opt_real)) then
    real_out = opt_real
  else
    real_out = real_default
  endif

end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine string_option (string_default, opt_string)
!
! Subroutine to retrun True or False dependending upon the state of an 
! optional string.
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   string_default -- Character(*): Default if opt_string is not present.
!   opt_string     -- Character(*), optional: Input string.
!
! Output:
!   string_option -- Character(*): 
!                       = opt_string       if present.
!                       = string_default   otherwise.
!-

subroutine string_option (string_out, string_default, opt_string) 

  implicit none

  character(*), intent(in) :: string_default
  character(*), intent(in), optional :: opt_string
  character(*) string_out

!

  if (present(opt_string)) then
    string_out = opt_string
  else
    string_out = string_default
  endif

end subroutine

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function eval_logical (word, err) result (this_logic)
!
! Function of convert a string into a logical value.
! Accepted possibilities are:
!   .TRUE.  .FALSE. 
!    TRUE    FALSE
!    T       F
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   word  -- Character(*): Input string.
!   err   -- Logical: True if there is an error. False otherwise.
!
! Output:
!   this_logic -- Logical: Result.
!-

function eval_logical (word, err) result (this_logic)

implicit none

character(*), intent(in) :: word
character(len(word)+8) :: wd
logical this_logic
logical, intent(out) :: err
integer ix

!

call str_upcase(wd, word)
call string_trim (wd, wd, ix)
err = .false.

if (wd(1:ix) == '.TRUE.' .or. wd(1:ix) == 'TRUE' .or. wd(1:ix) == 'T') then
  this_logic = .true.
elseif (wd(1:ix) == '.FALSE.' .or. wd(1:ix) == 'FALSE' .or. wd(1:ix) == 'F') then
  this_logic = .false.
else
  err = .true.
endif
  
end function

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function on_off_logic (logic) result (name)
!
! Function to return the string "ON" or "OFF".
!
! Modules needed:
!   use cesr_utils
!
! Input:
!   logic -- Logical: True or False.
!
! Output:
!   name  -- Character(3): "ON" or "OFF"
!-

function on_off_logic (logic) result (name)

  logical logic
  character(3) name

!

  if (logic) then
    name = 'ON'
  else
    name = 'OFF'
  endif

end function

end module
