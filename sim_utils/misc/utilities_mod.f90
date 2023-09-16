module utilities_mod

use precision_def

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function integer_option (integer_default, opt_integer)
!
! Function to return the value of an integer if it is present or a default value
! if it is not.
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

elemental function integer_option (integer_default, opt_integer) result (integer_out)

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

end function integer_option

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function logic_option (logic_default, opt_logic)
!
! Function to return True or False dependending upon the state of an 
! optional logical.
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

elemental function logic_option (logic_default, opt_logic) result (logic_out)

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

end function logic_option

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function real_option (real_default, opt_real)
!
! Function to return the value of a real if it is present or a default value
! if it is not.
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

elemental function real_option (real_default, opt_real) result (real_out)

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

end function real_option

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function string_option (string_default, opt_string) result (str_out)
!
! Function to return either opt_string if it is present, or string_default 
! if opt_string is not present.
!
! Input:
!   string_default -- Character(*): Default if opt_string is not present.
!   opt_string     -- Character(*), optional: Input string.
!
! Output:
!   string_out    -- Character(:), allocatable:
!                       = opt_string       if present.
!                       = string_default   otherwise.
!-

function string_option (string_default, opt_string) result (string_out)

implicit none

character(*) :: string_default
character(*), optional :: opt_string
character(:), allocatable :: string_out

!

if (present(opt_string)) then
  string_out = opt_string
else
  string_out = string_default
endif

end function string_option

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
  
end function eval_logical

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Function on_off_logic (logic, true_str, false_str) result (string)
!
! Function to encode the state of the logic argument in a string.
!
! Input:
!   logic     -- Logical: True or False.
!   true_str  -- character(*), optional: If present then use this instead of "ON".
!   false_str -- character(*), optional: If present then use this instead of "OFF"
!
! Output:
!   string  -- Character(8): "ON" or "OFF" except if true_str or false_str is present.
!-

function on_off_logic (logic, true_str, false_str) result (string)

logical logic
character(8) string
character(*), optional :: true_str, false_str

!

if (logic) then
  if (present(true_str)) then
    string = true_str
  else
    string = 'ON'
  endif

else
  if (present(false_str)) then
    string = false_str
  else
    string = 'OFF'
  endif
endif

end function on_off_logic

end module
