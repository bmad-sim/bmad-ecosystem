!+
! Function string_to_int (line, default, err_flag, err_print_flag) result (value)
!
! Function to convert a string to an integer.
!
! Input:
!   line            -- character*(*): String to decode
!   default         -- integer: Default value to use for a blank string or if there is an error.
!   err_print_flag  -- logical, optional: If present and False then suppress error message printing.
!
! Output:
!   line     -- character(*): String without integer substring and with leading blanks trimmed.
!   value    -- integer: Variable to hold value.
!   err_flag -- logical: Set .true. if there is a decoding error. A blank string is not an error.
!-

function string_to_int (line, default, err_flag, err_print_flag) result (value)

use output_mod, dummy => string_to_int

implicit none

integer default, value, ix, ios
character(*) line
logical err_flag
logical, optional :: err_print_flag
character(*), parameter :: r_name = 'string_to_int'

!

err_flag = .false.
value = default

call string_trim (line, line, ix)
if (ix == 0) return

read (line, *, iostat = ios) value
if (ios /= 0) then
  err_flag = .true.
  if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, 'ERROR DECODING NUMBER: ' // trim(line))
endif

call string_trim(line(ix+1:), line, ix)

end function
