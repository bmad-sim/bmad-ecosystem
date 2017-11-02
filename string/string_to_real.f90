!+
! Function string_to_real (line, default, err_flag, err_print_flag) result (value)
!
! Function to convert a string to an real.
!
! Input:
!   line            -- character*(*): String to decode
!   default         -- real: Default value to use for a blank string.
!   err_print_flag  -- logical, optional: If present and False then suppress error message printing.
!
! Output:
!   value    -- real: Variable to hold value
!   err_flag -- logical: Set .true. if there is a decoding error
!-

function string_to_real (line, default, err_flag, err_print_flag) result (value)

use output_mod, dummy => string_to_real

implicit none

real(rp) default, value
integer ix, ios
character(*) line
logical err_flag
logical, optional :: err_print_flag
character(*), parameter :: r_name = 'string_to_real'

!

err_flag = .false.

call string_trim (line, line, ix)
if (ix .eq. 0) then 
  value = default
else
  read (line, *, iostat = ios) value
  if (ios /= 0) then
    err_flag = .true.
    if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, 'ERROR DECODING NUMBER: ' // trim(line))
  endif
endif

end function
