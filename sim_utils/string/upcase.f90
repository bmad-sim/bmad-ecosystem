!+
! Function upcase (str_in) result (str_out)
!
! Routine to convert a string to upper case.
!
! Moule needed:
!   use sim_utils_interface
!
! Input:
!   str_in  -- Character(*): Input string.
!
! Output:
!   str_out -- Character(*): String with all letters converted to upper case.
!-

elemental function upcase(str_in) result (str_out)

use sim_utils_interface, only: str_upcase

implicit none

character(*), intent(in) :: str_in
character(len(str_in)) str_out

!

call str_upcase (str_out, str_in)

end function
