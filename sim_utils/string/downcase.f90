!+
! Function downcase (str_in) result (str_out)
!
! Routine to convert a string to lower case.
!
! Input:
!   str_in  -- Character(*): Input string.
!
! Output:
!   str_out -- Character(*): String with all letters converted to lower case.
!-

elemental function downcase(str_in) result (str_out)

use sim_utils_interface, only: str_downcase

implicit none

character(*), intent(in) :: str_in
character(len(str_in)) str_out

!

call str_downcase (str_out, str_in)

end function
