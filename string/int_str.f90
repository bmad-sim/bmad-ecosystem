!+
! Function int_str(int) result (str)
!
! Routine to return a string representation of an integer.
!
! Input:
!   int    -- integer:
!
! Output:
!   str       -- character(:), allocatable: String representation.


function int_str(int) result (str)

use sim_utils_interface, dummy => int_str

implicit none

integer int
character(:), allocatable :: str
character(20) string


write (string, '(i0)') int
allocate (character(len_trim(adjustl(string))):: str)
str = trim(adjustl(string))

end function int_str

