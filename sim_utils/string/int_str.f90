!+
! Function int_str(int, width) result (str)
!
! Routine to return a string representation of an integer.
!
! Input:
!   int     -- integer: Number to write
!   width   -- integer, optional: If present, width of output field. 
!               In this case the field will be zero padded. EG: int_str(12, 4) -> "0012"
!
! Output:
!   str       -- character(:), allocatable: String representation.


function int_str(int, width) result (str)

use sim_utils_interface, dummy => int_str

implicit none

integer int
integer, optional :: width
character(:), allocatable :: str
character(20) string, fmt

!

if (present(width)) then
  allocate (character(width):: str)
  write (fmt, '(3(a, i0))') '(i', width, '.', width, ')'
  write (str, fmt) int

else
  write (string, '(i0)') int
  allocate (character(len_trim(adjustl(string))):: str)
  str = trim(adjustl(string))
endif

end function int_str

