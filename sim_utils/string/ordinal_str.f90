!+
! Function ordinal_str(n) result (str)
!
! Routine to return a string representing the ordinal position of n.
! EG n = 1 -> "1st", n = 2 -> "2nd", etc.
! If n is negative, a negative sign prefix is added.
!
! Input:
!   n       -- integer: Cardinal number to convert to ordinal.
!
! Output:
!   str     -- character(:), allocatable: Ordinal string.
!-

function ordinal_str(n) result (str)

use sim_utils_interface, dummy => ordinal_str

implicit none

integer n
character(:), allocatable :: str
character(2) suffix

!

select case (modulo(abs(n), 10))
case (1);     suffix = 'st'
case (2);     suffix = 'nd'
case (3);     suffix = 'rd'
case default; suffix = 'th'
end select

allocate(character(4):: str)
str = int_str(n) // suffix

end function
