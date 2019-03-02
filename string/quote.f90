!+
! Function quote(str) result (q_str)
!
! Function to put double quote marks (") around a string.
! The string will be trimmed of trailing blanks.
!
! Input:
!   str     -- character(*): Input string
!
! Output:
!   q_str   -- character(:), allocatable: String with quote marks.
!-

function quote(str) result (q_str)
character(*) str
character(:), allocatable :: q_str
integer n

!

n = len_trim(str) + 2
allocate(character(n) :: q_str)
q_str = '"' // trim(str) // '"'

end function quote
