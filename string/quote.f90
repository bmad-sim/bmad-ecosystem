!+
! Function quote(str) result (q_str)
!
! Function to put double quote marks (") around a string.
! The string will not be modified if there are already quote marks of either kind.
! The string will always be trimmed of trailing blanks (but not leading blanks).
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

n = len_trim(str)
if (n > 0) then
  if ((str(1:1) == "'" .and. str(n:n) == "'") .or. (str(1:1) == '"' .and. str(n:n) == '"')) then
    allocate(character(n) :: q_str)
    q_str = trim(str)
    return
  endif
endif

allocate(character(n+2) :: q_str)
q_str = '"' // trim(str) // '"'

end function quote
