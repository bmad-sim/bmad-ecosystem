!+
! Subroutine query_string (query_str, do_upcase, return_str, ix, ios)
!
! Subroutine to ask a question and accept a string. Leading blanks will
! always be trimmed.
!
! Input:
!   query_str   -- Character*(*): Query string.
!   do_upcase   -- Logical: If True then return_str will up converted
!                     to upper case.
!
! Output:
!   return_str  -- Character*(*): String typed in by the user.
!   ix          -- Integer: Number of character in first word
!   ios         -- Integer: I/O status (= 0 if not an error like ^Z typed).
!-

subroutine query_string (query_str, do_upcase, return_str, ix, ios)

use sim_utils, except => query_string
use precision_def

implicit none

character(*) return_str
character(*) query_str
integer ix, ios
logical do_upcase

integer traceflag
common / traceback / traceflag

!

write (*, '(1x, 2a)', advance = 'NO') trim(query_str), ' '
read (*, '(a)', iostat = ios) return_str

if (ios /= 0) then
  ix = 0
  return_str = ""
  return
endif

call string_trim (return_str, return_str, ix)

if (do_upcase) call str_upcase (return_str, return_str)

end subroutine

