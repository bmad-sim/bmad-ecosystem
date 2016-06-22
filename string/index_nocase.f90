!+
! Function index_nocase (string, match_str) result (indx)
!
! Function to look for a sub-string of string that matches match_str.
! This routine is similar to the fortran INDEX function
! except it is case insensitive.
!
! Input:
!   string    -- Character(*): String to look for match.
!   match_str -- Character(*): String to match to.
!
! Output:
!   indx -- Integer:
!              = 0  -> no match.
!              > 0  -> index of where match starts.
!-

function index_nocase(string, match_str) result (indx)

implicit none

character(*) string, match_str
character(len(string)) string_1
character(len(match_str)) string_2
integer indx

!

call str_upcase(string_1, string)
call str_upcase(string_2, match_str)
indx = index(string_1, string_2)

end function
