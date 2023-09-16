!+
! Function str_first_not_in_set (line, set) result (ix_match)
!
! Function to locate the first character in "line" that
! does not match a character in "set". If all characters match,
! function returns 0.
!
! Input:
!   line, set -- Character(*)
!
! Output:
!   ix_match -- integer: Index of the first non-matched character in "line"
!-

function str_first_not_in_set(line, set) result (ix_match)

implicit none
character(*) line, set
integer ix_match

! return the first character in "set" not found in "line"

do ix_match = 1, len(line)
   if (index(set, line(ix_match:ix_match)) == 0) return
enddo

ix_match = 0 ! no match found

end function
