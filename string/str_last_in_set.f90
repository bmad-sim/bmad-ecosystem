!+
! Function str_last_in_set (line, set) result (ix_match)
!
! Function to locate the last character in "line" that matches
! a character in "set". If no match is found, function returns 0.
!
! Input:
!   line, set -- Character(*)
!
! Output:
!   ix_match -- integer: Index of the last matched character in "line"
!-

function str_last_in_set(line, set) result (ix_match)

implicit none
character(*) line, set
integer ix_match

! Find last character in "line" that is in "set"

do ix_match = len(line), 1, -1
   if (index(set,line(ix_match:ix_match)) /= 0) return
enddo

ix_match = 0 ! no match found

end function
