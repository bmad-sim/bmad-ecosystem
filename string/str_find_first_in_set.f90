!+
! Function str_find_first_in_set (line, set)
!
! Function to locate the first character in "line" that matches
! a character in "set". If no match is found, function returns 0.
!
! Input:
!   line, set -- Character(*)
!
! Output:
!   ix_match -- integer: Index of the first matched character in "line"
!-

function str_find_first_in_set(line, set) result (ix_match)

implicit none
character(*) line, set
integer ix_match

! Find first character in "line" that is in "set"

do ix_match = 1, len(line)
   if (index(set,line(ix_match:ix_match)) /= 0) return
enddo

ix_match = 0 ! no match found

end function
