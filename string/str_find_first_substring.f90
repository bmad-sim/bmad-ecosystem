!+
! Function str_find_first_substring(line, ix_match, substring) result (is_match)
!
! Function to return whether the "substring" is in "line".  
! The index of where the substring occurs in line is 
! returned in ix_match. If no match is found, ix_match is 0.
!
! Modules needed:
!   
!
! Input:
!   line, substring -- Character(*)
!
! Output:
!   ix_match -- integer: Index of the first matched character in "line"
!   is_match -- logical: returns whether a match was found
!-

function str_find_first_substring(line, ix_match, substring) result (is_match)

  implicit none
  character(*) line, substring
  integer ix_match
  logical is_match

  print *, line
  print *, substring
  ix_match = index(line, substring)

  if (ix_match <= 0) then ! no match found
    is_match = .false.
  else
    is_match = .true.
  end if


end function str_find_first_substring
