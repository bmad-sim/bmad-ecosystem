!+
! Function str_count(str, match) result (num)
!
! Count number of times match occurs in str. 
! Overlapping matches are ignored.
!
! Input:
!   str           -- character(*): String to match to.
!   match         -- character(*): Match string.
!
! Output:
!   num           -- integer: Number of matches.
!-

function str_count(str, match) result (num)

character(*) str, match
integer num, ix, ii

!

num = 0
ix = 1

do
  ii = index(str(ix:), match)
  if (ii == 0) return
  num = num + 1
  ix = ix + ii + len(match) - 1
enddo

end function
