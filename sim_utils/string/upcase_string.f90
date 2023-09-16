!+
! Subroutine upcase_string (string)
!
! Routine to convert a string to uppercase:
!
! Input:
!   string -- Character(*): Input string.
!
! Output:
!   string -- Character(*): String with all letters converted to uppercase.
!-

subroutine upcase_string(string)

implicit none

integer string_len, i, ix

character(*) string
character(26), parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
character(26), parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

!

string_len = len(string)

do i = 1, string_len
  ix = index(lower_case, string(i:i))
  if (ix /= 0) string(i:i) = upper_case(ix:ix)
enddo

end subroutine
