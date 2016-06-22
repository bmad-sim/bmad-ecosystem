!+
! Subroutine downcase_string (string)
!
! Routine to convert a string to lowercase:
!
! Input:
!   string -- Character(*): Input string.
!
! Output:
!   string -- Character(*): String with all letters converted to lower case.
!-

subroutine downcase_string(string)

implicit none

integer string_len, i, ix

character(*) string
character(26), parameter :: lower_case = 'abcdefghijklmnopqrstuvwxyz'
character(26), parameter :: upper_case = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

!

string_len = len(string)

do i = 1, string_len
  ix = index(upper_case, string(i:i))
  if (ix /= 0) string(i:i) = lower_case(ix:ix)
enddo

end subroutine
