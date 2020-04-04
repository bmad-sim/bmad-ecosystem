!+
! Function is_alphabetic (string) result (is_alpha)
!
! Function to tell if a string has all alphabetical characters.
! Spaces are counted as not alphabetic
!
! Input:
!   string  -- Character(*): Character string.
!
! Output:
!   is_alpha   -- Logical: True if is alphabetic 
!-

function is_alphabetic (string) result (is_alpha)

implicit none

character(*) string
integer i
logical is_alpha
character(*), parameter :: abc = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

is_alpha = .false.
do i = 1, len(string)
  if (index(abc, string(i:i)) == 0) return
enddo
is_alpha = .true.

end function

