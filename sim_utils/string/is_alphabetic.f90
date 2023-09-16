!+
! Function is_alphabetic (string, valid_chars) result (is_alpha)
!
! Function to tell if a string has all alphabetical characters.
! Spaces are counted as not alphabetic
!
! Input:
!   string      -- character(*): Character string.
!   valid_chars -- character(*), optional: If present, other characters that may be present and
!                     are considered valid.
!
! Output:
!   is_alpha   -- logical: True if is alphabetic.
!-

function is_alphabetic (string, valid_chars) result (is_alpha)

implicit none

character(*) string
integer i
logical is_alpha
character(*), optional :: valid_chars
character(*), parameter :: abc = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'

is_alpha = .false.

do i = 1, len(string)
  if (index(abc, string(i:i)) /= 0) cycle
  if (present(valid_chars)) then
    if (index(valid_chars, string(i:i)) /= 0) cycle
  endif
  return
enddo

is_alpha = .true.

end function

