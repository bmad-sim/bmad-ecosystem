!+
! Function is_integer (string, int, delims, ix_word) result (valid)
!
! Function to tell if the first word in a string is a valid integer.
!
! Leading blank characters in string are ignored.
! A blank string is considered to NOT be a valid integer.
!
! Input:
!   string        -- character(*): Character string.
!   delims        -- character(*), optional: List of end-of-word delimitors.
!
! Output:
!   int         -- integer, optional: Integer value. Default is -999 if not valid.
!   valid       -- logical: Is a valid integer.
!   ix_word     -- integer, optional: Index of first word. If not an integer then ix_word is 0.
!-

function is_integer (string, int, delims, ix_word) result (valid)

implicit none

character(*) string
integer i, i1
integer, optional :: int, ix_word
logical valid
character(*), optional :: delims

! ignore beginning spaces

valid = .false.
if (present(int)) int = -999
if (present(ix_word)) ix_word = 0

if (string == '') return

do i1 = 1, len(string)
  if (string(i1:i1) /= ' ') exit
enddo

! first char may be a + or -

if (index('+-', string(i1:i1)) /= 0) i1 = i1 + 1
if (i1 == len(string)+1) return

! need at least 1 digit

if (index('1234567890', string(i1:i1)) == 0) return
i1 = i1 + 1
valid = .true.

! check for a non-digit

do i = i1, len(string)
  if (string(i:i) == ' ') exit

  if (present(delims)) then
    if (index(delims, string(i:i)) /= 0) exit
  endif

  if (index('1234567890', string(i:i)) == 0) then
    valid = .false.
    return
  endif
enddo

if (i == len(string)+1) i = len(string)
if (present(int)) read(string(:i), *) int
if (valid .and. present(ix_word)) ix_word = i-1

end function
