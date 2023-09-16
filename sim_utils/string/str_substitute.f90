!+
! Subroutine str_substitute (string, str_match, str_replace, do_trim, ignore_escaped)
!
! Routine to substitute all instances of one sub-string for another in a string
!
! Input:
!   string          -- character(*): Character string.
!   str_match       -- character(*), optional: Sub-string to be replaced. 
!                        Default is the tab char.
!   str_replace     -- character(*): optional: String to substitute in. 
!                        Default is the space char.
!   do_trim         -- logical, optional: Default False. If true then substitute only
!                        in the region string(1:n) where n = len_trim(string).
!                        This argument only affects things if str_match is composed
!                        of blanks.
!   ignore_escaped  -- logical, optional: Default False. If true, ignore matches that are
!                        preceeded by the "\" escape character.
!
! Output:
!   string  -- character(*): String with all instances of str_match replaced for str_replace.
!-

subroutine str_substitute (string, str_match, str_replace, do_trim, ignore_escaped)

implicit none

character(*) string
character(*), optional :: str_match, str_replace

integer i, ixs, n_match
logical, optional :: do_trim, ignore_escaped
logical ignore

!

ignore = .false.
if (present(ignore_escaped)) ignore = ignore_escaped

n_match = 1
if (present(str_match)) n_match = len(str_match)

ixs = len(string)
if (present(do_trim)) then
  if (do_trim) ixs = len_trim(string)
endif
if(n_match.gt.1) ixs=ixs-(n_match-1)
i = ixs+1     !start looking where match interval just reaches end

do 
  i = i - 1
  if (i < 1) return

  if (present(str_match)) then
    if (string(i:i+n_match-1) /= str_match) cycle
  else
    if (string(i:i) /= char(9)) cycle
  endif

  if (ignore .and. i > 1) then
    if (string(i-1:i-1) == '\') cycle         ! '
  endif

  if (present(str_replace)) then
    string = string(:i-1) // str_replace // trim(string(i+n_match:))
    i = i - n_match + 1
  else
    string = string(:i-1) // ' ' // trim(string(i+n_match:))
  endif
enddo

end subroutine
