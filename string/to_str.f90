!+
! Function to_str(num, max_signif) result (string)
!
! Routine to return the string representation of a number.
!
! Input:
!   num         -- real(rp): Number to convert to a string.
!   max_signif  -- integer, optional: Maximum significant digits. Default is 14.
!
! Output:
!   string      -- character(:), allocatable: String representation.
!-

function to_str(num, max_signif) result (string)

use sim_utils, dummy => to_str

implicit none

real(rp) num
integer pl, n_signif
integer, optional :: max_signif
character(:), allocatable :: string
character(16) fmt

!

allocate (character(24) :: string)
n_signif = integer_option(14, max_signif)

if (num == 0) then
  string = '0'
  return
endif

pl = floor(log10(abs(num)))
if (abs(pl) >= 10) then
  fmt = '(2a, i0)'
elseif (pl < 0) then
  fmt = '(2a, i3.2)'
else
  fmt = '(2a, i2.2)'
endif

!

if (pl > 5 .or. n_signif < pl + 1) then
  write (string, fmt) trim(rchomp(num/10.0_rp**pl, 0, max_signif)), 'E+', pl

elseif (pl > -3) then
  string = rchomp(num, pl, max_signif)

else
  write (string, fmt) trim(rchomp(num*10.0_rp**(-pl), 0, max_signif)), 'E', pl
endif

string = trim(string)

!---------------------------------------------------
contains

function rchomp (rel, plc, max_signif) result (out)

implicit none

real(rp) rel
character(24) out
character(8) :: fmt = '(f24.xx)'
integer, optional :: max_signif
integer n_signif, it, plc, ix

!

n_signif = integer_option(14, max_signif)

write (fmt(6:7), '(i2.2)') n_signif - 1 - plc
write (out, fmt) rel
do it = len(out), 1, -1
  if (out(it:it) == ' ') cycle
  if (out(it:it) == '0') then
    out(it:it) = ' '
    cycle
  endif
  if (out(it:it) == '.') out(it:it) = ' '
  call string_trim(out, out, ix)
  return
enddo

end function rchomp

end function to_str
