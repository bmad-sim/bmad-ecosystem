!+
! Function real_to_string (num, fmt) result (str)
! 
! Routine to turn a real number into a string for printing.
!
! Input:
!   num -- real(rp): Input number.
!   fmt -- characher(*), optional: Currently not used.
!
! Output:
!   str -- character(24): String representation of num.
!-

function real_to_string (num, fmt) result (str)

use precision_def

implicit none

real(rp) num
integer pl
character(24) str
character(24) fmt2
character(*), optional :: fmt

if (num == 0) then
  str = '0'
  return
endif

pl = floor(log10(abs(num)))

if (pl > 5) then
  fmt2 = '(2a, i0)'
  write (str, fmt2) trim(rchomp(num/10.0**pl, 0)), 'E', pl

elseif (pl > -3) then
  str = rchomp(num, pl)

else
  fmt2 = '(2a, i0)'
  write (str, fmt2) trim(rchomp(num*10.0**(-pl), 0)), 'E', pl

endif

!--------------------------
contains

function rchomp (rel, plc) result (out)

implicit none

real(rp) rel
character(24) out
character(8) :: fmt = '(f24.xx)'
integer it, plc, ix

!

write (fmt(6:7), '(i2.2)') 10-plc
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

end function

end function
