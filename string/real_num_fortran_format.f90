!+
! Function real_num_fortran_format (number, width, n_blanks) result (fmt_str)
!
! Routine to find a "nice" edit descriptor format for a real number.
! The format will either be "es" for numbers that are very small or very large or "f".
!
! Also see: parse_fortran_format
!
! Example output:
!     "es12.4"
!     "f10.3"
!
! Input:
!   number    -- real(rp): Number to find format for
!   width     -- integer: Field width including leading blanks.
!   n_blanks  -- integer, optional: Number of leading blank spaces. Default is zero
!
! Output:
!   fmt_str   -- character(9): Format string.
!-

function real_num_fortran_format (number, width, n_blanks) result (fmt_str)

use sim_utils, dummy => real_num_fortran_format

implicit none

real(rp) number, num
integer, optional :: n_blanks
integer width, num_width, nb, nlog
character(9) str, fmt_str

! Differentiate between number = 0 and number = -0

nb = integer_option(0, n_blanks)

if (number == 0) then
  write (str, '(f8.2)') number
  if (index(str, '-') == 0) then
    write (fmt_str, '(a, i0, a, i0)') 'f', width, '.', width-2
  else
    write (fmt_str, '(a, i0, a, i0)') 'f', width, '.', width-3
  endif
  return
endif

!

num = abs(number)
num_width = width - nb 
if (number < 0) num_width = num_width - 1  ! Width excluding negative sign

nlog = log10(num) + 1d-14

if (nlog <= -99) then
  write (fmt_str, '(a, i0, a, i0, a)') 'es', width, '.', num_width-7, 'e3'

elseif (nlog < -3) then
  write (fmt_str, '(a, i0, a, i0)') 'es', width, '.', num_width-6

elseif (-3 <= nlog .and. nlog  < 0) then
  write (fmt_str, '(a, i0, a, i0)') 'f', width, '.', num_width-2

elseif (num_width-nlog-2 > -1) then
  write (fmt_str, '(a, i0, a, i0)') 'f', width, '.', num_width-nlog-2

elseif (nlog < 99) then
  write (fmt_str, '(a, i0, a, i0)') 'es', width, '.', num_width-6

else
  write (fmt_str, '(a, i0, a, i0, a)') 'es', width, '.', num_width-7, 'e3'
endif

end function
