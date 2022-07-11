!+
! Function real_str(r_num, n_signif, n_decimal) result (str)
!
! Routine to return a string representing a real number.
! Trailing zeros will be supressed.
!
! Also see: real_to_string
!
! Input:
!   r_num     -- real(rp): Real number
!   n_signif  -- integer: Number of significant places.
!   n_decimal -- integer, optional: If present, maximum number of places after the decimal point.
!
! Output:
!   str       -- character(:), allocatable: String representation with no tailing blanks or zeros.


function real_str(r_num, n_signif, n_decimal) result (str)

use sim_utils_interface, dummy => real_str

implicit none

real(rp) r_num
integer n_signif
integer, optional :: n_decimal
character(:), allocatable :: str
character(20) string

if (present(n_decimal)) then
  if (abs(r_num) > 10.0_rp**(n_signif - n_decimal)) then
    write (string, '(f20.' // int_str(n_decimal) // ')') r_num
  else
    string = real_to_string(r_num, 20, n_signif = n_signif)
  endif
else
  string = real_to_string(r_num, 20, n_signif = n_signif)
endif

allocate (character(len_trim(adjustl(string))):: str)
str = trim(adjustl(string))

end function real_str

