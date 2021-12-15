!+
! Function real_to_string (real_num, width, n_signif, n_decimal) result (str)
! 
! Routine to turn a real number into a string for printing. 
! Printing the number without an explicit exponent is preferred.
! 
! See also:
!   reals_to_string
!   reals_to_table_row
!   real_str
!
! Notes:
!   Keep in mind: With floating format the largest number of digits needed for the exponent is 5. EG: "-1.2e-123".
!   Thus with decimal point and possible negative sign, the field width should be at least n_signif + 7
!   If the field width is too small to display the number of significant digits wanted, "******" will be printed.
!
! Examples:
!   real_to_string(1.20001234e-4_rp, 12, 3)      => "1.20E-4     "
!   real_to_string(1.20001234e-2_rp, 12, 6)      => "0.0120001   "
!   real_to_string(1.20001234e-4_rp, 12, 7)      => "1.200012E-04"
!   real_to_string(1.20001234e-4_rp, 12, 10)     => "************"
!   real_to_string(1.20001234e-4_rp, 12, 2, 6)   => "    0.000120"
!   real_to_string(1.20001234e-8_rp, 12, 2, 6)   => "    0.000000"
!   real_to_string(1.20001234e+6_rp, 12, 2, 6)   => "     1.2E+06"
!
! Input:
!   real_num    -- real(rp): Real number to encode.
!   width       -- integer: width of number field. The output string length will also be equal to width.
!   n_signif    -- integer, optional: Nominal number of significant digits to display. Must be non-negative. Default is 15.
!                    If n_decimal >= 0: n_signif is ignored if the number can be encoded with fixed format.
!                           n_signif only affects encoding if the number is encoded with a floating point format.
!                    If n_decimal < 0: Number of significant digits displayed will be reduced if trailing zeros are suppressed.
!                    
!   n_decimal   -- integer, optional: 
!                     n_decimal = positive or 0: Number of digits after the decimal point to display if fixed format is used. 
!                         Fixed format is preferred and floating format will only be used if necessary.
!                         Numbers are right justified. This is good for tables.
!                     n_decimal = -1 (default): Trailing zeros replaced by blanks. Numbers are left justified.
!                         The more compact fixed/floating format will be used. This is good for numbers in text.
!
! Output:
!   str         -- character(width): String representation of real_num. The length will be equal to the value of the width argument.
!-

function real_to_string (real_num, width, n_signif, n_decimal) result (str)

use sim_utils, dummy => real_to_string

implicit none

real(rp) real_num, rnum, cut

integer width
integer, optional :: n_signif, n_decimal
integer ix, n_sig, n_dec, n_exp, n_digits, width_net, n_neg, e_width

character(width) str
character(30) fmt, str1, str2

!

rnum = real_num
n_sig = integer_option(15, n_signif)
n_dec = integer_option(-1, n_decimal)
if (n_sig > 17) n_sig = 17

if (n_sig <= 0) then
  str = 'ERROR!'
  return
endif

if (real_num == 0) then
  if (n_dec == -1) then
    str = '0'
  elseif (width < n_dec+2) then
    call x_string(str)
  else
    write (fmt, '(a, i0, a, i0, a)') '(f', width, '.', n_dec, ')'
    write (str, fmt) real_num
  endif
  return
endif

if (real_num < 0) then
  n_neg = 1
  width_net = width - 1     
else
  n_neg = 0
  width_net = width
endif

n_exp = floor( log10(abs(real_num)) )   ! real_num = a * 10^n_exp where a is in the range 1 <= a < 10
! EG: A real_num like 0.997 with n_sig = 2 is treated as 1.0 due to rounding.
cut = 10.0_rp**(n_exp+1) * (1 - 0.5_rp * 10.0_rp**(-n_sig))
if (n_sig > 0 .and. abs(real_num) > cut) then
  n_exp = n_exp + 1
  rnum = 10.0_rp**n_exp
endif

if (n_exp >= 100) then
  e_width = 5   ! "E+nnn"
elseif (n_exp >= 10) then
  e_width = 4   ! "E+nn"
elseif (n_exp >= 1) then
  e_width = 4   ! "E+nn"
elseif (n_exp == 0) then
  e_width = 0   ! No exponent needed
elseif (n_exp >= -9) then
  e_width = 4   ! "E-nn"
elseif (n_exp >= -99) then
  e_width = 4   ! "E-nn"
else
  e_width = 5   ! "E-nnn"
endif

!----------------
! Compact format wanted

if (n_dec == -1) then
  ! Construct floating point string
  call rchomp(rnum * 10.0_rp**(-n_exp), n_sig-1, str2)
  write (str2, '(2a, sp, i0.2)') trim(str2), 'E', n_exp

  if (n_exp > 6) then
    if (len_trim(str2) > width) then
      call x_string(str)
    else
      str = str2
    endif
    return
  endif

  ! Construct fixed point string
  if (n_exp+1 >= n_sig) then
    write (str1, '(i0)') nint(real_num)
  else
    call rchomp(rnum, n_sig-n_exp-1, str1)
  endif
  
  ! Choose the more compact
  if (min(len_trim(str1), len_trim(str2)) > width) then
    call x_string(str)
  elseif (len_trim(str2) < len_trim(str1)) then
    str = str2
  else
    str = str1
  endif
  return
endif

!----------------

! Case: Fixed format integer (no displayed decimal point).
! Example: real_to_string(9812.345, 12, 3, 0) => "       9812"

if (n_dec == 0 .and. n_exp+1 <= width_net) then
  write (fmt, '(a, i0, a)') '(f', width+1, '.0)'
  write (str1, '(f30.0)') real_num  
  if (str1(30:30) == '.') str1(30:30) = ' '
  str = adjustl(str1)
  str = adjustr(str)

! Case: Fixed format with decimal point

elseif (n_dec > 0 .and. n_exp+1+n_dec <= 17 .and. n_exp+2+n_dec <= width_net) then
  write (fmt, '(a, i0, a, i0, a)') '(f', width, '.', n_dec, ')'
  write (str, fmt) real_num

! Case: Floating format

else
  n_sig = min(n_sig-1, width_net-e_width-2)
  if (n_sig > 0) then
    if (n_exp >= 100) then
      write (fmt, '(a, i0, a, i0, a)') '(es', width, '.', n_sig, 'E3)'
    else
      write (fmt, '(a, i0, a, i0, a)') '(es', width, '.', n_sig, ')'
    endif
    write (str, fmt) real_num
  else
    call x_string(str)
  endif
endif

!--------------------------
contains

subroutine x_string(str)
character(*) str
character(30), parameter :: x_str = 'X*****************************'
str = x_str(1:len(str))
end subroutine

!--------------------------
! contains

subroutine rchomp (num, n_pl, str_out)

implicit none

real(rp) num
integer n_pl, it, ix
character(*) str_out
character(30) line, fmt

!

write (fmt, '(a, i0, a)') '(f25.', n_pl, ')'
write (line, fmt) num
do it = len(line), 1, -1
  if (line(it:it) == ' ') cycle
  if (line(it:it) == '0') then
    line(it:it) = ' '
    cycle
  endif
  if (line(it:it) == '.') line(it:it) = ' '
  call string_trim(line, str_out, ix)
  return
enddo

end subroutine

end function
