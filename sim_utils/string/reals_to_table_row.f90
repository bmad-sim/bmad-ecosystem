!+
! Function reals_to_table_row (real_arr, width, n_decimal, n_blank) result (str)
! 
! Routine to turn an array of real numbers into a string for printing tables.
! Fixed format is preferred and floating format will only be used if necessary.
! 
! See also:
!   real_to_string
!   reals_to_string
!
! Notes:
!   Keep in mind: With floating format the largest number of digits needed for the exponent is 5. EG: "-1.2e-123".
!   Thus with decimal point and possible negative sign, the field width should be at least n_signif + 7
!   If the field width is too small to display the number of significant digits wanted, "******" will be printed.
!
!
! Input:
!   real_arr(:) -- real(rp): Array of real numbers.
!   width       -- integer: width of number field including blanks. 
!                     The output string will have length size(real_arr) * width.
!                     Minimum width is n_decimal+n_blank+3 if real_arr can have negative numbers.
!                     Making width larger will increase the crossover point between when numbers are encoded as 
!                     floating point and when numbers are encoded as fixed point.
!   n_decimal   -- integer: Number of digits after the decimal point to display if fixed format is used. 
!   n_blank     -- integer, optional: Minimum number of blank spaces between numbers. Default is 2.
!
! Output:
!   str         -- character(width): String representation of real_arr. The length of str = size(real_arr) * width.
!-

function reals_to_table_row (real_arr, width, n_decimal, n_blank) result (str)

use sim_utils, dummy => reals_to_table_row

implicit none

real(rp) real_arr(:)

integer width, n_decimal
integer, optional :: n_blank
integer n_signif, i, nb

character(width*size(real_arr)) str

!

nb = integer_option(2, n_blank)
n_signif = width - nb - 7

str = ''

do i = 1, size(real_arr)
   str((i-1)*width+1+nb:i*width) = real_to_string(real_arr(i), width-nb, n_signif, n_decimal)
enddo

end function
