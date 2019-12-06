!+
! Function real_str(r_num, n_signif) result (str)
!
! Routine to return a string representing a real number.
! Trailing zeros will be supressed. 
!
! Input:
!   r_num     -- real(rp): Real number
!   n_signif  -- integer: Number of significant places.
!
! Output:
!   str       -- character(:), allocatable: String representation with no tailing blanks or zeros.


function real_str(r_num, n_signif) result (str)

use sim_utils_interface, dummy => real_str

implicit none

real(rp) r_num
integer n_signif
character(:), allocatable :: str
character(20) string

string = real_to_string(r_num, 20, n_signif = n_signif)
allocate (character(len_trim(adjustl(string))):: str)
str = trim(adjustl(string))

end function real_str

