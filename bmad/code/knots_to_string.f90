function knots_to_string (x_knot, y_knot) result (str)

use sim_utils

implicit none

real(rp) x_knot(:), y_knot(:)
integer ik
character(:), allocatable :: str

!

allocate (character(1) :: str)
str = ''

do ik = 1, size(x_knot)
  str = str // ', (' // real_str(x_knot(ik),12) // ', ' // real_str(y_knot(ik),12) // ')'
enddo

str = str(3:)

end function knots_to_string

