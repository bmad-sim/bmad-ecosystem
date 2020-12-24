!+
! Function cross_product (a, b) result (c)
! 
! Function to return the cross product of two vectors.
!
! Input:
!   a(3), b(3) -- Real(rp): Input vectors.
!
! Output:
!   c(3) -- Real(rp): Cross product: a X b.
!-

function cross_product (a, b) result (c)

use precision_def

implicit none

real(rp) a(:), b(:), c(3)

!

c = [a(2)*b(3) - a(3)*b(2), &
     a(3)*b(1) - a(1)*b(3), &
     a(1)*b(2) - a(2)*b(1)]

end function
