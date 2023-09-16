!+
! Function outer_product (a, b) result (c)
!
! Routine to form the outer product of two arrays.
!
! Input:
!   a(:), b(:)    -- real(rp): Two arrays.
!
! Output:
!   c(:,:)        -- real(rp): Outer product.
!-

function outer_product (a, b) result (c)

use precision_def

implicit none

real(rp) a(:), b(:)
real(rp) c(size(a), size(b))

integer i

!

do i = 1, size(a)
  c(i,:) = a(i) * b
enddo

end function
