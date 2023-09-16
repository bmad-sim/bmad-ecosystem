!+
! Function max_nonzero (lbnd, array1, array2) result (ix_max)
!
! Routine to return the index ix_max such that all entries in array1 and array2
! with index greater than ix_max are zero. and either array1(ix_max) or array2(ix_max) 
! is non zero.
!
! Inputs:
!   lbnd        -- integer: Lower bound of array1 and array2
!   array1(:)   -- real(rp): Array of numbers.
!   array2(:)   -- real(rp), optional: Array of numbers. 
!                   If present, array2 must be the same size as array1
!
! Output:
!   ix_max      -- integer: index of largest nonzero entry in array1 or array2.
!                   If both arrays are zero, lbnd-1 is returned.
!-

function max_nonzero (lbnd, array1, array2) result (ix_max)

use precision_def

implicit none

integer lbnd, ix_max

real(rp) array1(lbnd:)
real(rp), optional :: array2(lbnd:)

!

if (present(array2)) then
  do ix_max = ubound(array1, 1), lbnd, -1
    if (array1(ix_max) /= 0) return
    if (array2(ix_max) /= 0) return
  enddo
else
  do ix_max = ubound(array1, 1), lbnd, -1
    if (array1(ix_max) /= 0) return
  enddo
endif

ix_max = lbnd - 1

end function
