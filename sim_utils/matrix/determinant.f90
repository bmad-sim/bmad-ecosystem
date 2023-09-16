!+
! Function determinant (mat) result (det)
!
! Routine to take the determinant of a square matrix
! This routine is adapted from Numerical Recipes.
!
! Input: 
!   mat(:,:) -- Real(rp): Square matrix.
!
! Output: 
!   det      -- Real(rp): Deteminant.
!-

function determinant(mat) result (det)

use super_recipes_mod, except => determinant
use precision_def

implicit none

real(rp) mat(:, :), det

integer i, ns
integer :: indx(size(mat, 1))

real(rp) :: mat1(size(mat, 1), size(mat, 1)) 

logical err

! Special case

ns = size(mat, 1)

if (ns == 2) then
  det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
  return
endif

! singular case

do i = 1, ns
  if (all(mat(i,:) == 0)) then
    det = 0
    return
  endif
enddo

! normal case

mat1 = mat
call super_ludcmp(mat1, indx, det, err)
if (err .and. global_com%exit_on_error) call err_exit

do i = 1, ns
  det = det * mat1(i, i)
enddo


end function
