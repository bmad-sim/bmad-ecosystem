!+
! Subroutine mat_inverse (mat, mat_inv)
!
! Program to take the inverse of a square matrix.
! program is from NUMERICAL RECIPES
!
! Modules needed:
!   bmad_interface
!
! Input: 
!   mat(:,:) -- Real: Input matrix array
!
! Output: 
!   mat_inv(:,:) -- Real: inverse of mat1
!-

subroutine mat_inverse (mat, mat_inv)

  use nr

  implicit none

  real, intent(in) :: mat(:,:)
  real, intent(out) :: mat_inv(:,:)

  real, allocatable, save :: mat2(:,:)
  integer, allocatable, save :: indx(:)
  real d

  integer n, i

!

  n = size (mat, 1)

  if (.not. allocated(indx)) then
    allocate (mat2(n,n), indx(n))
  elseif (size(indx) /= n) then
    deallocate (mat2, indx)
    allocate (mat2(n,n), indx(n))
  endif

  mat2 = mat
  mat_inv = 0
  forall (i = 1:n) mat_inv(i,i) = 1

  call ludcmp (mat2, indx, d)

  do i = 1, n
    call lubksb (mat2, indx, mat_inv(1:n,i))
  enddo

end subroutine
