!+
! Subroutine mat_det (mat, det)
!
!=====================================================================
!=====================================================================
!============= DO NOT USE! USE THE ROUTINE DETERMINANT INSTEAD =======
!=====================================================================
!=====================================================================
!
! Subroutine to take the determinant of a square matrix
! This subroutine is adapted from NUMERICAL RECIPES
!
! Modules needed:
!   use sim_utils
!
! Input: 
!   mat(:,:) -- Real(rp): Square matrix.
!
! Output: 
!   det      -- Real(rp): Deteminant.
!-

#include "CESR_platform.inc"

subroutine mat_det(mat, det)

  use nr
  use precision_def

  implicit none

  integer i, ns
  integer, allocatable, save :: indx(:)

  real(rp) mat(:, :), det
  real(rp), allocatable, save :: mat1(:,:) 

! Special case

  ns = size(mat, 1)

  if (ns == 2) then
    det = mat(1,1) * mat(2,2) - mat(1,2) * mat(2,1)
    return
  endif

! Allocate space if needed

  if (allocated(indx)) then
    if (ns /= size(indx)) then
      deallocate (indx, mat1)
      allocate (indx(ns), mat1(ns, ns))
    endif
  else
    allocate (indx(ns), mat1(ns, ns))
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
  call ludcmp(mat1, indx, det)

  do i = 1, ns
    det = det * mat1(i, i)
  enddo

end subroutine
