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
!   mat(:,:) -- Real(rdef): Input matrix array
!
! Output: 
!   mat_inv(:,:) -- Real(rdef): inverse of mat1
!-

!$Id$
!$Log$
!Revision 1.4  2002/02/23 20:32:19  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/16 21:04:18  helms
!Fixed problem with passing optional arguments.
!
!Revision 1.2  2001/09/27 18:31:53  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine mat_inverse (mat, mat_inv)

  use precision_def
  use nr

  implicit none

  real(rdef), intent(in) :: mat(:,:)
  real(rdef), intent(out) :: mat_inv(:,:)

  real(rdef), allocatable, save :: mat2(:,:)
  integer, allocatable, save :: indx(:)
  real(rdef) d

  integer n, i, j

!

  n = size (mat, 1)

  if (.not. allocated(indx)) then
    allocate (mat2(n,n), indx(n))
  elseif (size(indx) /= n) then
    deallocate (mat2, indx)
    allocate (mat2(n,n), indx(n))
  endif

  mat2 = mat

! This seems not to work
!  mat_inv = 0
  do i = 1, n
     do j = 1, n
        mat_inv(i,j) = 0
     enddo
  enddo

  forall (i = 1:n) mat_inv(i,i) = 1

  call ludcmp (mat2, indx, d)

  do i = 1, n
    call lubksb (mat2, indx, mat_inv(1:n,i))
  enddo

end subroutine
