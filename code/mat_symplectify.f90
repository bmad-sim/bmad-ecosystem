!+
! Subroutine mat_symplectify (mat_in, mat_symp)
!
! Subroutine to form a symplectic matrix that is  "close" to the input matrix.
!
! Modules needed:
!   use bmad_interface
!
! Input:
!   mat_in(:,:) -- Real: Input matrix to symplectify.
!
! Output:
!   mat_symp(:,:) -- Real: Symplectic output matrix
!-

!$Id$
!$Log$
!Revision 1.4  2002/01/16 21:04:18  helms
!Fixed problem with passing optional arguments.
!
!Revision 1.3  2001/10/02 18:49:12  rwh24
!More compatibility updates; also added many explicit variable declarations.
!
!Revision 1.2  2001/09/27 18:31:53  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine mat_symplectify (mat_in, mat_symp)

  implicit none

  real, intent(in)  :: mat_in(:,:)
  real, intent(out) :: mat_symp(:,:)

  real, save, allocatable :: m1(:,:), m2(:,:), m3(:,:), m_symetric(:,:)

  integer n,i,j

!

  n = ubound(mat_in, 1)

  if (ubound(mat_symp, 1) /= n) then
    type *, 'ERROR IN MAT_SYMPLECTIFY: UNEQUAL MATRIX SIZES.'
    call err_exit
  endif

!

  if (.not. allocated (m1)) then
    allocate (m1(n,n), m2(n,n), m3(n,n), m_symetric(n,n))
  elseif (ubound(m1, 1) /= n) then
    deallocate (m1, m2, m3, m_symetric)
    allocate (m1(n,n), m2(n,n), m3(n,n), m_symetric(n,n))
  endif

! form the symmetrix matrix:
!         m1 = S * (1 - mat_in) * (1 + mat_in)^-1
!         m_symetric = (m1 + m1_transpose) / 2

  m1 =-mat_in
  forall (i = 1:n) m1(i,i) = m1(i,i) + 1

  m2 = mat_in
  forall (i = 1:n) m2(i,i) = m2(i,i) + 1
  call mat_inverse (m2, m2)

  m3 = matmul(m1, m2)

  do i = 2, n, 2
    m1(i-1,1:n) =  m3(i,  1:n) 
    m1(i,  1:n) = -m3(i-1,1:n) 
  enddo

  forall (i = 1:n, j = 1:n) m_symetric(i,j) = (m1(i,j) + m1(j,i)) / 2

! form the symplectic matrix
!         m_symp = (1 + S * m_symetric) * (1 - S * m_symetric)^-1

  do i = 2, n, 2
    m1(i-1,1:n) =  m_symetric(i,  1:n) 
    m1(i,  1:n) = -m_symetric(i-1,1:n) 
  enddo

  m2 = m1
  forall (i = 1:n) m2(i,i) = m2(i,i) + 1

  m3 = -m1
  forall (i = 1:n) m3(i,i) = m3(i,i) + 1
  call mat_inverse (m3, m3)
  
  mat_symp = matmul (m2, m3)

end subroutine
