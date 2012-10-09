!+
! Subroutine mat_symplectify (mat_in, mat_symp, f_scale)
!
! Subroutine to form the symplectic matrix that is as "close" as possible to the input matrix.
! If mat_in is symplectic, and f_length is 1, then mat_sym will be equal to mat_in.
! 
! The f_scale argument is used to "scale" mat_symp such that if mat_in is symplectic the
! following is true:
!   f_scale = 0   --> mat_symp = I (unit matrix, always true independent of mat_in)
!   f_scale = 1/N --> mat_symp^N = mat_in
!   f_scale = 1   --> mat_symp = mat_in
!   
! Module needed:
!   use sim_utils_interface
!
! Input:
!   mat_in(:,:) -- Real(rp): Input matrix to symplectify.
!   f_scale     -- Real(rp), optional: Scaling factor. Default is 1.
!
! Output:
!   mat_symp(:,:) -- Real(rp): Symplectic output matrix
!-

subroutine mat_symplectify (mat_in, mat_symp, f_scale)

use output_mod, except => mat_symplectify

implicit none

real(rp), intent(in)  :: mat_in(:,:)
real(rp), intent(out) :: mat_symp(:,:)
real(rp), intent(in), optional :: f_scale

real(rp), dimension(size(mat_in, 1), size(mat_in, 2)) :: m1, m2, m3, m_symmetric

integer n,i,j

!

n = ubound(mat_in, 1)

if (ubound(mat_symp, 1) /= n) then
  print *, 'ERROR IN MAT_SYMPLECTIFY: UNEQUAL MATRIX SIZES.'
  if (global_com%exit_on_error) call err_exit
endif

! Form the symmetrix matrix:
!         m1 = S * (1 - mat_in) * (1 + mat_in)^-1
!         m_symmetric = (m1 + m1_transpose) / 2

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

forall (i = 1:n, j = 1:n) m_symmetric(i,j) = (m1(i,j) + m1(j,i)) / 2

! Scale by f_scale

if (present(f_scale)) m_symmetric = f_scale * m_symmetric

! Form the symplectic matrix
!         m_symp = (1 + S * m_symmetric) * (1 - S * m_symmetric)^-1

do i = 2, n, 2
  m1(i-1,1:n) =  m_symmetric(i,  1:n)
  m1(i,  1:n) = -m_symmetric(i-1,1:n)
enddo

m2 = m1
forall (i = 1:n) m2(i,i) = m2(i,i) + 1

m3 = -m1
forall (i = 1:n) m3(i,i) = m3(i,i) + 1
call mat_inverse (m3, m3)

mat_symp = matmul (m2, m3)

end subroutine
