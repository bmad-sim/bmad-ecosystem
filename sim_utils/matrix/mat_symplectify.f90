!+
! Subroutine mat_symplectify (mat_in, mat_symp, p0_ratio, r_root)
!
! Subroutine to form the symplectic matrix that is as "close" as possible to the input matrix.
! If mat_in is symplectic, and r_root is 1, then mat_sym = mat_in.
! 
! When using coordinates like (x, Px/P0, y, Py/P0, z, dP/P0), where P0 is 
! the reference momentum, then:
!   p0_ratio = P0(exit coords) / P0(entrance coords). 
! Normally this is 1 except when dealing with something like an accelerating RF cavity.
! In the case where p0_ratio is not 1, the procedure is:
!   1) scale mat_in to a coordinate system where p0_ratio is 1.
!   2) symplecterize matrix.
!   3) scale back to a coordinate system with the given p0_ratio.
! Note: This also is usful with (-t, E) longitudinal coordinates
!
! The r_root argument is used to "scale" mat_symp such that if mat_in is symplectic the
! following is true:
!   r_root = 0   --> mat_symp = I (unit matrix, always true independent of mat_in)
!   r_root = 1/N --> mat_symp^N ~ mat_in
!   r_root = 1   --> mat_symp ~ mat_in
!   
! Input:
!   mat_in(:,:) -- Real(rp): Input matrix to symplectify.
!   p0_ratio    -- Real(rp): optional: Ratio of p0_exit / p0_entrance. Default is 1.
!                            If present, the size of mat_in must be 6.
!   r_root      -- Real(rp), optional: Scaling factor. Default is 1.
!
! Output:
!   mat_symp(:,:) -- Real(rp): Symplectic output matrix
!-

subroutine mat_symplectify (mat_in, mat_symp, p0_ratio, r_root)

use output_mod, except => mat_symplectify

implicit none

real(rp), intent(in)  :: mat_in(:,:)
real(rp), intent(out) :: mat_symp(:,:)
real(rp), intent(in), optional :: p0_ratio, r_root

real(rp), dimension(size(mat_in, 1), size(mat_in, 2)) :: m1, m2, m3, m_symmetric

integer n,i,j

!

n = ubound(mat_in, 1)

if (ubound(mat_symp, 2) /= n) then
  print *, 'ERROR IN MAT_SYMPLECTIFY: UNEQUAL MATRIX SIZES.'
  if (global_com%exit_on_error) call err_exit
endif

! Form the symmetrix matrix:
!         m1 = S * (1 - mat_in) * (1 + mat_in)^-1
!         m_symmetric = (m1 + m1_transpose) / 2

m3 = mat_scale_p0 (mat_in, p0_ratio, .true.)

m1 =-m3
forall (i = 1:n) m1(i,i) = m1(i,i) + 1

m2 = m3
forall (i = 1:n) m2(i,i) = m2(i,i) + 1
call mat_inverse (m2, m2)

m3 = matmul(m1, m2)

do i = 2, n, 2
  m1(i-1,1:n) =  m3(i,  1:n)
  m1(i,  1:n) = -m3(i-1,1:n)
enddo

forall (i = 1:n, j = 1:n) m_symmetric(i,j) = (m1(i,j) + m1(j,i)) / 2

! Scale by r_root

if (present(r_root)) m_symmetric = r_root * m_symmetric

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

mat_symp = mat_scale_p0 (matmul(m2, m3), p0_ratio)

end subroutine
