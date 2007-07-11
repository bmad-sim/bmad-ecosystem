#include "CESR_platform.inc"

module matrix_mod

  use precision_def
  use dcslib
  use bmad_basic_mod

  integer, parameter :: ok$              = 1
  integer, parameter :: in_stop_band$    = 2
  integer, parameter :: non_symplectic$  = 3
  integer, parameter :: unstable$        = 4

  character(16) :: status_name(5) = (/     'OK            ', &
      'IN_STOP_BAND  ', 'NON_SYMPLECTIC', 'UNSTABLE      ', '              ' /)

contains
      
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine mat_inverse (mat, mat_inv, err_flag, print_error)
!
! Takes the inverse of a square matrix using LU Decomposition from
! Numerical Recipes.
!
! Modules needed:
!   bmad_interface
!
! Input:
!   mat(:,:)     -- Real(rp): Input matrix array
!   print_err    -- Logical, optional: If True then the subroutine will type out
!                         a warning message. Default is False.
! Output:
!   mat_inv(:,:) -- Real(rp): inverse of mat1
!   err_flag     -- Logical, optional: Set true for a singular input matrix.
!-

subroutine mat_inverse (mat, mat_inv, err_flag, print_err)

  use nr

  implicit none

  real(rp) :: mat(:,:)
  real(rp) :: mat_inv(:,:)

  real(rp), allocatable, save :: mat2(:,:), vec(:)
  integer, allocatable, save :: indx(:)
  real(rp) d

  integer n, i
  logical, optional :: err_flag, print_err
  character(16) :: r_name = 'mat_inverse'

!

  n = size (mat, 1)

  if (.not. allocated(indx)) then
    allocate (mat2(n,n), indx(n), vec(n))
  elseif (size(indx) /= n) then
    deallocate (mat2, indx, vec)
    allocate (mat2(n,n), indx(n), vec(n))
  endif

  mat2 = mat  ! use temp mat so as to not change mat

  vec(1:n) = maxval(abs(mat), dim = 2)
  if (any(vec(1:n) == 0)) then
    if (logic_option(.false., print_err)) &
                                call out_io (s_error$, r_name, 'SINGULAR MATRIX.')
    if (present(err_flag)) err_flag = .false.
    return
  endif
  if (present(err_flag)) err_flag = .true.

  call ludcmp (mat2, indx, d)

  mat_inv(1:n,1:n) = 0
  forall (i = 1:n) mat_inv(i,i) = 1

  do i = 1, n
    call lubksb (mat2, indx, mat_inv(1:n,i))
  enddo

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine mat_symp_check (mat, error)
!
! Routine to check the symplecticity of a square matrix. The error is
! defined to via:
!     error = maxval (abs (Mat_transpose * S * Mat - S))
!
! Modules Needed:
!
! Input:
!   mat(:,:) -- Real(rp): Matrix to check
!
! Output:
!   error -- Real(rp): difference from symplecticity:
!             = 0    --> perfect.
!             = 1e-4 --> Fair, but I wouldn't use for long term tracking.
!             = 1    --> Terrible.
!-

subroutine mat_symp_check (mat, error)

  implicit none

  integer i, j, n

  real(rp), intent(in) :: mat(:,:)
  real(rp), allocatable, save :: m2(:,:)
  real(rp) error

  logical :: debug = .false.

!

  n = ubound(mat, 1)

  if (mod(n, 2) /= 0) then
    print *, 'ERROR IN MAT_SYMP_CHECK: MATRIX DOES NOT HAVE EVEN SIZE'
    call err_exit
  endif

! init

  if (.not. allocated (m2)) then
    allocate (m2(n,n))
  elseif (ubound(m2, 1) /= n) then
    deallocate (m2)
    allocate (m2(n,n))
  endif

!

  do j = 2, n, 2                  ! m2 = mat_transpose * S
    m2(1:n,j)   = -mat(j-1,1:n)
    m2(1:n,j-1) =  mat(j,  1:n)
  enddo

  m2 = matmul (m2, mat)

  do j = 2, n, 2
    m2(j, j-1) = m2(j, j-1) - 1
    m2(j-1, j) = m2(j-1, j) + 1
  enddo

  error = maxval(abs(m2))

  if (debug) then
    print *, 'MAT_SYMP_CHECK:', error
    do i = 1, n
      print '(5x, (6f11.5))', &
                        (max(min(m2(i, j), 999.0_rp), -999.0_rp), j = 1, n)
    enddo
  endif

end subroutine

!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Subroutine mat_symplectify (mat_in, mat_symp)
!
! Subroutine to form a symplectic matrix that is  "close" to the input matrix.
!
! Modules needed:
!
! Input:
!   mat_in(:,:) -- Real(rp): Input matrix to symplectify.
!
! Output:
!   mat_symp(:,:) -- Real(rp): Symplectic output matrix
!-

subroutine mat_symplectify (mat_in, mat_symp)

  implicit none

  real(rp), intent(in)  :: mat_in(:,:)
  real(rp), intent(out) :: mat_symp(:,:)

  real(rp), save, allocatable :: m1(:,:), m2(:,:), m3(:,:), m_symetric(:,:)

  integer n,i,j

!

  n = ubound(mat_in, 1)

  if (ubound(mat_symp, 1) /= n) then
    print *, 'ERROR IN MAT_SYMPLECTIFY: UNEQUAL MATRIX SIZES.'
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

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine mat6_dispersion (e_vec, mat6)
!
! Subroutine to put the dispersion into ele%mat6 given the eta vector e_vec
!
! Input:
!   e_vec(4) -- Real(rp): Dispersion vector
!   mat6(6,6) -- Real(rp): Matrix with 4x4 x-y submatrix already made.
!
! Output:
!   mat6(6,6) -- Real(rp): mat6(5, 1:4) components set for the dispersion.
!-

subroutine mat6_dispersion (e_vec, mat6)

  implicit none

  real(rp), intent(inout) :: mat6(:,:)
  real(rp), intent(in) :: e_vec(:)

  real(rp) e2_vec(4)

!

  mat6(1:4, 6) = e_vec(1:4)

  e2_vec(1) = -e_vec(2)
  e2_vec(2) =  e_vec(1)
  e2_vec(3) = -e_vec(4)
  e2_vec(4) =  e_vec(3)

  mat6(5,1:4) = matmul (e2_vec, mat6(1:4,1:4))

end subroutine

end module
