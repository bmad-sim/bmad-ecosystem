!+
! Subroutine mat_symp_check (mat, error)
!
! Routine to check the symplecticity of a square matrix. The error is
! defined to via:
!     error = maxval (abs (Mat_transpose * S * Mat - S))
!
! Modules Needed:
!   use bmad_interface
!
! Input:
!   mat(:,:) -- Real: Matrix to check
!
! Output:
!   error -- Real: difference from symplecticity:
!             = 0    --> perfect.
!             = 1e-4 --> Reasonably good.
!             = 1    --> Terrible.
!-


subroutine mat_symp_check (mat, error)

  use bmad_struct
  use bmad_interface

  implicit none

  integer i, j, n

  real, intent(in) :: mat(:,:)
  real, allocatable, save :: m2(:,:)
  real error
      
  logical :: debug = .false.

!

  n = ubound(mat, 1)

  if (mod(n, 2) /= 0) then
    type *, 'ERROR IN MAT_SYMP_CHECK: MATRIX DOES NOT HAVE EVEN SIZE'
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
    type *, 'MAT_SYMP_CHECK:', error
    do i = 1, n
      type '(5x, (<n>f11.5))', (max(min(m2(i, j), 999.0), -999.0), j = 1, n)
    enddo
  endif

end subroutine
