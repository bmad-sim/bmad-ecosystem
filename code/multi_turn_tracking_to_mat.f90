!+
! Subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, track0, chi)
!
! Subroutine to analyze multi-turn tracking data to find the 1-turn transfer
! matrix and the closed orbit offset at a given point in the ring.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!             
! Input:
!   track(:) -- Coord_struct: multi-turn tracking data to analyze.
!                track(i) is the particle position at a given point
!                in the ring on the i^th turn.
!   i_dim    -- Integer: Dimensionality of the data (2, 4, or 6).
!
! Output: 
!   mat1(:,:) -- Real: Calculated 1-turn matrix.
!   track0    -- Coord_struct: Closed orbit offset.
!   chi       -- Real: Figure of merit in the fitting:
!                   = 0 => perfect fit
!                   = 1 => terrible fit
!-

subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, track0, chi)

  use local_bmad_struct
  use bmad_interface
  use nr

  implicit none

  type (coord_struct), intent(in), target :: track(:)
  type (coord_struct), intent(out) :: track0
  real, intent(out) :: mat1(:,:)
  real, intent(out) :: chi
  integer, intent(in) :: i_dim

  real sum2, dsum2, chisq, dtrack(6), remainder(6)
  real, allocatable, save :: x(:), y(:), sig(:), v(:,:), w(:), a(:), m(:,:)
  type (coord_struct), allocatable, target, save :: d0track(:)
  integer i, n

  external multi_turn_func

! init

  n = size(track)

  if (.not. allocated (x)) then
    allocate (x(n-1), y(n-1), sig(n-1), d0track(n-1))
  elseif (size(x) /= n-1) then
    deallocate (x, y, sig)
    allocate (x(n-1), y(n-1), sig(n-1), d0track(n-1))
  endif

  if (.not. allocated (w)) then
    allocate (w(i_dim+1), v(i_dim+1, i_dim+1), a(i_dim+1), m(i_dim,i_dim))
  elseif (size(w) /= i_dim+1) then
    deallocate (w, v, a, m)
    allocate (w(i_dim+1), v(i_dim+1, i_dim+1), a(i_dim+1), m(i_dim,i_dim))
  endif

  x = (/ (i, i=1,n-1) /)
  sig = 1

! because of possible round-off errors we do the computation in two parts.
! first compute the closed orbit.

  do i = 1, i_dim
    y = track(2:n)%vec(i)
    multi_turn_func_common => track(1:n-1)
    call svdfit (x, y, sig, a, v, w, chisq, multi_turn_func)
    mat1(i,1:i_dim) = a(1:i_dim)
    remainder(i) = a(i_dim+1)
  enddo

  m = -mat1(1:i_dim,1:i_dim)
  forall (i = 1:i_dim) m(i,i) = 1 + m(i,i)
  call mat_inverse (m, m)
  track0%vec(1:i_dim) = matmul (m, remainder)

! second subtrack off the closed orbit and calculate the 1-turn matrix

  do i = 1, n-1
    d0track(i)%vec = track(i)%vec - track0%vec
  enddo

  do i = 1, i_dim
    y = track(2:n)%vec(i) - track0%vec(i)
    multi_turn_func_common => d0track
    call svdfit (x, y, sig, a, v, w, chisq, multi_turn_func)
    mat1(i,1:i_dim) = a(1:i_dim)
  enddo

! calculate chi -  the goodness of fit

  sum2 = 0
  dsum2 = 0

  do i = 1, n-1
    dtrack(1:i_dim) = track(i+1)%vec(1:i_dim) - &
                 matmul(mat1(1:i_dim,1:i_dim), track(i)%vec(1:i_dim)) - &
                 remainder(1:i_dim)
    dsum2 = dsum2 + sum(dtrack(1:i_dim)**2)
    sum2 = sum(track(i+1)%vec(1:i_dim)**2)
  enddo

  chi = sqrt(dsum2/sum2)

end subroutine


!-------------------------------------------------------------------------

function multi_turn_func (x, n)

  use local_bmad_struct

  implicit none                      

  real, intent(in) :: x
  integer, intent(in) :: n
  real, dimension(n) :: multi_turn_func

  multi_turn_func = (/ multi_turn_func_common(nint(x))%vec(1:n-1), 1.0 /)

end function
