!+
! Subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, track0, chi)
!
! Subroutine to analyze multi-turn tracking data to find the 1-turn transfer
! matrix and the closed orbit offset at a given point in the ring.
!
! Modules needed:
!   use bmad
!             
! Input:
!   track(:) -- Coord_struct: multi-turn tracking data to analyze.
!                track(i) is the particle position at a given point
!                in the ring on the i^th turn.
!   i_dim    -- Integer: Dimensionality of the data (2, 4, or 6).
!
! Output: 
!   mat1(:,:) -- Real(rp): Calculated 1-turn matrix.
!   track0    -- Coord_struct: Closed orbit offset.
!   chi       -- Real(rp): Figure of merit in the fitting:
!                   = 0 => perfect fit
!                   = 1 => terrible fit
!-

#include "CESR_platform.inc"

subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, track0, chi)

  use bmad_struct
  use bmad_interface, except => multi_turn_tracking_to_mat
  use nr

  implicit none

  type (coord_struct), intent(in), target :: track(:)
  type (coord_struct), intent(out) :: track0
  real(rp), intent(out) :: mat1(:,:)
  real(rp), intent(out) :: chi
  integer, intent(in) :: i_dim

  real(rp) sum2, dsum2, chisq, dtrack(6), remainder(6)
  real(rp), allocatable, save :: x(:), y(:), sig(:), v(:,:), &
                                   w(:), a(:), m(:,:)
  type (coord_struct), allocatable, target, save :: d0track(:)
  integer i, n

  interface
    function multi_turn_func (x, n)
      use precision_def
      real(rp), intent(in) :: x
      integer, intent(in) :: n
      real(rp) :: multi_turn_func(n)
    end function
  end interface

! init

  n = size(track)

  if (.not. allocated (x)) then
    allocate (x(n-1), y(n-1), sig(n-1), d0track(n-1))
  elseif (size(x) /= n-1) then
    deallocate (x, y, sig, d0track)
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

! Because of possible round-off errors we do the computation in two parts.
! First compute the closed orbit.
! Use the relation:
!            V_out = M V_in + (1 - M) C
! Now extend the dimensions of everything by 1:
!           VV_in  = (V_in,  1)
!           VV_out = (V_out, 1)
!           MM = | M  (1-M)C |
!                | 0    1    |
! Then:
!           VV_out = MM VV_in
! Which can be solved using svdfit.


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

! Second subtrack off the closed orbit and calculate the 1-turn matrix

  do i = 1, n-1
    d0track(i)%vec = track(i)%vec - track0%vec
  enddo

  do i = 1, i_dim
    y = track(2:n)%vec(i) - track0%vec(i)
    multi_turn_func_common => d0track
    call svdfit (x, y, sig, a, v, w, chisq, multi_turn_func)
    mat1(i,1:i_dim) = a(1:i_dim)
  enddo

! Calculate chi -  the goodness of fit

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
!-------------------------------------------------------------------------

function multi_turn_func (x, id)

  use bmad_struct
  use bmad_interface

  implicit none                      

  real(rp), intent(in) :: x
  integer, intent(in) :: id
  real(rp), dimension(id) :: multi_turn_func

! id = i_dim+1

  multi_turn_func = (/ multi_turn_func_common(nint(x))%vec(1:id-1), 1.0_rp /)

end function
