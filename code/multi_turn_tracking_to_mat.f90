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
!   mat1(:,:) -- Real(rdef): Calculated 1-turn matrix.
!   track0    -- Coord_struct: Closed orbit offset.
!   chi       -- Real(rdef): Figure of merit in the fitting:
!                   = 0 => perfect fit
!                   = 1 => terrible fit
!-

!$Id$
!$Log$
!Revision 1.5  2003/01/27 14:40:39  dcs
!bmad_version = 56
!
!Revision 1.4  2002/06/13 14:54:26  dcs
!Interfaced with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:20  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:54  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine multi_turn_tracking_to_mat (track, i_dim, mat1, track0, chi)

  use bmad_struct
  use bmad_interface
  use nr

  implicit none

  type (coord_struct), intent(in), target :: track(:)
  type (coord_struct), intent(out) :: track0
  real(rdef), intent(out) :: mat1(:,:)
  real(rdef), intent(out) :: chi
  integer, intent(in) :: i_dim

  real(rdef) sum2, dsum2, chisq, dtrack(6), remainder(6)
  real(rdef), allocatable, save :: x(:), y(:), sig(:), v(:,:), w(:), a(:), m(:,:)
  type (coord_struct), allocatable, target, save :: d0track(:)
  integer i, n

  external multi_turn_func

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

  use bmad

  implicit none                      

  real(rdef), intent(in) :: x
  integer, intent(in) :: n
  real(rdef), dimension(n) :: multi_turn_func

  multi_turn_func = (/ multi_turn_func_common(nint(x))%vec(1:n-1), 1.0_rdef /)

end function
