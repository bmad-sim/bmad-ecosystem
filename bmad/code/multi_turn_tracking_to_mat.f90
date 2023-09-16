!+
! Subroutine multi_turn_tracking_to_mat (track, n_var, map1, map0, track0, chi)
!
! Subroutine to analyze multi-turn tracking data to find the 1-turn transfer
! matrix and the closed orbit offset at a given point in the lat.
!             
! Input:
!   track(:) -- Coord_struct: multi-turn tracking data to analyze. track(i) is the 
!                particle position at a given point in the lat on the i^th turn.
!   n_var    -- Integer: Dimensionality of the data (2, 4, or 6).
!
! Output: 
!   map1(:,:) -- Real(rp): Calculated 1-turn matrix. 
!   map0(:)   -- Real(rp): 0th order part of the tranport map.
!   track0    -- Coord_struct: Closed orbit offset.
!   chi       -- Real(rp): Figure of merit in the fitting:
!                   = 0 => perfect fit
!                   = 1 => terrible fit
!-

subroutine multi_turn_tracking_to_mat (track, n_var, map1, map0, track0, chi)

use bmad_interface, except_dummy => multi_turn_tracking_to_mat

implicit none

type (coord_struct), intent(in), target :: track(:)
type (coord_struct), intent(out) :: track0

real(rp), intent(out) :: map1(:,:), map0(:)
real(rp), intent(out) :: chi
integer, intent(in) :: n_var

real(rp) sum2, dsum2, chisq, dtrack(6)
real(rp), allocatable :: a(:,:), a2(:,:), m(:,:), y(:), x(:)
integer i, n_dat

! init

n_dat = size(track)
allocate (y(n_dat-1), A(n_dat-1, n_var+1), A2(n_dat-1, n_var+1), m(n_var,n_var), x(n_var+1))

! Because of possible round-off errors the computation is done in two parts.
!
! First: Compute the closed orbit.
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

a = 0
do i = 1, n_dat-1
  a(i,:) = [track(i)%vec(1:n_var), 1.0_rp]
enddo

do i = 1, n_var
  a2 = a
  y = track(2:n_dat)%vec(i)
  call svd_fit (a2, y, 1.0e-5_rp, x)
  map1(i,1:n_var) = x(1:n_var)
  map0(i) = x(n_var+1)
enddo

m = -map1(1:n_var,1:n_var)
forall (i = 1:n_var) m(i,i) = 1 + m(i,i)
call mat_inverse (m, m)
track0%vec(1:n_var) = matmul (m, map0(1:n_var))

! Second: Subtract off the closed orbit and calculate the 1-turn matrix.

do i = 1, n_dat-1
  a(i,1:n_var) = track(i)%vec(1:n_var)-track0%vec(1:n_var)
enddo

do i = 1, n_var
  a2 = a
  y(1:n_dat-1) = track(2:n_dat)%vec(i) - track0%vec(i)
  call svd_fit (a2(:,1:n_var), y, 1.0e-5_rp, x(1:n_var))
  map1(i,1:n_var) = x(1:n_var)
enddo

! Calculate chi -  the goodness of fit

sum2 = 0
dsum2 = 0

do i = 1, n_dat-1
  dtrack(1:n_var) = track(i+1)%vec(1:n_var) - &
                        matmul(map1(1:n_var,1:n_var), track(i)%vec(1:n_var)) - map0(1:n_var)
  dsum2 = dsum2 + sum(dtrack(1:n_var)**2)
  sum2 = sum(track(i+1)%vec(1:n_var)**2)
enddo

chi = sqrt(dsum2/sum2)

end subroutine

