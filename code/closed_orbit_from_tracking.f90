!+
! Subroutine closed_orbit_from_tracking (ring, closed_orb_, i_dim, 
!                                                eps_rel, eps_abs, init_guess)
!
! Subroutine to find the closed orbit via tracking. 
!
! See closed_orbit_at_start for an alternative routine.
! This routine has the advantage of not using the 1-turn transfer matrix 
! (which is calculated from the transfer matrices of the elements).
! This routine is probably slower than closed_orbit_at_start (especially for 
! linear lattices) but (hopefully) will be better behaved when there
! are strong nonlinearities. Possible convergence porblems with this 
! routine occurs when the initial guess is far from the true closed orbit.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring     -- Ring_struct: Ring to track through.
!   i_dim    -- Integer: 
!       = 4  Transverse closed orbit at constant energy (dE/E = CO.Z.VEL)
!       = 6  Full closed orbit using the entire transfer 6x6 matrix.
!   eps_rel(:) -- Real: relative allowed error.
!   eps_abs(:) -- Real: absolute allowed error.
!   init_guess -- [Optional] Coord_struct: Starting guess for the closed 
!                orbit at the start of the ring. If not present then
!                the origin will be used. 
!   bmad_status -- BMAD status common block:
!         %exit_on_error -- If True then subroutine will terminate program
!                           if the orbit does not converge.
!         %type_out      -- If True then the subroutine will type out
!                           a warning message if the orbit does not converge.
!
!
! Output:
!   closed_orb_(0:n_ele_maxx) -- Coord_struct: closed orbit.
!   bmad_status -- BMAD status common block:
!       %ok         -- Set False if orbit does not converge.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:49  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine closed_orbit_from_tracking (ring, closed_orb_, i_dim, &
                                                 eps_rel, eps_abs, init_guess)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (coord_struct), intent(inout) :: closed_orb_(0:n_ele_maxx)
  type (coord_struct) :: orb_(0:n_ele_maxx)
  type (coord_struct), optional :: init_guess

  real, intent(in) :: eps_rel(:), eps_abs(:)

  real amp_co, amp_del, d_orb(6), orb_diff(6), amp(6)
  real mat6(6,6), mat6_inv(6,6), mat6_unit(6,6)

  real :: error, factor = 1

  integer i_dim, i, j, n_ele

  logical is_on(0:n_ele_maxx)

! make a unit matrix

  call mat_unit (mat6_unit, 6, 6)

  n_ele = ring%n_ele_ring

! Turn off RF voltage if i_dim == 4 (for constant delta_E)

  if (i_dim == 4) then
    is_on = .false.
    do i = 1, n_ele
      if (ring%ele_(i)%key == rfcavity$) then
        is_on(i) = ring%ele_(i)%is_on
        ring%ele_(i)%is_on = .false.
      endif
    enddo
  elseif (i_dim /= 6) then
    type *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: BAD "I_DIM":', i_dim
    call err_exit    
  endif

! determine initial guess

  if (present(init_guess)) then
    closed_orb_(0) = init_guess
  else
    closed_orb_(0)%vec = 0
  endif

!-------------------
! iterate until we find the closed orbit

  do j = 1, 100

! track the current guess

    call track_all (ring, closed_orb_)
    if (.not. bmad_status%ok) exit

! is this good enough? if so return.

    orb_diff = closed_orb_(n_ele)%vec - closed_orb_(0)%vec

!    amp_co = sum(abs(closed_orb_(0)%vec(1:i_dim)))
!    amp_del = sum(abs(orb_diff(1:i_dim)))
!    if (amp_del < amp_co*1.0e-5 + 1.0e-8) then  ! success

    amp = max(abs(closed_orb_(0)%vec), abs(closed_orb_(n_ele)%vec))
    if (all(orb_diff(1:i_dim) < eps_abs(1:i_dim)) .or. &
          all(orb_diff(1:i_dim) < eps_rel(1:i_dim) * amp(1:i_dim))) then
      if (i_dim == 4) where (is_on(1:n_ele)) ring%ele_(1:n_ele)%is_on = .true.
      return
    endif

! calculate a new guess based upon tracking nearby orbits.  
! the delta we use to offset the orbit is based upon the current guess.

    d_orb = 1e-4 * amp
    d_orb = max(d_orb, (/ 1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 1e-4 /))

    do i = 1, i_dim
      orb_(0) = closed_orb_(0)
      orb_(0)%vec(i) = orb_(0)%vec(i) + d_orb(i)
      call track_all (ring, orb_)
      if (.not. bmad_status%ok) exit
      mat6(1:6,i) = (orb_(n_ele)%vec - closed_orb_(n_ele)%vec) / d_orb(i)
    enddo

    call mat_symp_check (mat6(1:i_dim,1:i_dim), error)
    type *, 'error:', error
    type '(i4, a, 3p6f11.5)', j, ':', (closed_orb_(0)%vec(i), i = 1, i_dim)

    mat6 = mat6_unit - mat6
    call mat_inv (mat6, mat6_inv, i_dim, 6)
    closed_orb_(0)%vec(1:i_dim) = closed_orb_(0)%vec(1:i_dim) + &
                factor * matmul(mat6_inv(1:i_dim,1:i_dim), orb_diff(1:i_dim))

  enddo

! error

  if (bmad_status%type_out) then
    if (j == 101) then
      type *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: ORBIT NOT CONVERGING!'
    else
      type *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: PROBLEM TRACKING ORBIT.'
    endif
  endif

  if (bmad_status%exit_on_error) call err_exit
  bmad_status%ok = .false.

end subroutine
