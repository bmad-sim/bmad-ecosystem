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
!   use bmad
!
! Input:
!   ring     -- Ring_struct: Ring to track through.
!   i_dim    -- Integer: 
!       = 2,4  Transverse closed orbit at constant energy (dE/E = CO.Z.VEL)
!       = 6    Full closed orbit using the entire transfer 6x6 matrix.
!   eps_rel(6) -- Real(rp), optional: Relative allowed error.
!                   Default is 1e-5.
!   eps_abs(6) -- Real(rp), optional: Absolute allowed error.
!                   Default is 1e-8.
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
!   closed_orb_(0:) -- Coord_struct, allocatable: closed orbit. 
!                       This routine will allocate this array for you.
!   bmad_status     -- BMAD status common block:
!       %ok           -- Set False if orbit does not converge.
!-

#include "CESR_platform.inc" 

subroutine closed_orbit_from_tracking (ring, closed_orb_, i_dim, &
                                                 eps_rel, eps_abs, init_guess)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) ring
  type (coord_struct), allocatable, intent(inout) :: closed_orb_(:)
  type (coord_struct) :: start(100), end(100)
  type (coord_struct), optional :: init_guess

  real(rp), intent(in), optional :: eps_rel(:), eps_abs(:)

  real(rp) amp_co, amp_del, orb_diff(6), amp(6)
  real(rp) mat6(6,6), mat6_inv(6,6), mat6_unit(6,6)
  real(rp) start_mat(6,6), end_mat(6,6)
  real(rp) :: error, rel_err(6), abs_err(6)

  integer i_dim, i, i1, i2, j, k, jmax, n_ele, j0, jj, nd, nnd, msk(6)

  logical :: debug = .false., rf_on

! init

  rel_err = 1e-6
  if (present(eps_rel)) rel_err = eps_rel

  abs_err = 1e-6
  if (present(eps_abs)) abs_err = eps_abs

! make sure orb_ has the correct size

  call reallocate_coord (closed_orb_, ring%n_ele_max)

! make a unit matrix

  call mat_make_unit (mat6_unit)

  n_ele = ring%n_ele_ring
  nd = i_dim

! Turn off RF voltage if i_dim == 4 (for constant delta_E)
! Make sure RF is on if i_dim = 6

  if (nd == 2 .or. nd == 4) then
    ring%ele_(:)%internal_logic = ring%ele_(:)%is_on
    call set_on (rfcavity$, ring, .false.)
  elseif (nd == 6) then
    rf_on = .false.
    do i = 1, ring%n_ele_ring
      if (ring%ele_(i)%key == rfcavity$ .and. &
                        ring%ele_(i)%value(volt$) /= 0) rf_on = .true.
    enddo
    if (.not. rf_on) then
      print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: ', &
                                      'RF IS NOT ON FOR 6-DIM TRACKING!'
      call err_exit
    endif
  else
    print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: BAD "I_DIM":', nd
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

  jmax = 100
  do j = 1, jmax

! track the current guess

    call track_all (ring, closed_orb_)
    if (ring%param%lost) exit

! save start and end coords.

    start(j)%vec = closed_orb_(0)%vec
    end(j)%vec   = closed_orb_(n_ele)%vec

! is this good enough? if so return.

    orb_diff = closed_orb_(n_ele)%vec - closed_orb_(0)%vec
    amp = max(abs(closed_orb_(0)%vec), abs(closed_orb_(n_ele)%vec))

    if (all( abs(orb_diff(1:nd)) < abs_err(1:nd) + &
                                         rel_err(1:nd) * amp(1:nd) ) ) then
      if (nd == 2 .or. nd == 4) &
                          ring%ele_(:)%is_on = ring%ele_(:)%internal_logic
      return
    endif

! if we have enough data then form the transfer matrix and invert to get 
! the new closed orbit guess. 
! do not use the matrix if it is very non-symplectic.
!
! msk is used to mask out any plans where closed_orb_ has converged to
! the closed orbit. If we do not do this then the computation can blow up.

  if (j > nd) then

    nnd = 0  ! number of unmasked dimensions
    do i = 2, nd, 2
      i1 = i-1
      i2 = i
      if (.not. all( abs(orb_diff(i1:i2)) < abs_err(i1:i2) + &
                                     rel_err(i1:i2) * amp(i1:i2) )) then
        msk(nnd+1:nnd+2) = (/ i1, i2 /)
        nnd = nnd + 2
      endif
     enddo

    do i = 1, nd
      j0 = j - nd
      jj = j0 + i
      end_mat(1:nnd,i) = end(jj)%vec(msk(1:nnd)) - end(j0)%vec(msk(1:nnd))
      start_mat(1:nnd, i) = &
                    start(jj)%vec(msk(1:nnd)) - start(j0)%vec(msk(1:nnd))
    enddo

    call mat_inverse (start_mat(1:nnd,1:nnd), start_mat(1:nnd, 1:nnd))
    mat6 = matmul(end_mat(1:nnd,1:nnd), start_mat(1:nnd,1:nnd))
    call mat_symp_check (mat6(1:nnd,1:nnd), error)
    if (debug) print *, 'error:', error
    if (debug) print '(i4, a, 3p6f11.5)', j, ':', &
                                    (closed_orb_(0)%vec(i), i = 1, nd)

    if (error < 0.1) then  ! if error is low use matrix
      mat6(1:nnd,1:nnd) = mat6_unit(1:nnd,1:nnd) - mat6(1:nnd,1:nnd)
      call mat_inverse (mat6(1:nnd,1:nnd), mat6_inv(1:nnd,1:nnd))
      closed_orb_(0)%vec(msk(1:nnd)) = closed_orb_(0)%vec(msk(1:nnd)) + &
                        matmul(mat6_inv(1:nnd,1:nnd), orb_diff(msk(1:nnd)))
      cycle
    endif

  endif


! if we are here then we did not make a guess using the matrix.
! The new guess is the average of the start and end orbits.

    closed_orb_(0)%vec(1:nd) = &
              (closed_orb_(0)%vec(1:nd) + closed_orb_(n_ele)%vec(1:nd)) / 2

  enddo

! error

  if (bmad_status%type_out) then
    if (j == jmax+1) then
      print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: ORBIT NOT CONVERGING!'
    else
      print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: PROBLEM TRACKING ORBIT.'
    endif
  endif

  if (bmad_status%exit_on_error) call err_exit
  bmad_status%ok = .false.

end subroutine
