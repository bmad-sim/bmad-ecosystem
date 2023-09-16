!+
! Subroutine closed_orbit_from_tracking (lat, closed_orb, i_dim, eps_rel, eps_abs, init_guess, err_flag)
!
! Subroutine to find the closed orbit via tracking. 
!
! See closed_orbit_calc for an alternative routine.
! This routine has the advantage of not using the 1-turn transfer matrix 
! (which is calculated from the transfer matrices of the elements).
! This routine is probably slower than closed_orbit_calc (especially for 
! linear lattices) but (hopefully) will be better behaved when there
! are strong nonlinearities. Possible convergence porblems with this 
! routine occurs when the initial guess is far from the true closed orbit.
!
! Input:
!   lat     -- lat_struct: Lat to track through.
!   i_dim    -- Integer: 
!       = 2,4  Transverse closed orbit at constant energy.
!       = 6    Full closed orbit using the entire transfer 6x6 matrix.
!   eps_rel(6) -- Real(rp), optional: Relative allowed error.
!                   Default is bmad_com%rel_tol_tracking
!   eps_abs(6) -- Real(rp), optional: Absolute allowed error.
!                   Default is bmad_com%abs_tol_tracking
!   init_guess -- Coord_struct, optional: Starting guess for the closed orbit at the start of the lattice.
!                   Set init_guess%vec(6) to the appropriate value of pz when calculating off-energy orbits.
!                   If not present then the origin will be used. 
!
! Output:
!   closed_orb(0:) -- Coord_struct, allocatable: closed orbit. 
!                       This routine will allocate this array for you.
!   err_flag       -- Logical, optional: Set True if there is an error. False otherwise.
!-

subroutine closed_orbit_from_tracking (lat, closed_orb, i_dim, eps_rel, eps_abs, init_guess, err_flag)

use bmad_interface, except_dummy => closed_orbit_from_tracking

implicit none

type (lat_struct) lat
type (coord_struct), allocatable :: closed_orb(:)
type (coord_struct) :: start(100), end(100)
type (coord_struct), optional :: init_guess

real(rp), intent(in), optional :: eps_rel(:), eps_abs(:)

real(rp) orb_diff(6), amp(6)
real(rp) mat6(6,6), mat6_inv(6,6), mat6_unit(6,6)
real(rp) start_mat(6,6), end_mat(6,6)
real(rp) :: error, rel_err(6), abs_err(6)
real(rp), allocatable :: on_off_save(:)

integer i_dim, i, i1, i2, j, jmax, n_ele, j0, jj, nd, nnd, msk(6), track_state

logical, optional :: err_flag
logical :: debug = .false., rf_on, fluct_saved, aperture_saved

character(*), parameter :: r_name = 'closed_orbit_from_tracking'

! init

if (present(err_flag)) err_flag = .true.

rel_err = bmad_com%rel_tol_tracking
if (present(eps_rel)) rel_err = eps_rel

abs_err = bmad_com%abs_tol_tracking
if (present(eps_abs)) abs_err = eps_abs

! make sure orb has the correct size

call reallocate_coord (closed_orb, lat%n_ele_max)

fluct_saved = bmad_com%radiation_fluctuations_on
bmad_com%radiation_fluctuations_on = .false.  

aperture_saved = bmad_com%aperture_limit_on
bmad_com%aperture_limit_on = .false.

! make a unit matrix

call mat_make_unit (mat6_unit)

n_ele = lat%n_ele_track
nd = i_dim

! Turn off RF voltage if i_dim == 4 (for constant delta_E)
! Make sure RF is on if i_dim = 6

if (nd == 2 .or. nd == 4) then
  call set_on_off (rfcavity$, lat, off_and_save$, saved_values = on_off_save)
elseif (nd == 6) then
  rf_on = .false.
  do i = 1, lat%n_ele_track
    if (lat%ele(i)%key == rfcavity$ .and. &
                      lat%ele(i)%value(voltage$) /= 0) rf_on = .true.
  enddo
  if (.not. rf_on) then
    print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: ', 'RF IS NOT ON FOR 6-DIM TRACKING!'
    return
  endif
else
  print *, 'ERROR IN CLOSED_ORBIT_FROM_TRACKING: BAD "I_DIM":', nd
  return
endif

! determine initial guess

if (present(init_guess)) then
  closed_orb(0) = init_guess
else
  closed_orb(0)%vec = 0
endif

!-------------------
! iterate until we find the closed orbit

jmax = 100
do j = 1, jmax

  ! track the current guess

  call track_all (lat, closed_orb, track_state = track_state)
  if (track_state /= moving_forward$) exit

  ! save start and end coords.

  start(j)%vec = closed_orb(0)%vec
  end(j)%vec   = closed_orb(n_ele)%vec

  ! is this good enough? if so return.

  orb_diff = closed_orb(n_ele)%vec - closed_orb(0)%vec
  amp = max(abs(closed_orb(0)%vec), abs(closed_orb(n_ele)%vec))

  if (all( abs(orb_diff(1:nd)) < abs_err(1:nd) + rel_err(1:nd) * amp(1:nd) ) ) then
    if (nd == 2 .or. nd == 4) then
      call set_on_off (rfcavity$, lat, restore_state$, saved_values = on_off_save)
    endif
    bmad_com%radiation_fluctuations_on = fluct_saved  ! restore state
    bmad_com%aperture_limit_on = aperture_saved
    if (present(err_flag)) err_flag = .false.
    return
  endif

  ! If we have enough data then form the transfer matrix and invert to get 
  ! the new closed orbit guess. 
  ! do not use the matrix if it is very non-symplectic.
  !
  ! msk is used to mask out any plans where closed_orb has converged to
  ! the closed orbit. If we do not do this then the computation can blow up.

  if (j > nd) then

    nnd = 0  ! number of unmasked dimensions
    do i = 2, nd, 2
      i1 = i-1
      i2 = i
      if (.not. all( abs(orb_diff(i1:i2)) < abs_err(i1:i2) + rel_err(i1:i2) * amp(i1:i2) )) then
        msk(nnd+1:nnd+2) = [i1, i2 ]
        nnd = nnd + 2
      endif
     enddo

    do i = 1, nd
      j0 = j - nd
      jj = j0 + i
      end_mat(1:nnd,i) = end(jj)%vec(msk(1:nnd)) - end(j0)%vec(msk(1:nnd))
      start_mat(1:nnd, i) = start(jj)%vec(msk(1:nnd)) - start(j0)%vec(msk(1:nnd))
    enddo

    call mat_inverse (start_mat(1:nnd,1:nnd), start_mat(1:nnd, 1:nnd))
    mat6(1:nnd,1:nnd) = matmul(end_mat(1:nnd,1:nnd), start_mat(1:nnd,1:nnd))
    error = mat_symp_error (mat6(1:nnd,1:nnd))
    if (debug) print *, 'error:', error
    if (debug) print '(i4, a, 3p6f11.5)', j, ':', (closed_orb(0)%vec(i), i = 1, nd)

    if (error < 0.1) then  ! if error is low use matrix
      mat6(1:nnd,1:nnd) = mat6_unit(1:nnd,1:nnd) - mat6(1:nnd,1:nnd)
      call mat_inverse (mat6(1:nnd,1:nnd), mat6_inv(1:nnd,1:nnd))
      closed_orb(0)%vec(msk(1:nnd)) = closed_orb(0)%vec(msk(1:nnd)) + &
                        matmul(mat6_inv(1:nnd,1:nnd), orb_diff(msk(1:nnd)))
      cycle
    endif

  endif

  ! if we are here then we did not make a guess using the matrix.
  ! The new guess is the average of the start and end orbits.

  closed_orb(0)%vec(1:nd) = (closed_orb(0)%vec(1:nd) + closed_orb(n_ele)%vec(1:nd)) / 2

enddo

! error

if (j == jmax+1) then
  call out_io (s_error$, r_name, 'ORBIT NOT CONVERGING!')
else
  call out_io (s_error$, r_name, 'PROBLEM TRACKING ORBIT.')
endif

end subroutine
