#include "CESR_platform.inc"      

module runge_kutta_mod

  use em_field_mod

contains

!-----------------------------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------
!+
! Subroutine odeint_bmad (start, ele, param, end, s1, s2, &
!                            rel_tol, abs_tol, h1, hmin, local_ref_frame, track)
! 
! Subroutine to do Runge Kutta tracking. This routine is adapted from Numerical
! Recipes.  See the NR book for more details.
!
! Notice that this routine has an two tolerance arguments rel_tol and abs_tol.
! Odein only has 1. rel_tol (essentually equivalent to eps in odeint) 
! is scalled by the step size to to able to relate it to the final accuracy.
!
! Essentually (assuming random errors) one of these conditions holds:
!      %error in tracking < rel_tol
! or
!     absolute error in tracking < abs_tol
!
! Modules needed:
!   use bmad
!
! Input: 
!   start   -- Coord_struct: Starting coords: (x, x', y, y', z, delta).
!   ele     -- Ele_struct: Element to track through.
!     %tracking_method -- Determines which subroutine to use to calculate the 
!                         field. Note: BMAD does no supply em_field_custom.
!                           == custom$ then use em_field_custom
!                           /= custom$ then use em_field_standard
!   param   -- lat_param_struct: Beam parameters.
!     %enegy       -- Energy in GeV
!     %particle    -- Particle type [positron$, or electron$]
!   s1      -- Real: Starting point.
!   s2      -- Real: Ending point.
!   rel_tol -- Real: Same as eps for odeint scalled by sqrt(h/(s2-s1))
!               where h is the step size for the current step. rel_tol
!               sets the %error of the result
!   abs_tol -- Real: Sets the absolute error of the result
!   h1      -- Real: Initial guess for a step size.
!   h_min   -- Real: Minimum step size (can be zero).
!   local_ref_frame 
!           -- Logical: If True then take the 
!                input and output coordinates as being with 
!                respect to the frame of referene of the element. 
!
!   track   -- Track_struct: Structure holding the track information.
!     %save_track -- Logical: Set True if track is to be saved.
!
! Output:
!   end     -- Coord_struct: Ending coords: (x, x', y, y', z, delta).
!   track   -- Track_struct: Structure holding the track information.
!
! Common block:
!-

subroutine odeint_bmad (start, ele, param, end, s1, s2, &
                    rel_tol, abs_tol, h1, h_min, local_ref_frame, track)

  implicit none

  type (coord_struct), intent(in) :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct) ele
  type (lat_param_struct) param
  type (track_struct) track

  real(rp), intent(in) :: s1, s2, rel_tol, abs_tol, h1, h_min
  real(rp), parameter :: tiny = 1.0e-30_rp
  real(rp) :: h, h_did, h_next, s, s_sav, rel_tol_eff, abs_tol_eff, sqrt_N
  real(rp) :: dr_ds(6), r(6), r_scal(6)

  integer, parameter :: max_step = 10000
  integer :: n_step

  logical local_ref_frame

! init

  s = s1
  h = sign(h1, s2-s1)
  r = start%vec

  if (track%save_track) then
    s_sav = s - 2.0_rp * track%ds_save
    call allocate_saved_orbit (track, int(abs(s2-s1)/track%ds_save)+1)
  endif

! now track

  bmad_status%ok = .true.

  do n_step = 1, max_step

    call em_field_kick (ele, param, s, r, local_ref_frame, dr_ds)

    sqrt_N = sqrt(abs((s2-s1)/h))  ! number of steps we would take with this h
    rel_tol_eff = rel_tol / sqrt_N
    abs_tol_eff = abs_tol / sqrt_N
    r_scal(:) = abs(r(:)) + abs(h*dr_ds(:)) + TINY

    if (track%save_track .and. (abs(s-s_sav) > track%ds_save)) &
                             call save_a_step (track, ele, param, s, r, s_sav)
    if ((s+h-s2)*(s+h-s1) > 0.0) h = s2-s

    call rkqs_bmad (ele, param, r, dr_ds, s, h, rel_tol_eff, abs_tol_eff, &
                                               r_scal, h_did, h_next, local_ref_frame)
    if (.not. bmad_status%ok) return

    if (h_did == h) then
      track%n_ok = track%n_ok + 1
    else
      track%n_bad = track%n_bad + 1
    end if

    if ((s-s2)*(s2-s1) >= 0.0) then
      end%vec = r
      if (track%save_track) call save_a_step (track, ele, param, s, r, s_sav)
      return
    end if

    if (abs(h_next) < h_min) then
      bmad_status%ok = .false.
      if (bmad_status%type_out) print *, &
                'ERROR IN ODEINT_BMAD: STEPSIZE SMALLER THAN MINIMUM.' 
      if (bmad_status%exit_on_error) call err_exit
    endif

    h = h_next

  end do

  bmad_status%ok = .false.
  if (bmad_status%type_out) print *, 'ERROR IN ODEINT_BMAD: TOO MANY STEPS'
  if (bmad_status%exit_on_error) call err_exit

end subroutine odeint_bmad

!-----------------------------------------------------------------
!-----------------------------------------------------------------

subroutine rkqs_bmad (ele, param, r, dr_ds, s, h_try, rel_tol, abs_tol, &
                                       r_scal, h_did, h_next, local_ref_frame)

  implicit none

  type (ele_struct) ele
  type (lat_param_struct) param

  real(rp), intent(inout) :: r(6)
  real(rp), intent(in)    :: dr_ds(6), r_scal(6)
  real(rp), intent(inout) :: s
  real(rp), intent(in)    :: h_try, rel_tol, abs_tol
  real(rp), intent(out)   :: h_did, h_next

  real(rp) :: err_max, h, h_temp, s_new
  real(rp) :: r_err(6), r_temp(6)
  real(rp), parameter :: safety = 0.9_rp, p_grow = -0.2_rp
  real(rp), parameter :: p_shrink = -0.25_rp, err_con = 1.89e-4

  logical local_ref_frame

!

  h = h_try

  do

    call rkck_bmad (ele, param, r, dr_ds, s, h, r_temp, r_err, local_ref_frame)
    err_max = maxval(abs(r_err(:)/(r_scal(:)*rel_tol + abs_tol)))
    if (err_max <=  1.0) exit
    h_temp = safety * h * (err_max**p_shrink)
    h = sign(max(abs(h_temp), 0.1_rp*abs(h)), h)
    s_new = s + h

    if (s_new == s) then
      bmad_status%ok = .false.
      if (bmad_status%type_out) print *, &
                          'ERROR IN RKQS_BMAD: STEPSIZE UNDERFLOW'
      if (bmad_status%exit_on_error) call err_exit
      return
    endif

  end do

  if (err_max > err_con) then
    h_next = safety*h*(err_max**p_grow)
  else
    h_next = 5.0_rp * h
  end if

  h_did = h
  s = s+h
  r(:) = r_temp(:)

end subroutine rkqs_bmad

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

subroutine rkck_bmad (ele, param, r, dr_ds, s, h, r_out, r_err, local_ref_frame)

	implicit none

  type (ele_struct) ele
  type (lat_param_struct) param

	real(rp), intent(in) :: r(6), dr_ds(6)
	real(rp), intent(in) :: s, h
	real(rp), intent(out) :: r_out(6), r_err(6)
	real(rp) :: ak2(6), ak3(6), ak4(6), ak5(6), ak6(6), r_temp(6)
	real(rp), parameter :: a2=0.2_rp, a3=0.3_rp, a4=0.6_rp, &
    a5=1.0_rp, a6=0.875_rp, b21=0.2_rp, b31=3.0_rp/40.0_rp, &
    b32=9.0_rp/40.0_rp, b41=0.3_rp, b42=-0.9_rp, b43=1.2_rp, &
    b51=-11.0_rp/54.0_rp, b52=2.5_rp, b53=-70.0_rp/27.0_rp, &
    b54=35.0_rp/27.0_rp, &
		b61=1631.0_rp/55296.0_rp, b62=175.0_rp/512.0_rp, &
		b63=575.0_rp/13824.0_rp, b64=44275.0_rp/110592.0_rp, &
		b65=253.0_rp/4096.0_rp, c1=37.0_rp/378.0_rp, &
		c3=250.0_rp/621.0_rp, c4=125.0_rp/594.0_rp, &
		c6=512.0_rp/1771.0_rp, dc1=c1-2825.0_rp/27648.0_rp, &
		dc3=c3-18575.0_rp/48384.0_rp, dc4=c4-13525.0_rp/55296.0_rp, &
		dc5=-277.0_rp/14336.0_rp, dc6=c6-0.25_rp

  logical local_ref_frame

!

	r_temp=r+b21*h*dr_ds
	call em_field_kick(ele, param, s+a2*h, r_temp, local_ref_frame, ak2)
	r_temp=r+h*(b31*dr_ds+b32*ak2)
	call em_field_kick(ele, param, s+a3*h, r_temp, local_ref_frame, ak3) 
	r_temp=r+h*(b41*dr_ds+b42*ak2+b43*ak3)
	call em_field_kick(ele, param, s+a4*h, r_temp, local_ref_frame, ak4)
	r_temp=r+h*(b51*dr_ds+b52*ak2+b53*ak3+b54*ak4)
	call em_field_kick(ele, param, s+a5*h, r_temp, local_ref_frame, ak5)
	r_temp=r+h*(b61*dr_ds+b62*ak2+b63*ak3+b64*ak4+b65*ak5)
	call em_field_kick(ele, param, s+a6*h, r_temp, local_ref_frame, ak6)
	r_out=r+h*(c1*dr_ds+c3*ak3+c4*ak4+c6*ak6)
	r_err=h*(dc1*dr_ds+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6)

end subroutine rkck_bmad

end module
