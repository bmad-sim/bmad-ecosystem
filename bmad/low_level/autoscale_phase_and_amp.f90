!+
! Subroutine autoscale_phase_and_amp(ele, param, err_flag, scale_amp, scale_phase, call_bookkeeper)
!
! Routine to set the phase offset and amplitude scale of the accelerating field if
! this field is defined. This routine works on lcavity, rfcavity and e_gun elements.
!
! For e_gun elements there is no phase to calculate if the rf_frequency is zero.
!
! The "phase offset" is an addititive constant added to the RF phase so that ele%value(phi0$)
! is truely relative to the max accelerating phase for lcavities and relative to the accelerating
! zero crossing for rfcavities.
!
! The "amplitude scale" is a scaling factor for ele%value(voltage$) and ele%value(gradient$) so 
! that these quantities are reflect the actual on creast acceleration in volts and volts/meter.
!
! The amplitude scaling is done based upon the setting of:
!   Use the scale_amp arg if present.
!   If scale_amp arg is not present: Use ele%value(autoscale_amplitude$)
! A similar procedure is used for phase scaling.
!
! All calculations are done with a particle with the energy of the reference particle and 
! with z = 0.
!
! First: With the phase set for maximum acceleration, set the field_scale for the
! correct change in energy:
!     dE = ele%value(gradient$) * ele%value(l$) for lcavity elements.
!        = ele%value(voltage$)                  for rfcavity elements.
!
! Second:
! If the element is an lcavity then the RF phase is set for maximum acceleration.
! If the element is an rfcavity then the RF phase is set for zero acceleration and
! dE/dz will be negative (particles with z > 0 will be deaccelerated).
!
! Note: If |dE| is too small, this routine cannot scale and will do nothing.
!
! Tolerances use by the calculation are set by:
!   bmad_com        -- bmad_common_struct: Global parameters used by Bmad.
!     %autoscale_amp_abs_tol      -- Absolute amplitude tolerance. Default is 0.1 eV.
!     %autoscale_amp_rel_tol      -- Relative amplitude tolerance. Default is 1d-6.
!     %autoscale_phase_tol        -- Absolute phase tolerance. Default is 1d-5 rad/2pi.
!     %rf_phase_below_transition_ref -- Use phi0 = 0.5 as reference phase instead of 0?
!
! Input:
!   ele             -- ele_struct: RF element or e_gun.
!   param           -- lat_param_struct: lattice parameters
!   scale_amp       -- Logical, optional: Scale the amplitude? See above.
!   scale_phase     -- Logical, optional: Scale the phase? See above.
!   call_bookkeeper -- Logical, optional: Call lattice_bookkeeper at end? Default is True.
!
! Output:
!   ele      -- ele_struct: element with phase and amplitude adjusted. 
!   err_flag -- Logical, Set true if there is an error. False otherwise.
!-

subroutine autoscale_phase_and_amp(ele, param, err_flag, scale_phase, scale_amp, call_bookkeeper)

use super_recipes_mod, only: super_zbrent
use bmad_interface, dummy => autoscale_phase_and_amp

implicit none

type (ele_struct), target :: ele
type (lat_param_struct), target :: param
type (coord_struct) orbit0
type (em_field_struct) field1, field2
integer, parameter :: n_sample = 16

real(rp) pz, phi, pz_max, phi_max, e_tot, scale_correct, dE_peak_wanted, dE_cut
real(rp) dphi, e_tot_start, pz_plus, pz_minus, b, c, phi_tol, scale_tol, phi_max_old
real(rp) value_saved(num_ele_attrib$), phi0_autoscale_original, pz_arr(0:n_sample-1), pz_max1, pz_max2
real(rp) dE_max1, dE_max2, integral, int_tot, int_old, s

integer i, j, tracking_method_saved, num_times_lost, i_max1, i_max2
integer n_pts, n_pts_tot, n_loop, n_loop_max, status, sign_of_dE, tm
integer :: n_call ! Used for debugging.

logical :: is_lost
logical step_up_seen, err_flag, do_scale_phase, do_scale_amp, phase_scale_good, amp_scale_good
logical, optional :: scale_phase, scale_amp, call_bookkeeper
logical :: debug = .false.

character(*), parameter :: r_name = 'autoscale_phase_and_amp'

! Check if auto scale is needed.

err_flag = .false.
if (.not. ele%is_on) return

do_scale_phase = logic_option(is_true(ele%value(autoscale_phase$)), scale_phase)
do_scale_amp   = logic_option(is_true(ele%value(autoscale_amplitude$)), scale_amp)

if (ele%key == e_gun$ .and. ele%value(rf_frequency$) == 0) then
  do_scale_phase = .false.
endif

if (.not. do_scale_phase .and. .not. do_scale_amp) return

! Init.
! Note: phi0_autoscale is set in neg_pz_calc

! Autoscaling with linear tracking does not make sense so just do nothing.

tm = ele%tracking_method
if (tm == mad$) return
if (tm == linear$) return

! bmad_standard just needs to set e_tot$, p0c$, and phi0_autoscale$
! Zero length rfcavity elements use bmad_standard tracking

if (tm == bmad_standard$ .or. (ele%key == rfcavity$ .and. ele%value(l$) == 0 .and. (tm == runge_kutta$ .or. &
          tm == fixed_step_runge_kutta$ .or. tm == time_runge_kutta$ .or. tm == fixed_step_time_runge_kutta$))) then
  if (ele%key == lcavity$) then 
    ! Set e_tot$ and p0c$ 
    phi = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$)) 
    e_tot = ele%value(e_tot_start$) + ele%value(gradient$) * ele%value(l$) * cos(phi)
    call convert_total_energy_to (e_tot, param%particle, pc = ele%value(p0c$), err_flag = err_flag, print_err = .false.)
    if (err_flag) then
      call out_io (s_error$, r_name, 'REFERENCE ENERGY BELOW REST MASS AT EXIT END OF LCAVITY: ' // ele_full_name(ele))
      ! Unstable_factor is formulated to be usable for optimization when the lattice is not stable.
      param%unstable_factor = ele%ix_ele - e_tot / mass_of(param%particle)
      return
    endif
    ele%value(e_tot$) = e_tot
  endif 
  
  ele%value(phi0_autoscale$) = 0
  return
endif

phi0_autoscale_original = ele%value(phi0_autoscale$)

! Compute Energy gain at peak (zero phase)

select case (ele%key)
case (rfcavity$)
  dE_peak_wanted = ele%value(voltage$)
  e_tot_start = ele%value(e_tot$)
case (lcavity$)
  dE_peak_wanted = ele%value(gradient$) * ele%value(l$)
  e_tot_start = ele%value(e_tot_start$)
case (e_gun$)
  dE_peak_wanted = ele%value(gradient$) * ele%value(l$)
  e_tot_start = ele%value(e_tot_ref_init$)
case default
  call out_io (s_fatal$, r_name, 'CONFUSED ELEMENT TYPE!')
  if (global_com%exit_on_error) call err_exit ! exit on error.
  return
end select

sign_of_dE = sign_of(dE_peak_wanted)

! Auto scale amplitude when dE_peak_wanted is zero or very small is not possible.
! Therefore if dE_peak_wanted is less than dE_cut then do nothing.

! scale_tol is the tolerance for scale_correct.
! scale_tol = E_tol / dE_peak_wanted corresponds to a tolerance in dE_peak_wanted of E_tol. 

if (do_scale_amp) then
  dE_cut = 10 ! eV
  if (abs(dE_peak_wanted) < dE_cut) return
  scale_tol = max(bmad_com%autoscale_amp_rel_tol, bmad_com%autoscale_amp_abs_tol / abs(dE_peak_wanted))
else
  if (dE_peak_wanted == 0) sign_of_dE = 1    ! Assume want accelerating.
  scale_tol = 1
endif

if (ele%value(field_autoscale$) == 0) then
  ! Cannot autophase if not allowed to make the ele%value(field_autoscale$) non-zero.
  if (.not. do_scale_amp) then
    call out_io (s_fatal$, &
            r_name, 'CANNOT AUTO PHASE IF NOT ALLOWED TO MAKE THE FIELD_SCALE NON-ZERO FOR: ' // ele_full_name(ele))
    if (global_com%exit_on_error) call err_exit ! exit on error.
    return 
  endif
  ele%value(field_autoscale$) = 1.0_rp  ! Initial guess.
endif

! scale_correct is the correction factor applied to ele%value(field_autoscale$) on each iteration:
!  ele%value(field_autoscale$)(new) = ele%value(field_autoscale$)(old) * scale_correct

phi_tol = bmad_com%autoscale_phase_tol

!------------------------------------------------------
! zero frequency e_gun

if (ele%key == e_gun$ .and. ele%value(rf_frequency$) == 0) then
  tracking_method_saved = ele%tracking_method
  if (ele%tracking_method == taylor$) ele%tracking_method = symp_lie_ptc$
  value_saved = ele%value

  if (ele%tracking_method == bmad_standard$) ele%tracking_method = time_runge_kutta$
  ele%value(gradient_err$) = 0

  pz_max = pz_calc(0.0_rp, err_flag)
  if (err_flag) return
  if (is_lost) ele%value(field_autoscale$) = -ele%value(field_autoscale$) ! Maybe field in wrong direction?

  do i = 1, 100
    pz_max = pz_calc(0.0_rp, err_flag)
    if (err_flag) return
    scale_correct = dE_peak_wanted / dE_particle(pz_max)
    if (abs(scale_correct) > 1000) scale_correct = scale_correct / 10
    if (abs(scale_correct) > 1000) scale_correct = sign(1000.0_rp, scale_correct)
    ele%value(field_autoscale$) = ele%value(field_autoscale$) * scale_correct
    if (abs(scale_correct - 1) < scale_tol) exit
  enddo

  if (i == 101) then
    call out_io (s_fatal$, r_name, 'CANNOT FIND CORRECT AMPLITUDE SCALE FOR: ' // ele_full_name(ele))
    if (global_com%exit_on_error) call err_exit ! exit on error.
  endif

  ele%tracking_method = tracking_method_saved
  ele%value = value_saved

  if (logic_option(.true., call_bookkeeper) .and. ele%ix_ele > 0) call lattice_bookkeeper(ele%branch%lat)
  return
endif

!------------------------------------------------------
! Set error fields to zero

value_saved = ele%value
ele%value(phi0$) = 0
ele%value(phi0_multipass$) = 0
ele%value(phi0_err$) = 0
if (ele%key == lcavity$ .or. ele%key == e_gun$) ele%value(gradient_err$) = 0

tracking_method_saved = ele%tracking_method

if (ele%tracking_method == taylor$) ele%tracking_method = symp_lie_ptc$

if (ele%tracking_method == bmad_standard$) then
  if (ele%key == e_gun$) then
    ele%tracking_method = time_runge_kutta$
  else
    ele%tracking_method = runge_kutta$
  endif
endif

phi_max = ele%value(phi0_autoscale$)   ! Init guess
if (ele%key == rfcavity$) phi_max = ele%value(phi0_max$)

phi_max_old = 100 ! Number far from unity

! See if ele%value(phi0_autoscale$) and ele%value(field_autoscale$) are already set correctly.
! If so we can quit.

phase_scale_good = .true.
amp_scale_good = .true. 

pz_max = pz_calc(phi_max, err_flag)
if (err_flag) return

if (.not. is_lost) then
  if (do_scale_phase) then
    pz_plus  = pz_calc(phi_max + 2 * phi_tol, err_flag); if (err_flag) return
    pz_minus = pz_calc(phi_max - 2 * phi_tol, err_flag); if (err_flag) return
    phase_scale_good = (pz_max >= pz_plus .and. pz_max >= pz_minus)
  endif

  if (do_scale_amp) then
    scale_correct = dE_peak_wanted / dE_particle(pz_max) 
    amp_scale_good = (abs(scale_correct - 1) < 2 * scale_tol)
  endif

  if (phase_scale_good .and. amp_scale_good) then
    ele%value(phi0_autoscale$) = phi0_autoscale_original
    call cleanup_this()
    return
  endif
endif

! The ele%value(field_autoscale$) may be orders of magnitude off so do an initial guess
! based upon the integral of abs(voltage) through the element.

if (do_scale_amp .and. ele%value(l$) /= 0) then
  n_pts = 1
  int_tot = 0
  n_pts_tot = 0

  do 
    integral = 0
    do i = 1, n_pts
      s = ele%value(l$) * (2*i - 1.0) / (2*n_pts)
      ! Sample field at two phases and take the max. This is crude but effective.
      ele%value(phi0_autoscale$) = 0
      call em_field_calc (ele, param, s, orbit0, .true., field1)
      ele%value(phi0_autoscale$) = pi/2
      call em_field_calc (ele, param, s, orbit0, .true., field2)
      integral = integral + max(abs(field1%e(3)), abs(field2%e(3))) * ele%value(l$) / n_pts
    enddo

    n_pts_tot = n_pts_tot + n_pts
    int_old = int_tot
    int_tot = ((n_pts_tot - n_pts) * int_tot + n_pts * integral) / n_pts_tot
    if (n_pts_tot > 16) then
      if (int_tot == 0 .and. int_old == 0) then
        call out_io (s_error$, r_name, 'FIELD IS ZERO FOR: ' // ele_full_name(ele)) 
        call cleanup_this()
        return
      endif

      if (abs(int_tot - int_old) <= 0.2 * (int_tot + int_old)) then
        ele%value(field_autoscale$) = ele%value(field_autoscale$) * abs(dE_peak_wanted) / integral
        exit
      endif
    endif

    n_pts = 2 * n_pts

  enddo
endif

! OK so the input ele%value(phi0_autoscale$) or ele%value(field_autoscale$) are not set correctly...
! First choose a starting phi_max by finding an approximate phase for max acceleration.
! We start by testing n_sample phases.
! pz_max1 gives the maximal acceleration. pz_max2 gives the second largest.

pz_arr(0) = pz_max
dphi = 1.0_rp / n_sample

do i = 1, n_sample - 1
  pz_arr(i) = pz_calc(phi_max + i*dphi, err_flag); if (err_flag) return
enddo

i_max1 = maxloc(sign_of_dE*pz_arr, 1) - 1
pz_max1 = pz_arr(i_max1)
dE_max1 = dE_particle(pz_max1)

pz_arr(i_max1) = -sign_of_dE  ! To find next max
i_max2 = maxloc(sign_of_dE*pz_arr, 1) - 1
pz_max2 = pz_arr(i_max2)
dE_max2 = dE_particle(pz_max2)

! If we do not have any phase that shows acceleration this generally means that the
! initial particle energy is low and the ele%value(field_autoscale$) is much to large.

if (sign_of_dE*dE_max1 <= 0) then
  call out_io (s_error$, r_name, 'CANNOT FIND ACCELERATING PHASE REGION FOR: ' // ele_full_name(ele))
  err_flag = .true.
  return
endif

! If dE_max1 is large compared to dE_max2 then just use the dE_max1 phase. 
! Otherwise take half way between dE_max1 and dE_max2 phases.

if (2*abs(dE_max2) < abs(dE_max1)) then  ! Just use dE_max1 point
  phi_max = phi_max + dphi * i_max1
  pz_max = pz_max1
! wrap around case when i_max1 = 0 and i_max2 = n_sample-1 or vice versa.
elseif (2*abs(i_max1 - i_max2) > n_sample) then   
  phi_max = phi_max + dphi * (i_max1 + i_max2 - n_sample) / 2.0
  pz_max = pz_calc(phi_max, err_flag); if (err_flag) return
else
  phi_max = phi_max + dphi * (i_max1 + i_max2) / 2.0
  pz_max = pz_calc(phi_max, err_flag); if (err_flag) return
endif

! Now adjust ele%value(field_autoscale$) for the correct acceleration at the phase for maximum acceleration. 

n_call = 0  ! For debug purposes.
n_loop_max = 10000
num_times_lost = 0
dphi = 0.05

if (debug) print *, '------------------------------------------'
if (debug) print *, 'TOL: ', scale_tol, phi_tol

main_loop: do n_loop = 1, n_loop_max

  ! Find approximately the phase for maximum acceleration.
  ! First go in +phi direction until pz decreases.

  step_up_seen = .false.

  do i = 1, 100
    phi = phi_max + dphi
    pz = pz_calc(phi, err_flag); if (err_flag) return
    if (debug) print *, 'FWD:', i, trim(ele%name), '  ', phi, pz

    if (is_lost) then
      do j = -19, 20
        if (debug) print *, j, phi_max+j/40.0, pz_calc(phi_max + j / 40.0, err_flag)
      enddo
      call out_io (s_error$, r_name, 'CANNOT STABLY TRACK PARTICLE FOR ELEMENT: ' // ele_full_name(ele))
      err_flag = .true.
      return
    endif

    if (sign_of_dE*pz < sign_of_dE*pz_max) then
      pz_plus = pz
      exit
    endif

    pz_minus = pz_max
    pz_max = pz
    phi_max = phi
    step_up_seen = .true.
  enddo

  ! If needed: Now go in -phi direction until pz decreases

  if (.not. step_up_seen) then
    do
      phi = phi_max - dphi
      pz = pz_calc(phi, err_flag); if (err_flag) return
      if (debug) print *, 'REV:', i, trim(ele%name), '  ', phi, pz
      if (sign_of_dE*pz < sign_of_dE*pz_max) then
        pz_minus = pz
        exit
      endif
      pz_plus = pz_max
      pz_max = pz
      phi_max = phi
    enddo
  endif

  ! Quadradic interpolation to get the maximum phase.
  ! Formula: pz = a + b*dt + c*dt^2 where dt = (phi-phi_max) / dphi

  b = (pz_plus - pz_minus) / 2
  c = pz_plus - pz_max - b

  phi_max = phi_max - b * dphi / (2 * c)
  pz_max = pz_calc(phi_max, err_flag); if (err_flag) return

  if (debug) print '(a, 3es18.10, f10.6)', 'MAX:', phi_max, pz_max, ele%value(field_autoscale$), dphi

  ! Now scale ele%value(field_autoscale$)
  ! scale_correct = dE(design) / dE (from tracking)
  ! Can overshoot so if scale_correct is too large then scale back by a factor of 10

  if (do_scale_amp) then
    scale_correct = dE_peak_wanted / dE_particle(pz_max)
    if (scale_correct > 1000) scale_correct = max(1000.0_rp, scale_correct / 10)
    ele%value(field_autoscale$) = ele%value(field_autoscale$) * scale_correct
  else
    scale_correct = 1
  endif

  if (debug) print '(a, i5, 2es12.3)', 'TEST:', n_loop, abs(scale_correct - 1), abs(phi_max-phi_max_old) 

  if (abs(scale_correct - 1) < scale_tol .and. abs(phi_max-phi_max_old) < phi_tol) exit
  phi_max_old = phi_max

  dphi = 0.05
  if (abs(scale_correct - 1) < 0.1) dphi = max(phi_tol, 0.1*sqrt(2*abs(scale_correct - 1))/twopi)

  if (do_scale_phase) then
    pz_max = pz_calc(phi_max, err_flag); if (err_flag) return
  endif

  if (n_loop == n_loop_max) then
    call out_io (s_warn$, r_name, 'AUTO SCALING NOT CONVERGING FOR ELEMENT: ' // ele_full_name(ele))
  endif

enddo main_loop

! For an rfcavity now find the zero crossing with negative slope which is
! about 90deg away from max acceleration.

if (ele%key == rfcavity$) then
  value_saved(phi0_max$) = ele%value(phi0_autoscale$)  ! Save for use with OPAL
  if (do_scale_phase) then
    dphi = 0.1
    if (bmad_com%rf_phase_below_transition_ref) dphi = -dphi

    phi_max = phi_max - dphi
    do
      phi = phi_max - dphi
      pz = pz_calc(phi, err_flag); if (err_flag) return
      if (pz < 0) exit
      phi_max = phi
    enddo
    ele%value(phi0_autoscale$) = modulo2(super_zbrent(neg_pz_calc, phi_max-dphi, phi_max, &
                                                                1e-15_rp, 1d-9, status), 0.5_rp)
  endif
endif

! Cleanup

call cleanup_this()

if (associated (ele%branch)) then
  if (do_scale_amp)   call set_flags_for_changed_attribute (ele, ele%value(field_autoscale$))
  if (do_scale_phase) call set_flags_for_changed_attribute (ele, ele%value(phi0_autoscale$))
endif

if (logic_option(.true., call_bookkeeper) .and. ele%ix_ele > 0) call lattice_bookkeeper(ele%branch%lat)

!------------------------------------
contains

subroutine cleanup_this ()

value_saved(field_autoscale$) = ele%value(field_autoscale$) 
value_saved(phi0_autoscale$) = ele%value(phi0_autoscale$)

ele%value = value_saved
if (.not. do_scale_phase) ele%value(phi0_autoscale$) = phi0_autoscale_original

ele%tracking_method = tracking_method_saved

end subroutine cleanup_this

!----------------------------------------------------------------
! contains

function pz_calc (phi, err_flag) result (pz)

implicit none

type (coord_struct) start_orb, end_orb
real(rp), intent(in) :: phi
real(rp) pz
logical err_flag

! 


time_runge_kutta_com%print_too_many_step_err = .false.
ele%value(phi0_autoscale$) = phi
call attribute_bookkeeper(ele, .true.)
if (ele%tracking_method == linear$) call make_mat6(ele, param)

call init_coord (start_orb, ele = ele, element_end = upstream_end$)
call track1 (start_orb, ele, param, end_orb, err_flag = err_flag, ignore_radiation = .true.)
time_runge_kutta_com%print_too_many_step_err = .true.

pz = end_orb%vec(6)
is_lost = .not. particle_is_moving_forward(end_orb)
if (is_lost) pz = -2

n_call = n_call + 1

end function pz_calc

!------------------------------------
! contains
! Function returns the energy gain of a particle given final pz

function dE_particle(pz) result (de)

real(rp) pz, e_tot, de

!

if (pz < -1) then
  e_tot = -1  ! If particle has been lost.
else
  call convert_pc_to ((1 + pz) * ele%value(p0c$), param%particle, e_tot = e_tot)
endif

de = e_tot - e_tot_start

end function dE_particle

!----------------------------------------------------------------
! contains

function neg_pz_calc (phi, status) result (neg_pz)

implicit none

real(rp), intent(in) :: phi
real(rp) neg_pz
integer status
logical err_flag

! brent finds minima so need to flip the final energy

neg_pz = -pz_calc(phi, err_flag)

end function neg_pz_calc

end subroutine autoscale_phase_and_amp
