module rf_mod

use runge_kutta_mod

real(rp), pointer, private :: field_scale, dphi0_ref
type (lat_param_struct), pointer, private :: param_com
type (ele_struct), pointer, private :: ele_com

integer, private, save :: n_loop ! Used for debugging.
logical, private, save :: is_lost

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Subroutine rf_auto_scale_phase_and_amp(ele, param, err_flag)
!
! Routine to set the reference phase and amplitude of the accelerating rf field mode if
! this mode is present. The accelerating mode is defined to be ele%em_field%mode(1) if
! mode(1)%m = 0. 
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
! Note: For an rfcavity if ele%lat%rf_auto_scale_phase = F and lat%rf_auto_scale_amp = T then
! the first step is done above but then the phase is reset to the original phase in step 2.
!
! Modules needed
!   use rf_mod
!
! Input:
!
!   ele   -- ele_struct: RF element. Either lcavity or rfcavity.
!     %value(gradient$) -- Accelerating gradient to match to if an lcavity.
!     %value(voltage$)  -- Accelerating voltage to match to if an rfcavity.
!     %lat%rf_auto_scale_phase ! Scale phase? Default is True if ele%lat is not associated.
!     %lat%rf_auto_scale_amp   ! Scale amp?   Default is True if ele%lat is not associated.
!   param -- lat_param_struct: lattice parameters
!
! Output:
!   ele      -- ele_struct: element with phase and amplitude adjusted.
!   err_flag -- Logical, Set true if there is an error. False otherwise.
!-

subroutine rf_auto_scale_phase_and_amp(ele, param, err_flag)

use super_recipes_mod
use nr, only: zbrent

implicit none

type (ele_struct), target :: ele
type (lat_param_struct), target :: param

real(rp) pz, phi, pz_max, phi_max, e_tot, scale_correct, dE_peak, dE_cut, E_tol
real(rp) dphi, e_tot_start, pz_plus, pz_minus, b, c, phi_tol, scale_tol, phi_max_old
real(rp) value_saved(n_attrib_maxx), dphi0_ref_original

integer i, tracking_method_saved, num_times_lost

logical step_up_seen, err_flag, has_been_lost, scale_phase, scale_amp, adjust_phase, adjust_amp

character(28), parameter :: r_name = 'rf_auto_scale_phase_and_amp'

! Check if auto scale is needed.
! Adjust_phase is used to determine if the phase can be adjusted when scaling the amplitude.

err_flag = .false.

scale_phase = .true.
scale_amp   = .true.
if (associated (ele%lat)) then
  scale_phase = ele%lat%rf_auto_scale_phase
  scale_amp   = ele%lat%rf_auto_scale_amp
  if (.not. scale_phase .and. .not. scale_amp) return
endif

adjust_phase = (scale_phase .or. ele%key == rfcavity$)
adjust_amp   = scale_amp

! Init.
! Note: dphi0_ref is set in neg_pz_calc

if (ele%tracking_method == bmad_standard$ .or. ele%tracking_method == mad$) return

nullify(field_scale)

select case (ele%field_calc)
case (bmad_standard$) 
  field_scale => ele%value(field_scale$)
  dphi0_ref => ele%value(dphi0_ref$)
case (grid$, map$, custom$)
  do i = 1, size(ele%em_field%mode)
    if (ele%em_field%mode(i)%freq /= 0 .and. ele%em_field%mode(i)%m == 0) then
      field_scale => ele%em_field%mode(i)%field_scale
      dphi0_ref => ele%em_field%mode(i)%dphi0_ref
      exit
    endif
  enddo
end select

dphi0_ref_original = dphi0_ref

if (.not. associated(field_scale)) return

! Compute Energy gain at peak (zero phase)

ele_com => ele
param_com => param

if (.not. ele%is_on) return
select case (ele%key)
case (rfcavity$)
  dE_peak = ele%value(voltage$)
  e_tot_start = ele%value(e_tot$)
case (lcavity$)
  dE_peak = ele%value(gradient$) + ele%value(l$)
  e_tot_start = ele%value(e_tot_start$)
case default
  call out_io (s_fatal$, r_name, 'CONFUSED ELEMENT TYPE!')
  if (bmad_status%exit_on_error) call err_exit ! exit on error.
end select

! Autophasing when dE_peak is zero or very small is not possible.
! Therefore if dE_peak is less than dE_cut then do nothing.

dE_cut = 10 ! eV
if (abs(dE_peak) < dE_cut) return

! Set error fields to zero

value_saved = ele%value
ele%value(phi0$) = 0
ele%value(dphi0$) = 0
ele%value(phi0_err$) = 0
if (ele%key == lcavity$) ele%value(gradient_err$) = 0

tracking_method_saved = ele%tracking_method
if (ele%tracking_method == bmad_standard$) ele%tracking_method = runge_kutta$

phi_max = dphi0_ref   ! Init guess
if (ele%key == rfcavity$) phi_max = ele%value(dphi0_max$)

phi_max_old = 100 ! Number far from unity
dphi = 0.05

! scale_correct is the correction factor applied to field_scale on each iteration:
!  field_scale(new) = field_scale(old) * scale_correct
! scale_tol is the tolerance for scale_correct.
! scale_tol = E_tol / dE_peak corresponds to a tolerance in dE_peak of E_tol. 

E_tol = 0.1 ! eV
scale_tol = max(1d-7, E_tol / dE_peak) ! tolerance for scale_correct
phi_tol = 1d-5

! See if %dphi0_ref and %field_scale are already set correctly

pz_max   = -neg_pz_calc(phi_max)

if (adjust_phase) then
  pz_plus  = -neg_pz_calc(phi_max + 2 * phi_tol)
  pz_minus = -neg_pz_calc(phi_max - 2 * phi_tol)
else
  pz_plus  = -100
  pz_minus = -100
endif

if (adjust_amp) then
  call convert_pc_to ((1 + pz_max) * ele%value(p0c$), param%particle, e_tot = e_tot)
  scale_correct = dE_peak / (e_tot - e_tot_start)
else
  scale_correct = 1
endif

if (pz_max > pz_plus .and. pz_max > pz_minus .and. abs(scale_correct - 1) < 2 * scale_tol) then
  call cleanup_this()
  dphi0_ref = dphi0_ref_original
  return
endif

! Now adjust %field_scale for the correct acceleration at the phase for maximum accelleration. 

n_loop = 0  ! For debug purposes.
num_times_lost = 0

main_loop: do

  ! Find approximately the phase for maximum acceleration.
  ! First go in +phi direction until pz decreases.

  if (adjust_phase) then
    has_been_lost = .false.
    step_up_seen = .false.

    do i = 1, 100
      phi = phi_max + dphi
      pz = -neg_pz_calc(phi)
      if (is_lost) then  ! field too strong so reduce
        field_scale = field_scale / 10
        pz_max = -neg_pz_calc(phi_max)
        step_up_seen = .false.
        has_been_lost = .true.
        cycle
      endif
      if (pz < pz_max) exit
      pz_max = pz
      phi_max = phi
      step_up_seen = .true.
    enddo

    ! Can get into an infinite loop where, if the starting energy is too low,
    ! the algorithm cannot converge. In this case bail out.

    if (has_been_lost) num_times_lost = num_times_lost + 1
    if (num_times_lost == 3) then
      call out_io (s_error$, r_name, 'CANNOT STABLY TRACK PARTICLE!', &
                                     '[Can happen when particle energy is near E_rest_mass]')
      err_flag = .true.
      return
    endif

    ! If needed: Now go in -phi direction until pz decreases

    pz_plus = pz
    if (.not. step_up_seen) then
      do
        phi = phi_max - dphi
        pz = -neg_pz_calc(phi)
        if (pz < pz_max) exit
        pz_max = pz
        phi_max = phi
      enddo
    endif

    pz_minus = pz

    ! Quadradic interpolation to get the maximum phase.
    ! Formula: pz = a + b*dt + c*dt^2 where dt = (phi-phi_max) / dphi

    b = (pz_plus - pz_minus) / 2
    c = pz_plus - pz_max - b

    phi_max = phi_max - b * dphi / (2 * c)
    pz_max = -neg_pz_calc(phi_max)

  endif

  ! Now scale %field_scale
  ! scale_correct = dE(design) / dE (from tracking)
  ! Can overshoot so if scale_correct is too large then scale back by a factor of 10

  if (adjust_amp) then
    call convert_pc_to ((1 + pz_max) * ele%value(p0c$), param%particle, e_tot = e_tot)
    scale_correct = dE_peak / (e_tot - e_tot_start)
    if (scale_correct > 1000) scale_correct = max(1000.0_rp, scale_correct / 10)
    field_scale = field_scale * scale_correct
  else
    scale_correct = 1
  endif

  if (abs(scale_correct - 1) < scale_tol .and. abs(phi_max-phi_max_old) < phi_tol) exit
  phi_max_old = phi_max

  dphi = 0.05
  if (abs(scale_correct - 1) < 0.1) dphi = max(phi_tol, 0.1*sqrt(2*abs(scale_correct - 1))/twopi)

  if (adjust_phase) then
    pz_max = -neg_pz_calc(phi_max)
  endif

enddo main_loop

! For an rfcavity now find the zero crossing with negative slope which is
! about 90deg away from max acceleration.

if (ele%key == rfcavity$) then
  ele%value(dphi0_max$) = dphi0_ref  ! Save for use with OPAL
  if (scale_phase) then
    dphi = 0.1
    do
      phi = phi_max + dphi
      pz = -neg_pz_calc(phi)
      if (pz < 0) exit
      phi_max = phi
    enddo
    dphi0_ref = modulo2 (zbrent(neg_pz_calc, phi_max, phi_max+dphi, 1d-9), 0.5_rp)
  else
    dphi0_ref = dphi0_ref_original
  endif
endif

! Cleanup

call cleanup_this()

!------------------------------------
contains

subroutine cleanup_this ()

select case (ele%field_calc)
case (bmad_standard$) 
  value_saved(field_scale$) = field_scale 
  value_saved(dphi0_ref$) = dphi0_ref
end select

ele%value = value_saved
ele%tracking_method = tracking_method_saved

end subroutine cleanup_this

end subroutine rf_auto_scale_phase_and_amp

!----------------------------------------------------------------

function neg_pz_calc (phi) result (neg_pz)

implicit none

type (coord_struct) start_orb, end_orb
real(rp), intent(in) :: phi
real(rp) neg_pz

! brent finds minima so need to flip the final energy

dphi0_ref = phi
call track1 (start_orb, ele_com, param_com, end_orb)

neg_pz = -end_orb%vec(6)
if (param_com%lost) neg_pz = 1

is_lost = param_com%lost

n_loop = n_loop + 1

end function

end module
