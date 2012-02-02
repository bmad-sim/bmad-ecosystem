module rf_mod

use runge_kutta_mod

real(rp), pointer, private :: field_scale, dphi0_ref
type (lat_param_struct), pointer, private :: param_com
type (ele_struct), pointer, private :: ele_com

integer, private, save :: n_loop ! Used for debugging.

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Subroutine rf_auto_scale_phase_and_amp(ele, param)
!
! Routine to set the reference phase and amplitude of the accelerating rf field mode if
! this mode is present. The accelerating mode is defined to be ele%em_field%mode(1) if
! mode(1)%m = 0. 
!
! All calculations are done with a particle with the energy of the reference particle and 
! with z = 0.
!
! First: With the phase set for maximum acceleration, set the field_scale for the
! correct acceleration:
!     acceleration = ele%value(gradient$) * ele%value(l$) for lcavity elements.
!                  = ele%value(voltage$)                  for rfcavity elements.
!
! Second:
! If the element is an lcavity then the RF phase is set for maximum acceleration.
! If the element is an rcavity then the RF phase is set for zero acceleration and
! dE/dz will be negative (particles with z > 0 will be deaccelerated).
!
! Modules needed
!   use rf_mod
!
! Input:
!   ele   -- ele_struct: RF element. Either lcavity or rfcavity.
!     %value(gradient$) -- Accelerating gradient to match to if an lcavity.
!     %value(voltage$)  -- Accelerating voltage to match to if an rfcavity.
!   param -- lat_param_struct: lattice parameters
!
! Output:
!   ele -- ele_struct: element with phase and amplitude adjusted.
!-

subroutine rf_auto_scale_phase_and_amp(ele, param)

use super_recipes_mod
use nr, only: zbrent

implicit none

type (ele_struct), target :: ele
type (lat_param_struct), target :: param

real(rp) pz, phi, pz_max, phi0, phi_max, e_tot, f_correct, wanted_de
real(rp) dphi, e_tot_start, pz_plus, pz_minus, b, c, phi_tol, pz_tol, phi_max_old
real(rp) value_saved(n_attrib_maxx)

integer i, tracking_method_saved

logical step_up_seen

! Check if auto scale is needed

if (associated (ele%lat)) then
  if (.not. ele%lat%rf_auto_scale_phase_and_amp) return
endif

if (ele%tracking_method == bmad_standard$ .or. ele%tracking_method == mad$) return

! Init.
! Note: dphi0_ref is set in neg_pz_calc

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

if (.not. associated(field_scale)) return

! Compute wanted energy gain. 
! Autophasing when the wanted energy gain is zero or very small is not possible
! Therefore if the wanted energy gain is less than 1 eV then do nothing.

ele_com => ele
param_com => param

if (.not. ele%is_on) return
select case (ele%key)
case (rfcavity$)
  wanted_de = ele%value(voltage$)
  e_tot_start = ele%value(e_tot$)
case (lcavity$)
  wanted_de = (ele%value(gradient$) + ele%value(gradient_err$) + &
                            gradient_shift_sr_wake (ele, param)) * ele%value(l$)
                                        
  e_tot_start = ele%value(e_tot_start$)
case default
  if (bmad_status%exit_on_error) call err_exit ! exit on error.
end select

if (abs(wanted_de) < 1) return

! Set error fields to zero

value_saved = ele%value
ele%value(phi0$) = 0
ele%value(dphi0$) = 0
ele%value(phi0_err$) = 0
if (ele%key == lcavity$) ele%value(gradient_err$) = 0

tracking_method_saved = ele%tracking_method
if (ele%tracking_method == bmad_standard$) ele%tracking_method = runge_kutta$

phi0 = dphi0_ref
phi_max = dphi0_ref   ! Init guess
if (ele%key == rfcavity$) phi_max = ele%value(dphi0_max$)

phi_max_old = 100 ! Number far from unity
dphi = 0.05
phi_tol = 1d-5
pz_tol = 1d-7

! See if %dphi0_ref and %field_scale are already set correctly

pz_plus  = -neg_pz_calc(phi_max + 2 * phi_tol)
pz_minus = -neg_pz_calc(phi_max - 2 * phi_tol)
pz_max = -neg_pz_calc(phi_max)

call convert_pc_to ((1 + pz_max) * ele%value(p0c$), param%particle, e_tot = e_tot)
f_correct = wanted_de / (e_tot - e_tot_start)

if (pz_max > pz_plus .and. pz_max > pz_minus .and. abs(f_correct - 1) < 2 * pz_tol) then
  call cleanup_this()
  dphi0_ref = phi0
  return
endif

! Now adjust %field_scale for the correct acceleration at the phase for maximum accelleration. 

n_loop = 0  ! For debug purposes.

main_loop: do

  ! Find approximately the phase for maximum acceleration.
  ! First go in +phi direction until pz decreases.

  step_up_seen = .false.
  do i = 1, 10
    phi = phi_max + dphi
    pz = -neg_pz_calc(phi)
    if (pz < pz_max) exit
    pz_max = pz
    phi_max = phi
    step_up_seen = .true.
    if (i == 10) then  ! field too strong and always loosing particles
      field_scale = field_scale / 10
      pz = -neg_pz_calc(phi)
      cycle main_loop
    endif
  enddo

  pz_plus = pz

  ! If needed: Now go in -phi direction until pz decreases

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

  ! Now scale %field_scale
  ! f_correct = dE(design) / dE (from tracking)
  ! Can overshoot so if f_correct is too large then scale back by a factor of 10

  call convert_pc_to ((1 + pz_max) * ele%value(p0c$), param%particle, e_tot = e_tot)
  f_correct = wanted_de / (e_tot - e_tot_start)
  if (f_correct > 1000) f_correct = max(1000.0_rp, f_correct / 10)
  field_scale = field_scale * f_correct

  if (abs(f_correct - 1) < pz_tol .and. abs(phi_max-phi_max_old) < phi_tol) exit
  phi_max_old = phi_max

  dphi = 0.05
  if (abs(f_correct - 1) < 0.1) dphi = max(phi_tol, 0.1*sqrt(2*abs(f_correct - 1))/twopi)

  pz_max = -neg_pz_calc(phi_max)

enddo main_loop

! For an rfcavity now find the zero crossing with negative slope which is
! about 90deg away from max acceleration.

if (ele%key == rfcavity$) then
  ele%value(dphi0_max$) = dphi0_ref  ! Save for use with OPAL
  dphi = 0.1
  do
    phi = phi_max + dphi
    pz = -neg_pz_calc(phi)
    if (pz < 0) exit
    phi_max = phi
  enddo
  dphi0_ref = modulo2 (zbrent(neg_pz_calc, phi_max, phi_max+dphi, 1d-9), 0.5_rp)
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

n_loop = n_loop + 1

end function

end module
