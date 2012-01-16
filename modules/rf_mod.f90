module rf_mod

use runge_kutta_mod

real(rp), pointer, private :: field_scale, phi0_ref
type (lat_param_struct), pointer, private :: param_com
type (ele_struct), pointer, private :: ele_com

integer, private :: n_loop

contains

!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------------
!+
! Subroutine rf_accel_mode_adjust_phase_and_amp (ele, param)
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

subroutine rf_accel_mode_adjust_phase_and_amp (ele, param)

use super_recipes_mod
use nr, only: zbrent

implicit none

type (ele_struct), target :: ele
type (lat_param_struct), target :: param

real(rp) pz, theta, pz_max, theta0, theta_max, e_tot, f_correct, wanted_de
real(rp) dtheta, e_tot_start, pz_plus, pz_minus, b, c, theta_tol, pz_tol, theta_max_old
real(rp) phi0_saved, dphi0_saved, phi0_err_saved

integer i, tracking_method_saved

logical step_up_seen

! Init

if (ele%tracking_method == bmad_standard$ .or. ele%tracking_method == mad$) return

nullify(field_scale)

select case (ele%field_calc)
case (bmad_standard$) 
  field_scale => ele%value(field_scale$)
  phi0_ref => ele%value(phi0_ref$)
case (grid$, map$, custom$)
  do i = 1, size(ele%em_field%mode)
    if (ele%em_field%mode(i)%freq /= 0 .and. ele%em_field%mode(i)%m == 0) then
      field_scale => ele%em_field%mode(i)%field_scale
      phi0_ref => ele%em_field%mode(i)%phi0_ref
      exit
    endif
  enddo
end select

if (.not. associated(field_scale)) return

! 

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
  call err_exit ! exit on error.
end select

n_loop = 0  ! For debug purposes.

if (wanted_de == 0) then
  field_scale = 0
  return
endif

phi0_saved     = ele%value(phi0$) 
dphi0_saved    = ele%value(dphi0$)
phi0_err_saved = ele%value(phi0_err$) 
ele%value(phi0$) = 0
ele%value(dphi0$) = 0
ele%value(phi0_err$) = 0

tracking_method_saved = ele%tracking_method
if (ele%tracking_method == bmad_standard$) ele%tracking_method = runge_kutta$

theta_max = phi0_ref   ! Init guess
if (ele%key == rfcavity$) theta_max = theta_max - 0.25

theta_max_old = 100 ! Number far from unity
dtheta = 0.05
theta_tol = 1d-5
pz_tol = 1d-7

! See if %phi0_ref and %field_scale are already set correctly

pz_plus  = -neg_pz_calc(theta_max + 2 * theta_tol)
pz_minus = -neg_pz_calc(theta_max - 2 * theta_tol)
pz_max = -neg_pz_calc(theta_max)

call convert_pc_to ((1 + pz_max) * ele%value(p0c$), param%particle, e_tot = e_tot)
f_correct = wanted_de / (e_tot - e_tot_start)

if (pz_max > pz_plus .and. pz_max > pz_minus .and. abs(f_correct - 1) < 2 * pz_tol) return

! Now adjust %field_scale for the correct acceleration at the phase for maximum accelleration. 

main_loop: do

  ! Find approximately the phase for maximum acceleration.
  ! First go in +theta direction until pz decreases.

  step_up_seen = .false.
  do i = 1, 10
    theta = theta_max + dtheta
    pz = -neg_pz_calc(theta)
    if (pz < pz_max) exit
    pz_max = pz
    theta_max = theta
    step_up_seen = .true.
    if (i == 10) then  ! field too strong and always loosing particles
      field_scale = field_scale / 10
      pz = -neg_pz_calc(theta)
      cycle main_loop
    endif
  enddo

  pz_plus = pz

  ! If needed: Now go in -theta direction until pz decreases

  if (.not. step_up_seen) then
    do
      theta = theta_max - dtheta
      pz = -neg_pz_calc(theta)
      if (pz < pz_max) exit
      pz_max = pz
      theta_max = theta
    enddo
  endif

  pz_minus = pz

  ! Quadradic interpolation to get the maximum phase.
  ! Formula: pz = a + b*dt + c*dt^2 where dt = (theta-theta_max) / dtheta

  b = (pz_plus - pz_minus) / 2
  c = pz_plus - pz_max - b

  theta_max = theta_max - b * dtheta / (2 * c)
  pz_max = -neg_pz_calc(theta_max)

  ! Now scale %field_scale
  ! f_correct = dE(design) / dE (from tracking)
  ! Can overshoot so if f_correct is too large then scale back by a factor of 10

  call convert_pc_to ((1 + pz_max) * ele%value(p0c$), param%particle, e_tot = e_tot)
  f_correct = wanted_de / (e_tot - e_tot_start)
  if (f_correct > 1000) f_correct = max(1000.0_rp, f_correct / 10)
  field_scale = field_scale * f_correct

  if (abs(f_correct - 1) < pz_tol .and. abs(theta_max-theta_max_old) < theta_tol) exit
  theta_max_old = theta_max

  dtheta = 0.05
  if (abs(f_correct - 1) < 0.1) dtheta = max(theta_tol, 0.1*sqrt(2*abs(f_correct - 1))/twopi)

  pz_max = -neg_pz_calc(theta_max)

enddo main_loop

! For an rfcavity now find the zero crossing with negative slope which is
! about 90deg away from max acceleration.

if (ele%key == rfcavity$) then
  ele%value(phi0_max$) = phi0_ref  ! Save for use with OPAL
  dtheta = 0.1
  do
    theta = theta_max + dtheta
    pz = -neg_pz_calc(theta)
    if (pz < 0) exit
    theta_max = theta
  enddo
  phi0_ref = modulo2 (zbrent(neg_pz_calc, theta_max, theta_max+dtheta, 1d-9), 0.5_rp)
endif

! Cleanup

ele%value(phi0$) = phi0_saved
ele%value(dphi0$) = dphi0_saved
ele%value(phi0_err$) = phi0_err_saved
ele%tracking_method = tracking_method_saved

end subroutine rf_accel_mode_adjust_phase_and_amp

!----------------------------------------------------------------

function neg_pz_calc (theta) result (neg_pz)

implicit none

type (coord_struct) start_orb, end_orb
real(rp), intent(in) :: theta
real(rp) neg_pz

! brent finds minima so need to flip the final energy

phi0_ref = theta
call track1 (start_orb, ele_com, param_com, end_orb)
neg_pz = -end_orb%vec(6)
if (param_com%lost) neg_pz = 1

n_loop = n_loop + 1

end function

end module
