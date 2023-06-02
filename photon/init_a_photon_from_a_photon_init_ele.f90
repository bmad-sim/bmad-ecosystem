!+
! Subroutine init_a_photon_from_a_photon_init_ele (ele, param, orbit)
!
! Routine to initialize a photon from an photon_init element.
! This routine is called by init_coord and is not meant to be called directly.
!
! Input:
!   ele           -- ele_struct: patch element.
!   param         -- lat_param_struct
!
! Output:
!   orbit         -- coord_struct: Coords after applying a patch transformation.
!-

Subroutine init_a_photon_from_a_photon_init_ele (ele, param, orbit)

use track1_photon_mod, except_dummy => init_a_photon_from_a_photon_init_ele
use super_recipes_mod, only: super_zbrent

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct) orbit
type (photon_element_struct), pointer :: ph
type (spline_struct) spl

real(rp) r(3), rv(2), rr, drr, p2, f, de
integer n, ix, status

character(*), parameter :: r_name = 'init_a_photon_from_a_photon_init_ele'

! Spatial position

if (nint(ele%value(spatial_distribution$)) == uniform$) then
  call ran_uniform(r)
  r = 2 * r - 1
elseif (nint(ele%value(spatial_distribution$)) == gaussian$) then
  call ran_gauss(r)
else
  call out_io (s_fatal$, r_name, 'BAD ELE SPATIAL_DISTRIBUTION: ' // &
                                             distribution_name(nint(ele%value(spatial_distribution$))))
  if (global_com%exit_on_error) call err_exit
  orbit%state = lost$
  return
endif

orbit%vec(1:5:2) = [ele%value(sig_x$) * r(1), ele%value(sig_y$) * r(2), ele%value(sig_z$) * r(3)]

! Direction of photon

select case (nint(ele%value(velocity_distribution$)))
case (spherical$)
  call point_photon_emission (ele, param, orbit, +1, twopi)  ! Includes targeting

case (uniform$)
  call ran_uniform(rv)
  rv = 2 * rv - 1
  orbit%vec(2:4:2) = [ele%value(x_pitch$), ele%value(y_pitch$)] + rv * [ele%value(sig_vx$), ele%value(sig_vy$)]
  orbit%vec(6) = sqrt(1 - orbit%vec(2)**2 - orbit%vec(4)**2)

case (gaussian$)
  call ran_gauss(rv)
  orbit%vec(2:4:2) = [ele%value(x_pitch$), ele%value(y_pitch$)] + rv * [ele%value(sig_vx$), ele%value(sig_vy$)]
  orbit%vec(6) = sqrt(1 - orbit%vec(2)**2 - orbit%vec(4)**2)
end select

! Energy of photon

select case (nint(ele%value(energy_distribution$)))
case (uniform$, gaussian$)
  if (nint(ele%value(energy_distribution$)) == uniform$) then
    call ran_uniform(rr)
    rr = 2 * rr - 1
  else
    call ran_gauss(rr)
  endif

  p2 = 1
  if (ele%value(E2_probability$) /= 0) call ran_uniform(p2)

  if (p2 < ele%value(E2_probability$)) then
    orbit%p0c = ele%value(sig_E2$) * rr + ele%value(E2_center$)
  else
    orbit%p0c = ele%value(sig_E$)  * rr + ele%value(E_center$)
  endif  

case (curve$)
  call ran_uniform(rr)
  if (.not. allocated(ele%photon%init_energy_prob)) then
    call out_io (s_fatal$, r_name, 'NO ENERGY_PROBABILITY_CURVE SPECIFIED FOR: ' // ele%name)
  endif
  ph => ele%photon
  n = ubound(ph%init_energy_prob, 1)
  f = ph%integrated_init_energy_prob(n)
  rr = rr * ph%integrated_init_energy_prob(n)
  ix = bracket_index(rr, ph%integrated_init_energy_prob, 1)
  drr = rr - ph%integrated_init_energy_prob(ix)
  spl = ph%init_energy_prob(ix)
  orbit%p0c = super_zbrent(e_curve_func, spl%x0, spl%x1, 0.0_rp, 1.0e-14_rp*ph%init_energy_prob(n)%x1, status)

case default
  call out_io (s_error$, r_name, 'BAD ELE%VALUE(ENERGY_DISTRIBUTION SETTING: ' // &
                                             distribution_name(nint(ele%value(energy_distribution$))))
  orbit%state = lost$
end select

!

if (is_true(ele%value(E_center_relative_to_ref$))) orbit%p0c = orbit%p0c + ele%value(p0c$) 

orbit%s = orbit%vec(5) + orbit%s + ele%value(z_offset_tot$)
orbit%t = 0

! Translate from element to lab coordinates
! and track to entrance end of ele

call offset_photon (ele, orbit, unset$)

call track_a_drift_photon (orbit, -orbit%s, .true.)

!---------------------------------------------------------------------
contains

function e_curve_func(energy, status) result (dprob)

real(rp), intent(in) :: energy
real(rp) dprob
integer status

!

dprob = spline1(spl, energy, -1) - drr

end function e_curve_func

end subroutine init_a_photon_from_a_photon_init_ele

