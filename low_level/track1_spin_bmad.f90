!+
! subroutine track1_spin_bmad (start_orb, ele, param, end_orb)
!
! Bmad_standard particle spin tracking through a single element.
!
! Note: spin tracking through a patch element is handled in track_a_patch since
! this is needed by runge_kutta tracking.
!
! Input :
!   start_orb  -- Coord_struct: Starting coords.
!   ele        -- Ele_struct: Element to track through.
!   param      -- lat_param_struct: Beam parameters.
!   end_orb    -- Coord_struct: Ending coords.
!
! Output:
!   end_orb    -- Coord_struct:
!     %spin(3)       -- Ending spin
!-

subroutine track1_spin_bmad (start_orb, ele, param, end_orb)

use bmad_routine_interface, dummy => track1_spin_bmad

implicit none

type (coord_struct) :: start_orb, end_orb
type (coord_struct) :: temp_start, temp_end, ave_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (fringe_field_info_struct) fringe_info
type (em_field_struct) field

real(rp) spline_x(0:3), spline_y(0:3), omega(3), s_edge_track, s_end_lab
real(rp) voltage, k_rf, phase
integer key, dir

character(*), parameter :: r_name = 'track1_spin_bmad'

! Spin tracking handled by track_a_patch for patch elements.

dir = start_orb%direction * start_orb%time_dir

if (ele%key == patch$) return

! crab_cavity

if (ele%key == crab_cavity$) then
  ave_orb = start_orb
  ave_orb%vec  = (start_orb%vec  + end_orb%vec) / 2
  ave_orb%t    = (start_orb%t    + end_orb%t) / 2
  ave_orb%beta = (start_orb%beta + end_orb%beta) / 2

  voltage = e_accel_field(ele, voltage$) * rel_tracking_charge_to_mass(ave_orb, param%particle)
  k_rf = twopi * ele%value(rf_frequency$) / c_light
  phase = twopi * (ele%value(phi0$) + ele%value(phi0_multipass$) - &
          (particle_rf_time (ave_orb, ele, .false.) - rf_ref_time_offset(ele)) * ele%value(rf_frequency$))

  field%B(2) = -voltage * sin(phase) / c_light
  field%E(3) = voltage * k_rf * ave_orb%beta * ave_orb%vec(1) * cos(phase)

  omega = spin_omega(field, ave_orb, start_orb%direction * ele%orientation) * start_orb%time_dir
  call rotate_spin(omega, end_orb%spin)

  return
endif

! A slice_slave may or may not span a fringe. calc_next_fringe_edge will figure this out.

if (dir == 1) then
  s_end_lab = ele%value(l$)
else
  s_end_lab = 0
endif

temp_start = start_orb
call calc_next_fringe_edge (ele, s_edge_track, fringe_info, temp_start, .true.)

call offset_particle (ele, set$, temp_start, set_hvkicks = .false., set_spin = .true.)

if (fringe_info%has_fringe .and. fringe_info%particle_at == first_track_edge$) then
  if (fringe_info%ds_edge /= 0) call track_a_drift (temp_start, fringe_info%ds_edge)
  call apply_element_edge_kick (temp_start, fringe_info, ele, param, .true.)
  call calc_next_fringe_edge (ele, s_edge_track, fringe_info, end_orb, .false.)
endif

temp_end  = end_orb

call offset_particle (ele, set$, temp_end, set_hvkicks = .false., s_pos = s_end_lab)

if (fringe_info%has_fringe .and. fringe_info%particle_at == second_track_edge$) then
  if (fringe_info%ds_edge /= 0) call track_a_drift (temp_end, fringe_info%ds_edge)
  temp_end%species = antiparticle(temp_end%species)  ! To reverse element edge kick
  call apply_element_edge_kick (temp_end, fringe_info, ele, param, .true.)
  temp_end%species = end_orb%species
endif

temp_end%spin = temp_start%spin

! temp_start%s = temp_end%s when, for example, in an RF cavity when ele is a super_slave covering
! a region outside the active region.

if (ele%value(l$) == 0 .or. ele%key == multipole$ .or. ele%key == ab_multipole$ .or. ele%key == taylor$) then
  temp_end%vec = (temp_end%vec + temp_start%vec) / 2
  call multipole_spin_tracking (ele, param, temp_end)
elseif (temp_start%s /= temp_end%s) then
  call spline_fit_orbit (ele, temp_start, temp_end, spline_x, spline_y)
  omega = trapzd_omega (ele, dir, spline_x, spline_y, temp_start, temp_end, param)
  if (ele%key == sbend$) omega = omega + [0.0_rp, ele%value(g$)*ele%value(l$)*dir*ele%orientation, 0.0_rp]
  call rotate_spin(omega, temp_end%spin)
endif

!----------

if (fringe_info%has_fringe .and. fringe_info%particle_at == second_track_edge$) then
  call apply_element_edge_kick (temp_end, fringe_info, ele, param, .true.)
endif

call offset_particle (ele, unset$, temp_end, set_hvkicks = .false., set_spin = .true.)

end_orb%spin = temp_end%spin

!---------------------------------------------------------------------------------------------
contains

function trapzd_omega (ele, dir, spline_x, spline_y, start_orb, end_orb, param) result (omega)

use super_recipes_mod, only: super_polint

implicit none

type q_array_struct
  real(rp) h
  real(rp) omega(3)
end type

integer, parameter ::  j_max = 10

type (q_array_struct) q_array(j_max), z(0:512)
type (ele_struct) ele
type (coord_struct) start_orb, end_orb, orb
type (lat_param_struct) param

real(rp) s0, s1, del_s, s, spline_x(0:3), spline_y(0:3), omega(3)
real(rp) dint, eps, quat(0:3)
real(rp), parameter :: eps_rel = 1d-5, eps_abs = 1d-8

integer dir
integer j, k, n, n_pts

!

s0 = (0             + bmad_com%significant_length/10)
s1 = (ele%value(l$) - bmad_com%significant_length/10)

q_array(1)%h = 1
z(0)%omega = omega_func(s0, dir, spline_x, spline_y, start_orb, end_orb, ele, param)
z(1)%omega = omega_func(s1, dir, spline_x, spline_y, start_orb, end_orb, ele, param)

del_s = s1 - s0
q_array(1)%omega = quat_to_omega(quat_mul(omega_to_quat(z(1)%omega * abs(del_s) / 2), omega_to_quat(z(0)%omega * abs(del_s) / 2)))

do j = 2, j_max
  ! This is trapzd from NR
  n_pts = 2**(j-2)
  del_s = (s1 - s0) / (2 * n_pts)
  quat = omega_to_quat(z(0)%omega * abs(del_s) / 2)

  z(2:2*n_pts:2) = z(1:n_pts)

  do n = 1, n_pts
    s = s0 + del_s * (2*n - 1)
    z(2*n-1)%omega = omega_func(s, dir, spline_x, spline_y, start_orb, end_orb, ele, param)
    quat = quat_mul(omega_to_quat(z(2*n-1)%omega * abs(del_s)), quat)
    if (n == n_pts) del_s = del_s / 2
    quat = quat_mul(omega_to_quat(z(2*n)%omega * abs(del_s)), quat)
  enddo
  
  q_array(j)%omega = quat_to_omega(quat)
  q_array(j)%h = q_array(j-1)%h / 4

  eps = eps_abs + eps_rel * sum(abs(q_array(j)%omega))

  select case (ele%key)
  case (wiggler$, rfcavity$, lcavity$);  if (j < 5) cycle  ! Cannot trust until have enough points
  end select

  do k = 1, 3
    call super_polint (q_array(1:j)%h, q_array(1:j)%omega(k), 0.0_rp, omega(k), dint)
    if (abs(dint) > eps .and. j < j_max) exit ! Failed test. Note: Last loop with j = j_max -> no test.
    if (k == 3) return                        ! Passed all tests or last loop
  enddo

enddo

end function trapzd_omega

!-----------------------------------------------------------------------------------
! contains

function omega_func (s_eval, dir, spline_x, spline_y, start_orb, end_orb, ele, param) result (omega)

implicit none

type (coord_struct) start_orb, end_orb, orb
type (ele_struct) ele
type (em_field_struct) field
type (lat_param_struct) param

real(rp) s_eval, spline_x(0:3), spline_y(0:3), omega(3), B(3)
real(rp) ds, dss, s_tot, s2

integer dir

!

ds = s_eval
s_tot = abs(end_orb%s - start_orb%s)

orb = end_orb

orb%vec(5) = start_orb%vec(5) * (s_tot - ds) / s_tot + end_orb%vec(5) * ds / s_tot
orb%vec(6) = start_orb%vec(6) * (s_tot - ds) / s_tot + end_orb%vec(6) * ds / s_tot

dss = ds * start_orb%time_dir
orb%vec(1) =                     spline_x(0) + spline_x(1) * dss + spline_x(2) * dss**2 + spline_x(3) * dss**3
orb%vec(2) = (1 + orb%vec(6)) * (spline_x(1) + 2 * spline_x(2) * dss + 3 * spline_x(3) * dss**2)
orb%vec(3) =                     spline_y(0) + spline_y(1) * dss + spline_y(2) * dss**2 + spline_y(3) * dss**3
orb%vec(4) = (1 + orb%vec(6)) * (spline_y(1) + 2 * spline_y(2) * dss + 3 * spline_y(3) * dss**2)

orb%t      = start_orb%t      * (s_tot - ds) / s_tot + end_orb%t      * ds / s_tot
orb%beta   = start_orb%beta   * (s_tot - ds) / s_tot + end_orb%beta   * ds / s_tot

if (dir == 1) then
  s2 = s_eval
else
  s2 = ele%value(l$) - s_eval
endif

call em_field_calc (ele, param, s2, orb, .true., field)

! 1 + g*x term comes from the curved coordinates.

omega = spin_omega (field, orb, start_orb%direction * ele%orientation) * start_orb%time_dir
if (ele%key == sbend$) omega = (1 + ele%value(g$) * orb%vec(1)) * omega

end function omega_func

end subroutine track1_spin_bmad

