!+
! Subroutine track_a_sad_mult (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking of a particle through a sad_mult element.
!
! Input:
!   ele          -- Ele_struct: Sad_mult element.
!   param        -- lat_param_struct: Lattice parameters.
!   orbit        -- coord_struct: Starting position.
!   mat6(6,6)    -- real(rp), optional: Transfer matrix up to the sad_mult.
!   make_matrix  -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit        -- coord_struct: End position.
!   mat6(6,6)    -- real(rp), optional: Transfer matrix. 
!-

subroutine track_a_sad_mult (orbit, ele, param, mat6, make_matrix)

use fringe_mod, dummy1 => track_a_sad_mult

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) rel_pc, dz4_coef(4,4), mass, e_tot
real(rp) ks, k1, tilt1, length, z_start, t_start, charge_dir, kx, ky, sol_center(2)
real(rp) xp_start, yp_start, mat4(4,4), mat1(6,6), f1, f2, ll, r_scale
real(rp) a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx), a2_pole(0:n_pole_maxx), b2_pole(0:n_pole_maxx)
real(rp) :: vec0(6), kmat(6,6)

integer n, nd, orientation, n_div, np_max, physical_end, fringe_at, ix_pole_max
logical, optional :: make_matrix

character(*), parameter :: r_name = 'track_a_sad_mult'

!

if (ele%value(rf_frequency$) /= 0) then
  call out_io (s_fatal$, r_name, 'RF CAVITY NOT YET IMPLEMENTED FOR SAD_MULT ELEMENTS!')
  if (global_com%exit_on_error) call err_exit
  return
endif

!

z_start = orbit%vec(5)
t_start = orbit%t
length = ele%value(l$)
rel_pc = 1 + orbit%vec(6)
n_div = nint(ele%value(num_steps$))

rel_pc = 1 + orbit%vec(6)
orientation = ele%orientation * orbit%direction * orbit%time_dir
charge_dir = rel_tracking_charge_to_mass(orbit, param%particle) * orientation

call multipole_ele_to_ab (ele, .false., ix_pole_max, a_pole, b_pole)

! If element has zero length then the SAD ignores f1 and f2.

if (length == 0) then
  call offset_particle (ele, set$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)
  call ab_multipole_kicks (a_pole, b_pole, ix_pole_max, ele, orbit, magnetic$, 1.0_rp, mat6, make_matrix)
  call offset_particle (ele, unset$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

  orbit%s = ele%s
  orbit%location = downstream_end$
  return
endif

! 

call multipole1_ab_to_kt(a_pole(1), b_pole(1), 1, k1, tilt1)
k1 = charge_dir * k1 / length
ks = rel_tracking_charge_to_mass(orbit, param%particle) * ele%value(ks$)
sol_center = rot_2d ([ele%value(x_offset_mult$), ele%value(y_offset_mult$)], -ele%value(tilt$))

call offset_particle (ele, set$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

orbit%vec(2) = orbit%vec(2) + 0.5_rp * sol_center(2) * ks
orbit%vec(4) = orbit%vec(4) - 0.5_rp * sol_center(1) * ks

! Entrance edge kicks
! The multipole hard edge routine takes care of the quadrupole hard edge.

call hard_multipole_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix)
if (orbit_too_large (orbit, param)) return
call soft_quadrupole_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix)
call sad_mult_hard_bend_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix)
if (orbit%state /= alive$) return
call sad_soft_bend_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix)

! Body

r_scale = 1.0_rp / n_div
a2_pole = a_pole;  a2_pole(1) = 0  ! Quad term taken care of by sol_quad_mat6_calc.
b2_pole = b_pole;  b2_pole(1) = 0

do nd = 0, n_div

  ll = length / n_div
  if (nd == 0 .or. nd == n_div) ll = ll / 2

  ! 

  if (abs(k1) < 1d-40) then
    call solenoid_track_and_mat (ele, ll, param, orbit, orbit, mat6, make_matrix)
  else
    call sol_quad_mat6_calc (ks, k1, tilt1, ll, ele, orbit, mat6, make_matrix)
  endif

  ! multipole kicks

  if (nd == n_div) exit
  call ab_multipole_kicks (a2_pole, b2_pole, ix_pole_max, ele, orbit, magnetic$, r_scale, mat6, make_matrix)
  if (orbit_too_large (orbit)) return
enddo



! Exit edge kicks

call sad_soft_bend_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix)
call sad_mult_hard_bend_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix)
if (orbit%state /= alive$) return
call soft_quadrupole_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix)
call hard_multipole_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix)
if (orbit_too_large (orbit, param)) return

! End stuff

orbit%vec(2) = orbit%vec(2) - 0.5_rp * sol_center(2) * ks
orbit%vec(4) = orbit%vec(4) + 0.5_rp * sol_center(1) * ks

call offset_particle (ele, unset$, orbit, set_hvkicks = .false., mat6 = mat6, make_matrix = make_matrix)

orbit%t = t_start + length * ele%value(E_tot$) / (c_light * ele%value(p0c$)) - (orbit%vec(5) - z_start) / (orbit%beta * c_light)

!

if (orbit%direction*orbit%time_dir == 1) then
  orbit%s = ele%s
else
  orbit%s = ele%s_start
endif
orbit%location = downstream_end$

end subroutine track_a_sad_mult 
