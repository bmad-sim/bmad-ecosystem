!+
! Subroutine track_a_sad_mult (orbit, ele, param, mat6)
!
! Routine to track a particle through a sad_mult element.
!
! Input:
!   ele          -- Ele_struct: Sad_mult element.
!   param        -- lat_param_struct: Lattice parameters.
!   orbit        -- Coord_struct: Starting position.
!   mat6(6,6)    -- real(rp), optional: Transfer matrix up to the sad_mult.
!
! Output:
!   orbit        -- Coord_struct: End position.
!   mat6(6,6)    -- real(rp), optional: Transfer matrix. 
!-

subroutine track_a_sad_mult (orbit, ele, param, mat6)

use track1_mod, dummy1 => track_a_sad_mult
use make_mat6_mod, dummy => track_a_sad_mult

implicit none

type (coord_struct) :: orbit
type (ele_struct), target :: ele, ele2, ele3
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) rel_pc, dz4_coef(4,4), mass, e_tot
real(rp) ks, k1, length, z_start, t_start, charge_dir, kx, ky
real(rp) xp_start, yp_start, mat4(4,4), mat1(6,6), f1, f2, ll, k0
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) :: vec0(6), kmat(6,6)

integer n, nd, orientation, n_div, np_max, physical_end, fringe_at, ix_pole_max

logical fringe_here

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
orientation = ele%orientation * orbit%direction
charge_dir = rel_tracking_charge_to_mass(orbit, param) * orientation

knl = 0; tilt = 0
call multipole_ele_to_kt (ele, .true., ix_pole_max, knl, tilt)

! Setup ele2 which is used in offset_particle

call transfer_ele(ele, ele2)
ele2%value(x_pitch_tot$) = ele%value(x_pitch_tot$) + ele%value(x_pitch_mult$)
ele2%value(y_pitch_tot$) = ele%value(y_pitch_tot$) + ele%value(y_pitch_mult$)
ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) + ele%value(x_offset_mult$)
ele2%value(y_offset_tot$) = ele%value(y_offset_tot$) + ele%value(y_offset_mult$)

! If element has zero length then the SAD ignores f1 and f2.

if (length == 0) then
  call offset_particle (ele2, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

  if (ix_pole_max > -1) then
    call multipole_kicks (knl, tilt, param%particle, ele, orbit)
    if (present(mat6)) then
      call multipole_kick_mat (knl, tilt, param%particle, ele, orbit, 1.0_rp, mat6)
    endif
  endif

  call offset_particle (ele2, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

  orbit%s = ele%s
  orbit%location = downstream_end$
  return
endif

! Go to frame of reference of the multipole quad component

ks = rel_tracking_charge_to_mass(orbit, param) * ele%value(ks$)
k1 = charge_dir * knl(1) / length

if (ele%value(x_pitch_mult$) /= 0 .or. ele%value(y_pitch_mult$) /= 0) then
  kx = knl(0) * cos(tilt(0)) - ks * ele%value(x_pitch_mult$)
  ky = knl(0) * sin(tilt(0)) + ks * ele%value(y_pitch_mult$)
  knl(0) = norm2([kx, ky])
  tilt(0) = atan2(ky, kx)
endif

if (ele%value(x_offset_mult$) /= 0 .or. ele%value(y_offset_mult$) /= 0) then
  orbit%vec(2) = orbit%vec(2) + ele%value(y_offset_mult$) * ks / 2
  orbit%vec(4) = orbit%vec(4) - ele%value(x_offset_mult$) * ks / 2
endif

ele2%value(tilt_tot$) = tilt(1) 
tilt = tilt - tilt(1)

call multipole_kt_to_ab (knl, tilt, a_pole, b_pole)
knl(1) = 0 ! So multipole_kicks does not conflict with sol_quad calc. 

call offset_particle (ele2, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

! Entrance edge kicks
! The multipole hard edge routine takes care of the quadrupole hard edge.

k0 = knl(0)/length

call hard_multipole_edge_kick (ele, param, first_track_edge$, orbit, mat6, present(mat6), a_pole, b_pole)
if (orbit_too_large (orbit, param)) return
call soft_quadrupole_edge_kick (ele, param, first_track_edge$, orbit, mat6, present(mat6))
call no_edge_angle_hard_bend_edge_kick (ele, param, first_track_edge$, orbit, mat6, present(mat6), k0, tilt(0))
if (orbit%state /= alive$) return
call soft_bend_edge_kick (ele, param, first_track_edge$, orbit, mat6, present(mat6), k0, tilt(0))

! Body

knl = knl / n_div
call transfer_ele(ele, ele3)

do nd = 0, n_div

  ll = length / n_div
  if (nd == 0 .or. nd == n_div) ll = ll / 2
  ele3%value(l$) = ll

  ! 

  if (abs(k1) < 1d-40) then
    call solenoid_track_and_mat (ele3, param, orbit, orbit, mat1)
    if (present(mat6)) mat6 = matmul(mat1, mat6)

  else
    if (present(mat6)) then
      call sol_quad_mat6_calc (ks, k1, ll, orbit%vec, mat1)
      mat6 = matmul(mat1, mat6)
    endif
    vec0 = 0
    vec0(6) = orbit%vec(6)
    call sol_quad_mat6_calc (ks, k1, ll, vec0, mat1, dz4_coef)
    orbit%vec(5) = orbit%vec(5) + sum(orbit%vec(1:4) * matmul(dz4_coef, orbit%vec(1:4))) 
    orbit%vec(1:4) = matmul (mat1(1:4,1:4), orbit%vec(1:4))
  endif

  ! multipole kicks

  if (nd == n_div) exit

  call multipole_kicks (knl, tilt, param%particle, ele, orbit)

  if (present(mat6)) then
    call multipole_kick_mat (knl, tilt, param%particle, ele, orbit, 1.0_rp, mat1)
    mat6(2,:) = mat6(2,:) + mat1(2,1) * mat6(1,:) + mat1(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat1(4,1) * mat6(1,:) + mat1(4,3) * mat6(3,:)
  endif

  ! Check for orbit too large to prevent infinities.

  if (orbit_too_large (orbit)) return
enddo



! Exit edge kicks

call soft_bend_edge_kick (ele, param, second_track_edge$, orbit, mat6, present(mat6), k0, tilt(0))
call no_edge_angle_hard_bend_edge_kick (ele, param, second_track_edge$, orbit, mat6, present(mat6), k0, tilt(0))
if (orbit%state /= alive$) return
call soft_quadrupole_edge_kick (ele, param, second_track_edge$, orbit, mat6, present(mat6))
call hard_multipole_edge_kick (ele, param, second_track_edge$, orbit, mat6, present(mat6), a_pole, b_pole)
if (orbit_too_large (orbit, param)) return

! End stuff

call offset_particle (ele2, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

if (ele%value(x_offset_mult$) /= 0 .or. ele%value(y_offset_mult$) /= 0) then
  ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) + ele%value(x_offset_mult$)
  ele2%value(y_offset_tot$) = ele%value(y_offset_tot$) + ele%value(y_offset_mult$)
  orbit%vec(2) = orbit%vec(2) - ele%value(y_offset_mult$) * ks / 2
  orbit%vec(4) = orbit%vec(4) + ele%value(x_offset_mult$) * ks / 2
endif

orbit%vec(5) = orbit%vec(5) + low_energy_z_correction (orbit, ele2, ele2%value(l$), mat6, present(mat6))
orbit%t = t_start + length * ele%value(E_tot$) / (c_light * ele%value(p0c$)) - (orbit%vec(5) - z_start) / (orbit%beta * c_light)

!

if (orbit%direction == 1) then
  orbit%s = ele%s
else
  orbit%s = ele%s_start
endif
orbit%location = downstream_end$

end subroutine track_a_sad_mult 
