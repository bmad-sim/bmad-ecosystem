module sad_mod

use track1_mod
use make_mat6_mod

contains

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine sad_mult_track_and_mat (ele, param, start_orb, end_orb, mat6)
!
! Routine to track a particle through a sad_mult element.
!
! Module needed:
!   use sad_mod
!
! Input:
!   ele          -- Ele_struct: Sad_mult element.
!   param        -- lat_param_struct: Lattice parameters.
!   start_orb    -- Coord_struct: Starting position.
!
! Output:
!   end_orb      -- Coord_struct: End position.
!   mat6(6,6)    -- real(rp), optional: Transfer matrix. 
!-

subroutine sad_mult_track_and_mat (ele, param, start_orb, end_orb, mat6)

implicit none

type (coord_struct) :: start_orb, end_orb
type (ele_struct), target :: ele, ele2, ele3
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) rel_pc, dz4_coef(4,4), mass, e_tot
real(rp) ks, k1, length, z_start, charge_dir, kx, ky
real(rp) xp_start, yp_start, mat4(4,4), mat1(6,6), f1, f2, ll, k0
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp) :: vec0(6), kmat(6,6)

integer n, nd, orientation, n_div, np_max, physical_end, fringe_at, ix_pole_max

logical fringe_here

character(*), parameter :: r_name = 'sad_mult_track_and_mat'

!

if (ele%value(rf_frequency$) /= 0) then
  call out_io (s_fatal$, r_name, 'RF CAVITY NOT YET IMPLEMENTED FOR SAD_MULT ELEMENTS!')
  if (global_com%exit_on_error) call err_exit
  return
endif

!

end_orb = start_orb 
length = ele%value(l$)
rel_pc = 1 + end_orb%vec(6)
n_div = nint(ele%value(num_steps$))

rel_pc = 1 + end_orb%vec(6)
orientation = ele%orientation * end_orb%direction
charge_dir = rel_tracking_charge_to_mass(end_orb, param) * orientation

if (present(mat6)) call mat_make_unit(mat6)

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
  call offset_particle (ele2, param, set$, end_orb, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

  if (ix_pole_max > -1) then
    call multipole_kicks (knl, tilt, param%particle, ele, end_orb)
    if (present(mat6)) then
      call multipole_kick_mat (knl, tilt, param%particle, ele, end_orb, 1.0_rp, mat6)
    endif
  endif

  call offset_particle (ele2, param, unset$, end_orb, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

  if (present(mat6)) then
    ele%vec0 = end_orb%vec - matmul(mat6, start_orb%vec)
  endif

  end_orb%s = ele%s
  end_orb%location = downstream_end$
  return
endif

! Go to frame of reference of the multipole quad component

ks = rel_tracking_charge_to_mass(end_orb, param) * ele%value(ks$)
k1 = charge_dir * knl(1) / length

if (ele%value(x_pitch_mult$) /= 0 .or. ele%value(y_pitch_mult$) /= 0) then
  kx = knl(0) * cos(tilt(0)) - ks * ele%value(x_pitch_mult$)
  ky = knl(0) * sin(tilt(0)) + ks * ele%value(y_pitch_mult$)
  knl(0) = norm2([kx, ky])
  tilt(0) = atan2(ky, kx)
endif

if (ele%value(x_offset_mult$) /= 0 .or. ele%value(y_offset_mult$) /= 0) then
  end_orb%vec(2) = end_orb%vec(2) + ele%value(y_offset_mult$) * ks / 2
  end_orb%vec(4) = end_orb%vec(4) - ele%value(x_offset_mult$) * ks / 2
endif

ele2%value(tilt_tot$) = tilt(1) 
tilt = tilt - tilt(1)

call multipole_kt_to_ab (knl, tilt, a_pole, b_pole)
knl(1) = 0 ! So multipole_kicks does not conflict with sol_quad calc. 

call offset_particle (ele2, param, set$, end_orb, set_multipoles = .false., set_hvkicks = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

! Entrance edge kicks
! The multipole hard edge routine takes care of the quadrupole hard edge.

k0 = knl(0)/length

call hard_multipole_edge_kick (ele, param, first_track_edge$, end_orb, mat6, present(mat6), a_pole, b_pole)
if (orbit_too_large (end_orb, param)) return
call soft_quadrupole_edge_kick (ele, param, first_track_edge$, end_orb, mat6, present(mat6))
call no_edge_angle_hard_bend_edge_kick (ele, param, first_track_edge$, end_orb, mat6, present(mat6), k0, tilt(0))
if (end_orb%state /= alive$) return
call soft_bend_edge_kick (ele, param, first_track_edge$, end_orb, mat6, present(mat6), k0, tilt(0))

! Body

knl = knl / n_div
call transfer_ele(ele, ele3)

do nd = 0, n_div

  ll = length / n_div
  if (nd == 0 .or. nd == n_div) ll = ll / 2
  ele3%value(l$) = ll

  ! 

  if (abs(k1) < 1d-40) then
    call solenoid_track_and_mat (ele3, param, end_orb, end_orb, mat1)
    if (present(mat6)) mat6 = matmul(mat1, mat6)

  else
    if (present(mat6)) then
      call sol_quad_mat6_calc (ks, k1, ll, end_orb%vec, mat1)
      if (present(mat6)) mat6 = matmul(mat1, mat6)
    endif
    vec0 = 0
    vec0(6) = end_orb%vec(6)
    call sol_quad_mat6_calc (ks, k1, ll, vec0, mat1, dz4_coef)
    end_orb%vec(5) = end_orb%vec(5) + sum(end_orb%vec(1:4) * matmul(dz4_coef, end_orb%vec(1:4))) 
    end_orb%vec(1:4) = matmul (mat1(1:4,1:4), end_orb%vec(1:4))
  endif

  ! multipole kicks

  if (nd == n_div) exit

  call multipole_kicks (knl, tilt, param%particle, ele, end_orb)

  if (present(mat6)) then
    call multipole_kick_mat (knl, tilt, param%particle, ele, end_orb, 1.0_rp, mat1)
    mat6(2,:) = mat6(2,:) + mat1(2,1) * mat6(1,:) + mat1(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat1(4,1) * mat6(1,:) + mat1(4,3) * mat6(3,:)
  endif

  ! Check for end_orb too large to prevent infinities.

  if (orbit_too_large (end_orb)) return
enddo



! Exit edge kicks

call soft_bend_edge_kick (ele, param, second_track_edge$, end_orb, mat6, present(mat6), k0, tilt(0))
call no_edge_angle_hard_bend_edge_kick (ele, param, second_track_edge$, end_orb, mat6, present(mat6), k0, tilt(0))
if (end_orb%state /= alive$) return
call soft_quadrupole_edge_kick (ele, param, second_track_edge$, end_orb, mat6, present(mat6))
call hard_multipole_edge_kick (ele, param, second_track_edge$, end_orb, mat6, present(mat6), a_pole, b_pole)
if (orbit_too_large (end_orb, param)) return

! End stuff

call offset_particle (ele2, param, unset$, end_orb, set_multipoles = .false., set_hvkicks = .false., &
                                                                                       mat6 = mat6, make_matrix = present(mat6))

if (ele%value(x_offset_mult$) /= 0 .or. ele%value(y_offset_mult$) /= 0) then
  ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) + ele%value(x_offset_mult$)
  ele2%value(y_offset_tot$) = ele%value(y_offset_tot$) + ele%value(y_offset_mult$)
  end_orb%vec(2) = end_orb%vec(2) - ele%value(y_offset_mult$) * ks / 2
  end_orb%vec(4) = end_orb%vec(4) + ele%value(x_offset_mult$) * ks / 2
endif

end_orb%vec(5) = end_orb%vec(5) + low_energy_z_correction (end_orb, ele2, ele2%value(l$), mat6, present(mat6))

!

end_orb%t = start_orb%t + (length + start_orb%vec(5) - end_orb%vec(5)) / (end_orb%beta * c_light)
end_orb%s = ele%s
end_orb%location = downstream_end$

end subroutine sad_mult_track_and_mat 

end module
