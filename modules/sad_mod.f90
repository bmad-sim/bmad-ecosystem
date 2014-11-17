module sad_mod

use track1_mod
use make_mat6_mod

contains

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine sad_mult_track_and_mat (ele, param, start_orb, end_orb, end_in, make_matrix)
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
!   end_in       -- Logical: If True then end_orb will be taken as input. Not output as normal.
!   make_matrix  -- Logical: If True then make the transfer matrix.
!
! Output:
!   ele          -- Ele_struct: Sad_mult element.
!     %mat6(6,6)   -- Transfer matrix. 
!   end_orb      -- Coord_struct: End position.
!-

subroutine sad_mult_track_and_mat (ele, param, start_orb, end_orb, end_in, make_matrix)

implicit none

type (coord_struct) :: orbit, start_orb, end_orb
type (ele_struct), target :: ele, ele2
type (lat_param_struct) :: param

real(rp) rel_pc, dz4_coef(4,4), mass, e_tot
real(rp) ks, k1, length, z_start, charge_dir, kx, ky
real(rp) xp_start, yp_start, mat4(4,4), mat1(6,6), f1, f2, ll, k0
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx), a_pole(0:n_pole_maxx), b_pole(0:n_pole_maxx)
real(rp), pointer :: mat6(:,:)
real(rp) :: vec0(6), kmat(6,6)

integer n, nd, orientation, n_div, np_max, physical_end, fringe_at

logical make_matrix, end_in, has_nonzero_pole, fringe_here

character(*), parameter :: r_name = 'sad_mult_track_and_mat'

!

if (ele%value(rf_frequency$) /= 0) then
  call out_io (s_fatal$, r_name, 'RF CAVITY NOT YET IMPLEMENTED FOR SAD_MULT ELEMENTS!')
  if (global_com%exit_on_error) call err_exit
  return
endif

!

orbit = start_orb 
length = ele%value(l$)
rel_pc = 1 + orbit%vec(6)
n_div = nint(ele%value(num_steps$))

rel_pc = 1 + orbit%vec(6)
orientation = ele%orientation * orbit%direction
charge_dir = relative_tracking_charge(orbit, ele, param) * orientation
mat6 => ele%mat6

if (make_matrix) call mat_make_unit(mat6)

knl = 0; tilt = 0
call multipole_ele_to_kt (ele, .true., has_nonzero_pole, knl, tilt)
knl = knl * charge_dir

! Setup ele2 which is used in offset_particle

call transfer_ele(ele, ele2)
ele2%value(x_pitch_tot$) = ele%value(x_pitch_tot$) + ele%value(x_pitch_mult$)
ele2%value(y_pitch_tot$) = ele%value(y_pitch_tot$) + ele%value(y_pitch_mult$)
ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) + ele%value(x_offset_mult$)
ele2%value(y_offset_tot$) = ele%value(y_offset_tot$) + ele%value(y_offset_mult$)

! If element has zero length then the SAD ignores f1 and f2.

if (length == 0) then
  call offset_particle (ele2, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false.)

  if (has_nonzero_pole) then
    call multipole_kicks (knl, tilt, orbit)
    if (make_matrix) then
      call multipole_kick_mat (knl, tilt, orbit, 1.0_rp, mat6)
    endif
  endif

  call offset_particle (ele2, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false.)

  if (make_matrix) then
    call mat6_add_pitch (ele2%value(x_pitch_tot$), ele2%value(y_pitch_tot$), ele2%orientation, mat6)
    ele%vec0 = orbit%vec - matmul(mat6, start_orb%vec)
  endif

  orbit%s = ele%s
  if (.not. end_in) then
    end_orb = orbit
    end_orb%location = downstream_end$
  endif
  return
endif

! Go to frame of reference of the multipole quad component

ks = relative_tracking_charge(orbit, ele2, param) * ele%value(ks$)
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

call offset_particle (ele2, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false.)

! Entrance edge kicks
! The multipole hard edge routine takes care of the quadrupole hard edge.

k0 = knl(0)/length

call hard_multipole_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix, a_pole, b_pole)
call soft_quadrupole_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix)
call no_edge_angle_hard_bend_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix, k0, tilt(0))
call soft_bend_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix, k0, tilt(0))

! Body

knl = knl / n_div

do nd = 0, n_div

  ll = length / n_div
  if (nd == 0 .or. nd == n_div) ll = ll / 2

  ! Matrix step

  if (make_matrix) then
    if (abs(k1) < 1d-40) then
      call solenoid_mat6_calc (ks, ll, 0.0_rp, orbit, mat1)
    else
      call sol_quad_mat6_calc (ks, k1, ll, orbit%vec, mat1)
    endif
    mat6 = matmul(mat1, mat6)
  endif

  ! track step

  if (abs(k1) < 1d-40) then
    xp_start = orbit%vec(2) + ks * orbit%vec(3) / 2 
    yp_start = orbit%vec(4) - ks * orbit%vec(1) / 2
    call solenoid_mat4_calc (ks, ll, rel_pc, mat4)
    orbit%vec(5) = orbit%vec(5) - ll * (xp_start**2 + yp_start**2 ) / (2 * rel_pc**2)
    orbit%vec(1:4) = matmul (mat4, orbit%vec(1:4))
  else
    vec0 = 0
    vec0(6) = orbit%vec(6)
    call sol_quad_mat6_calc (ks, k1, ll, vec0, mat1, dz4_coef)
    orbit%vec(5) = orbit%vec(5) + sum(orbit%vec(1:4) * matmul(dz4_coef, orbit%vec(1:4))) 
    orbit%vec(1:4) = matmul (mat1(1:4,1:4), orbit%vec(1:4))
  endif

  ! multipole kicks

  if (nd == n_div) exit

  call multipole_kicks (knl, tilt, orbit)

  if (make_matrix) then
    call multipole_kick_mat (knl, tilt, orbit, 1.0_rp, mat1)
    mat6(2,:) = mat6(2,:) + mat1(2,1) * mat6(1,:) + mat1(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat1(4,1) * mat6(1,:) + mat1(4,3) * mat6(3,:)
  endif

  ! Check for orbit too large to prevent infinities.

  if (orbit_too_large (orbit)) then
    if (.not. end_in) end_orb = orbit
    return
  endif

enddo

! Exit edge kicks

call soft_bend_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix, k0, tilt(0))
call no_edge_angle_hard_bend_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix, k0, tilt(0))
call soft_quadrupole_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix)
call hard_multipole_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix, a_pole, b_pole)

! End stuff

call offset_particle (ele2, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false.)

if (ele%value(x_offset_mult$) /= 0 .or. ele%value(y_offset_mult$) /= 0) then
  ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) + ele%value(x_offset_mult$)
  ele2%value(y_offset_tot$) = ele%value(y_offset_tot$) + ele%value(y_offset_mult$)
  orbit%vec(2) = orbit%vec(2) - ele%value(y_offset_mult$) * ks / 2
  orbit%vec(4) = orbit%vec(4) + ele%value(x_offset_mult$) * ks / 2
endif

call track1_low_energy_z_correction (orbit, ele2, param)

if (make_matrix) then
  if (ele2%value(tilt_tot$) /= 0) call tilt_mat6 (mat6, ele2%value(tilt_tot$))

  call mat6_add_pitch (ele2%value(x_pitch_tot$), ele2%value(y_pitch_tot$), ele2%orientation, mat6)

  ! 1/gamma^2 m56 correction

  mass = mass_of(orbit%species)
  e_tot = ele%value(p0c$) * (1 + orbit%vec(6)) / orbit%beta
  mat6(5,6) = mat6(5,6) + length * mass**2 * ele%value(e_tot$) / e_tot**3

  ele%vec0 = orbit%vec - matmul(mat6, start_orb%vec)
endif

!

orbit%t = start_orb%t + (length + start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)
orbit%s = ele%s

if (.not. end_in) then
  end_orb = orbit
  end_orb%location = downstream_end$
endif

end subroutine sad_mult_track_and_mat 

end module
