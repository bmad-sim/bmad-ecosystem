module sad_mod

use track1_mod
use make_mat6_mod

contains

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine track_a_sad_mult (orbit, ele, param)
!
! Routine to track a particle through a sad_mult element.
!
! Module needed:
!   use sad_mod
!
! Input:
!   orbit      -- Coord_struct: Starting position.
!   ele        -- Ele_struct: Sad_mult element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   orbit      -- Coord_struct: End position.
!-

subroutine track_a_sad_mult (orbit, ele, param)

implicit none

type (coord_struct) :: orbit, start_orb
type (ele_struct)  :: ele, ele2
type (lat_param_struct) :: param

real(rp) rel_pc, vec0(6), dz4_coef(4,4), mat6(6,6)
real(rp) ks, k1, length, z_start, charge_dir
real(rp) xp_start, yp_start, mat4(4,4)
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

integer orientation

logical has_nonzero

!

start_orb = orbit

length = ele%value(l$)
rel_pc = 1 + orbit%vec(6)
orientation = ele%orientation * start_orb%direction
charge_dir = param%rel_tracking_charge * orientation

call multipole_ele_to_kt (ele, param, .true., has_nonzero, knl, tilt)

ks = param%rel_tracking_charge * ele%value(ks$) / rel_pc
k1 = charge_dir * knl(1) / rel_pc / length

call transfer_ele(ele, ele2)
ele2%value(tilt_tot$) = tilt(1)
ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) - ele%value(x_offset_sol$)

call offset_particle (ele2, orbit, param, set$, set_multipoles = .false., set_hvkicks = .false.)

if (k1 == 0) then
  xp_start = orbit%vec(2) + ks * orbit%vec(3) / 2
  yp_start = orbit%vec(4) - ks * orbit%vec(1) / 2
  orbit%vec(5) = orbit%vec(5) - length * (xp_start**2 + yp_start**2 ) / 2

  call solenoid_mat4_calc (ks, length, mat4)
  orbit%vec(1:4) = matmul (mat4, orbit%vec(1:4))
else
  vec0 = 0
  call sol_quad_mat6_calc (ks, k1, length, mat6, vec0, dz4_coef)
  orbit%vec(5) = orbit%vec(5) + sum(orbit%vec(1:4) * matmul(dz4_coef, orbit%vec(1:4)))   
  orbit%vec(1:4) = matmul (mat6(1:4,1:4), orbit%vec(1:4))
endif

call offset_particle (ele2, orbit, param, unset$, set_multipoles = .false., set_hvkicks = .false.)

call track1_low_energy_z_correction (orbit, ele2, param)

orbit%t = start_orb%t + (ele%value(l$) + start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)
orbit%s = ele%s

end subroutine track_a_sad_mult

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine make_mat6_sad_mult (ele, param, c0, c1)
!
! Routine to make the transfer matrix through a sad_mult element.
!
! Module needed:
!   use sad_mod
!
! Input:
!   ele           -- ele_struct: sad_mult element.
!   param         -- lat_param_struct: Lattice parameters.
!   c0            -- coord_struct: Coordinates at beginning of element.
!   c1            -- coord_struct: Coordinates at end of element.
!
! Output:
!   ele           -- ele_struct: sad_mult element.
!     %mat6(6,6)
!     %vec0(6)
!-

subroutine make_mat6_sad_mult (ele, param, c0, c1)

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) param
type (coord_struct) :: c0, c1, c00

real(rp) length, rel_pc, charge_dir, ks, k1, e_tot, mass
real(rp), pointer :: mat6(:,:)
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

integer orientation

logical has_nonzero

!

length = ele%value(l$)
rel_pc = 1 + c0%vec(6)
orientation = ele%orientation * c0%direction
charge_dir = param%rel_tracking_charge * orientation
mat6 => ele%mat6

call multipole_ele_to_kt (ele, param, .true., has_nonzero, knl, tilt)

ks = param%rel_tracking_charge * ele%value(ks$)
k1 = charge_dir * knl(1) / length

ele%value(tilt_tot$) = tilt(1)
c00 = c0
call offset_particle (ele, c00, param, set$, set_multipoles = .false., set_hvkicks = .false.)

if (k1 == 0) then
  call solenoid_mat6_calc (ks, length, ele%value(tilt_tot$), c00, mat6)
else
  call sol_quad_mat6_calc (ks, k1, length, mat6, c00%vec)
endif

if (ele%value(tilt_tot$) /= 0) then
  call tilt_mat6 (mat6, ele%value(tilt_tot$))
endif

call mat6_add_pitch (ele%value(x_pitch_tot$), ele%value(y_pitch_tot$), ele%orientation, mat6)

! 1/gamma^2 m56 correction

mass = mass_of(param%particle)
e_tot = ele%value(p0c$) * (1 + c0%vec(6)) / c0%beta
mat6(5,6) = mat6(5,6) + length * mass**2 * ele%value(e_tot$) / e_tot**3

ele%vec0 = c1%vec - matmul(mat6, c0%vec)
ele%value(tilt_tot$) = 0

end subroutine make_mat6_sad_mult

end module
