module sad_mod

use track1_mod
use make_mat6_mod

contains

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine track_a_sad_mult (start_orb, ele, param, end_orb)
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
type (ele_struct)  :: ele
type (lat_param_struct) :: param

real(rp) rel_pc, vec0(6), dz4_coef(4,4), mat6(6,6)
real(rp) ks, k1, length, z_start, charge_dir

integer orientation

!

start_orb = orbit

length = ele%value(l$)
rel_pc = 1 + orbit%vec(6)
orientation = ele%orientation * start_orb%direction
charge_dir = param%rel_tracking_charge * orientation

call offset_particle (ele, orbit, param, set$, set_multipoles = .false., set_hvkicks = .false.)

ks = param%rel_tracking_charge * ele%value(ks$) / rel_pc
k1 = charge_dir * ele%value(k1$) / rel_pc

vec0 = 0
call sol_quad_mat6_calc (ks, k1, length, mat6, vec0, dz4_coef)
orbit%vec(5) = orbit%vec(5) + sum(orbit%vec(1:4) * matmul(dz4_coef, orbit%vec(1:4)))   
orbit%vec(1:4) = matmul (mat6(1:4,1:4), orbit%vec(1:4))

call offset_particle (ele, orbit, param, unset$, set_multipoles = .false., set_hvkicks = .false.)

call track1_low_energy_z_correction (orbit, ele, param)

orbit%t = start_orb%t + (ele%value(l$) + start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)
orbit%s = ele%s

end subroutine track_a_sad_mult

end module
