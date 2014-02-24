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
real(rp) ks, k1, length, z_start, charge_dir
real(rp) xp_start, yp_start, mat4(4,4), mat1(6,6), f1, f2, ll
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(rp), pointer :: mat6(:,:)
real(rp) :: vec0(6)

integer n, nd, orientation, n_div, np_max, physical_end, fringe_at

logical make_matrix, end_in, has_nonzero

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
charge_dir = param%rel_tracking_charge * orientation

if (make_matrix) then
  mat6 => ele%mat6
  call mat_make_unit(mat6)
endif

call multipole_ele_to_kt (ele, param, .true., has_nonzero, knl, tilt)

! If element has zero length then the SAD ignores f1 and f2.

if (length == 0) then
  call offset_particle (ele, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false.)
  call multipole_kicks (knl, tilt, orbit)
  call offset_particle (ele, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false.)
  if (make_matrix) then
    call multipole_kick_mat (knl, tilt, orbit%vec, 1.0_rp, mat6)
  endif
  return
endif

!

ks = param%rel_tracking_charge * ele%value(ks$)
k1 = charge_dir * knl(1) / length

f1 = -k1 * ele%value(f1$) * abs(ele%value(f1$)) / (24 * rel_pc)
f2 =  k1 * ele%value(f2$) / rel_pc

call transfer_ele(ele, ele2)
ele2%value(tilt_tot$) = tilt(1)
ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) - ele%value(x_offset_sol$)

tilt = tilt - tilt(1)
knl(1) = 0
knl = knl / n_div

call offset_particle (ele2, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false.)

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (first_track_edge$, orbit%direction, ele%orientation)
if (at_this_ele_end(physical_end, fringe_at)) then
  if (ele%value(fringe_type$) == full_sad$ .or. ele%value(fringe_type$) == linear_sad$) then
    call linear_fringe_kick (f1, f2)
  endif
  if (ele%value(fringe_type$) == full_sad$ .or. ele%value(fringe_type$) == nonlin_only_sad$) then
    call nonlinear_fringe_kick ()
  endif
endif

! Body

do nd = 0, n_div

  ll = length / n_div
  if (nd == 0 .or. nd == n_div) ll = ll / 2

  ! Matrix step

  if (make_matrix) then
    if (k1 == 0) then
      call solenoid_mat6_calc (ks, ll, 0.0_rp, orbit, mat1)
    else
      call sol_quad_mat6_calc (ks, k1, ll, orbit%vec, mat1)
    endif
    mat6 = matmul(mat1, mat6)
  endif

  ! track step

  if (k1 == 0) then
    xp_start = orbit%vec(2) + ks * orbit%vec(3) / (2 * rel_pc)
    yp_start = orbit%vec(4) - ks * orbit%vec(1) / (2 * rel_pc)
    call solenoid_mat4_calc (ks, ll, rel_pc, mat4)
    orbit%vec(5) = orbit%vec(5) - ll * (xp_start**2 + yp_start**2 ) / 2
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

  call multipole_kicks (knl/rel_pc, tilt, orbit)

  if (make_matrix) then
    call multipole_kick_mat (knl, tilt, orbit%vec, 1.0_rp, mat1)
    mat6(2,:) = mat6(2,:) + mat1(2,1) * mat6(1,:) + mat1(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat1(4,1) * mat6(1,:) + mat1(4,3) * mat6(3,:)
  endif

enddo

! End stuff

physical_end = physical_ele_end (second_track_edge$, orbit%direction, ele%orientation)
if (at_this_ele_end(physical_end, fringe_at)) then
  if (ele%value(fringe_type$) == full_sad$ .or. ele%value(fringe_type$) == linear_sad$) then
    call linear_fringe_kick (-f1, f2)
  endif
  if (ele%value(fringe_type$) == full_sad$ .or. ele%value(fringe_type$) == nonlin_only_sad$) then
    call nonlinear_fringe_kick ()
  endif
endif

call offset_particle (ele2, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false.)

if (make_matrix) then
  if (ele2%value(tilt_tot$) /= 0) then
    call tilt_mat6 (mat6, ele2%value(tilt_tot$))
  endif

  call mat6_add_pitch (ele2%value(x_pitch_tot$), ele2%value(y_pitch_tot$), ele2%orientation, mat6)

  ! 1/gamma^2 m56 correction

  mass = mass_of(param%particle)
  e_tot = ele%value(p0c$) * (1 + orbit%vec(6)) / orbit%beta
  mat6(5,6) = mat6(5,6) + length * mass**2 * ele%value(e_tot$) / e_tot**3

  ele%vec0 = orbit%vec - matmul(mat6, start_orb%vec)
endif

!

call track1_low_energy_z_correction (orbit, ele2, param)

orbit%t = start_orb%t + (length + start_orb%vec(5) - orbit%vec(5)) / (orbit%beta * c_light)
orbit%s = ele%s

if (.not. end_in) end_orb = orbit

!----------------------------------------------------------------------------------------------
contains

subroutine linear_fringe_kick (f1, f2)

real(rp) f1, f2, ef1, mat1(6,6), vec(4)

! Notice that orbit%vec(:) is in (x', y') units.

if (f1 == 0 .and. f2 == 0) return

ef1 = exp(f1)
vec = orbit%vec(1:4)


orbit%vec(5) = orbit%vec(5) - (f1 * vec(1) + f2 * (1 + f1/2) * vec(2) / ef1) * vec(2) + &
                              (f1 * vec(3) + f2 * (1 - f1/2) * vec(4) * ef1) * vec(4)

orbit%vec(1:2) = [vec(1) * ef1 + vec(2) * f2, vec(2) / ef1]
orbit%vec(3:4) = [vec(3) / ef1 - vec(4) * f2, vec(4) * ef1]

!

if (make_matrix) then
  mat1 = 0
  mat1(1,1) = ef1
  mat1(1,2) = f2 / rel_pc 
  mat1(1,6) = -vec(1) * f1 * ef1 / rel_pc - 2 * vec(2) * f2 / rel_pc
  mat1(2,2) = 1 / ef1
  mat1(2,6) = vec(2) * f1 / ef1 
  mat1(3,3) = 1 / ef1
  mat1(3,4) = -f2 / rel_pc
  mat1(3,6) =  vec(3) * f1 / ef1 / rel_pc + 2 * vec(4) * f2 / rel_pc
  mat1(4,4) = ef1
  mat1(4,6) = -vec(4) * f1 * ef1
  mat1(5,1) = -f1 * vec(2)
  mat1(5,2) = -(f1 * vec(1) + f2 * (2 + f1) * vec(2) / ef1) / rel_pc
  mat1(5,3) =  f1 * vec(4)
  mat1(5,4) =  (f1 * vec(3) + f2 * (2 - f1) * vec(4) * ef1) / rel_pc
  mat1(5,5) = 1
  mat1(5,6) = 2 * f1 * vec(1) * vec(2) / rel_pc + &
              f2 * (1 + f1/2) * vec(2)**2 * (3 - f1) / (ef1 * rel_pc) + & 
              f2 * f1 * vec(2)**2 / (2 * ef1 * rel_pc) - &
              2 * f1 * vec(3) * vec(4) / rel_pc - &
              f2 * (1 - f1/2) * vec(4)**2 * ef1 * (3 + f1) / rel_pc + & 
              f2 * f1 * vec(4)**2 * ef1 / (2 * rel_pc)

  mat1(6,6) = 1

  mat6 = matmul(mat1, mat6)
endif

end subroutine linear_fringe_kick

!----------------------------------------------------------------------------------------------
! contains

subroutine nonlinear_fringe_kick ()

end subroutine nonlinear_fringe_kick

end subroutine sad_mult_track_and_mat 

end module
