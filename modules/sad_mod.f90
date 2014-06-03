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
real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
real(rp), pointer :: mat6(:,:)
real(rp) :: vec0(6), kmat(6,6)

integer n, nd, orientation, n_div, np_max, physical_end, fringe_at

logical make_matrix, end_in, has_nonzero, fringe_here

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
mat6 => ele%mat6

if (make_matrix) then
  call mat_make_unit(mat6)
endif

call multipole_ele_to_kt (ele, param, .true., has_nonzero, knl, tilt)

! Setup ele2 which is used in offset_particle

call transfer_ele(ele, ele2)
ele2%value(x_pitch_tot$) = ele%value(x_pitch_tot$) + ele%value(x_pitch_mult$)
ele2%value(y_pitch_tot$) = ele%value(y_pitch_tot$) + ele%value(y_pitch_mult$)
ele2%value(x_offset_tot$) = ele%value(x_offset_tot$) + ele%value(x_offset_mult$)
ele2%value(y_offset_tot$) = ele%value(y_offset_tot$) + ele%value(y_offset_mult$)

! If element has zero length then the SAD ignores f1 and f2.

if (length == 0) then
  call offset_particle (ele2, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false.)
  call multipole_kicks (knl, tilt, orbit)
  if (make_matrix) then
    call multipole_kick_mat (knl, tilt, orbit%vec, 1.0_rp, mat6)
  endif
  call offset_particle (ele2, param, unset$, orbit, set_multipoles = .false., set_hvkicks = .false., set_tilt = .false.)

  if (make_matrix) then
    call mat6_add_pitch (ele2%value(x_pitch_tot$), ele2%value(y_pitch_tot$), ele2%orientation, mat6)
    ele%vec0 = orbit%vec - matmul(mat6, start_orb%vec)
  endif

  orbit%s = ele%s
  if (.not. end_in) end_orb = orbit

  return
endif

! Go to frame of reference of the multipole

ks = param%rel_tracking_charge * ele%value(ks$)
k1 = charge_dir * knl(1) / length
knl(1) = 0

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

call offset_particle (ele2, param, set$, orbit, set_multipoles = .false., set_hvkicks = .false.)

! Quadrupole edge kick

if (make_matrix) then
  call quadrupole_edge_mat6 (ele, first_track_edge$, orbit, kmat, fringe_here)
  if (fringe_here) mat6 = kmat
endif
call quadrupole_edge_kick (ele, first_track_edge$, orbit)

! Dipole edge kick

k0 = knl(0)/length
call sad_linear_dipole_edge (ele, k0, tilt(0), first_track_edge$, orbit, mat6, make_matrix)

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
    call multipole_kick_mat (knl, tilt, orbit%vec, 1.0_rp, mat1)
    mat6(2,:) = mat6(2,:) + mat1(2,1) * mat6(1,:) + mat1(2,3) * mat6(3,:)
    mat6(4,:) = mat6(4,:) + mat1(4,1) * mat6(1,:) + mat1(4,3) * mat6(3,:)
  endif

  ! Check for orbit too large to prevent infinities.

  if (orbit_too_large (orbit)) then
    if (.not. end_in) end_orb = orbit
    return
  endif

enddo

! End stuff

! Dipole edge kick

call sad_linear_dipole_edge (ele, k0, tilt(0), second_track_edge$, orbit, mat6, make_matrix)

! Quadrupole edge kick

if (make_matrix) then
  call quadrupole_edge_mat6 (ele, second_track_edge$, orbit, kmat, fringe_here)
  if (fringe_here) mat6 = matmul(kmat, mat6)
endif
call quadrupole_edge_kick (ele, second_track_edge$, orbit)

!

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

if (.not. end_in) end_orb = orbit

end subroutine sad_mult_track_and_mat 

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine sad_linear_dipole_edge (ele, g_bend, tilt, particle_at, orbit, mat6, make_matrix)
!
! Routine to track through the "linear" dipole fringe field.
!
! Input:
!   ele         -- ele_struct: Element with fringe.
!   g_bend      -- real(rp): Dipole bend strength.
!   tilt        -- real(rp): field rotation.
!   particle_at -- integer: Either first_track_edge$ or second_track_edge$.
!   orbit       -- coord_struct: Starting coordinates.
!   mat6(6,6)   -- real(rp): Transfer matrix up to the fringe.
!   make_matrix -- logical: Make the transfer matrix?
!
! Output:
!   orbit       -- coord_struct: Ending coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the fringe.
!-

subroutine sad_linear_dipole_edge (ele, g_bend, tilt, particle_at, orbit, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orbit

real(rp) :: mat6(6,6)
real(rp) g_bend, tilt
real(rp) g, px, y, y2, rel_p, p_zy, yg, kmat(6,6)

integer fringe_type, fringe_at, physical_end, particle_at
logical make_matrix

! Fringe here?

fringe_type = nint(ele%value(fringe_type$))
if (fringe_type /= sad_linear$ .and. fringe_type /= sad_full$) return

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end(physical_end, fringe_at)) return

! Rotate

if (tilt /= 0) call tilt_coords (tilt, orbit%vec)

! edge kick

g = g_bend
if (particle_at == second_track_edge$) g = -g

px = orbit%vec(2)
y = orbit%vec(3)
y2 = y**2
rel_p = 1 + orbit%vec(6)
p_zy = sqrt(rel_p**2 - px**2)
yg = y2 * g**2 / 12

orbit%vec(1) = orbit%vec(1) + g * y2 * (1 - yg) * rel_p**2 / (2 * p_zy**3)
orbit%vec(4) = orbit%vec(4) - g * px * y * (1 - 2 * yg) / p_zy
orbit%vec(5) = orbit%vec(5) - g * y2 * px * (1 - yg) * rel_p / (2 * p_zy**3)

if (make_matrix) then
  call mat_make_unit(kmat)
  kmat(1,2) = 3 * g * px * y2 * (1 - yg) * rel_p**2 / (2 * p_zy**5)
  kmat(1,3) =  g * (y - 2 * y * yg) * rel_p**2 / p_zy**3
  kmat(1,6) = -g * y2 * (1 - yg) * rel_p * (rel_p**2 + 2 * px**2) / (2 * p_zy**5)
  kmat(4,2) = -g * y * (1 - 2 * yg) * rel_p**2 / p_zy**3
  kmat(4,3) = -g * px * (1 - 6 * yg) / p_zy
  kmat(4,6) =  g * px * y * (1 - 2 * yg) * rel_p / p_zy**3
  kmat(5,2) = -g * y2 * (1 - yg) * (rel_p**2 + 2 * px**2) * rel_p / (2 * p_zy**5)
  kmat(5,3) = -g * px * (y - 2 * y * yg) * rel_p / p_zy**3
  kmat(5,6) =  g * px * y2 * (1 - 2 * yg) * (2 * rel_p**2 + px**2) / (2 * p_zy**5)
  call tilt_mat6(kmat, tilt)
  mat6 = matmul(kmat, mat6)
endif

! Rotate

if (tilt /= 0) call tilt_coords (-tilt, orbit%vec)

end subroutine sad_linear_dipole_edge

end module
