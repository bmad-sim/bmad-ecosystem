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
charge_dir = param%rel_tracking_charge * orientation
mat6 => ele%mat6

if (make_matrix) call mat_make_unit(mat6)

knl = 0; tilt = 0
call multipole_ele_to_kt (ele, param, .true., has_nonzero_pole, knl, tilt)

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
      call multipole_kick_mat (knl, tilt, orbit%vec, 1.0_rp, mat6)
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

ks = param%rel_tracking_charge * ele%value(ks$)
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

call multipole_hard_edge_kick (ele, a_pole, b_pole, first_track_edge$, orbit, mat6, make_matrix)
call quadrupole_soft_edge_kick (ele, first_track_edge$, orbit, mat6, make_matrix)
call sad_bend_linear_edge_kick (ele, k0, tilt(0), first_track_edge$, orbit, mat6, make_matrix)
call sad_bend_soft_edge_kick (ele, param, first_track_edge$, orbit, mat6, make_matrix, k0, tilt(0))

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

! Exit edge kicks

call sad_bend_soft_edge_kick (ele, param, second_track_edge$, orbit, mat6, make_matrix, k0, tilt(0))
call sad_bend_linear_edge_kick (ele, k0, tilt(0), second_track_edge$, orbit, mat6, make_matrix)
call quadrupole_soft_edge_kick (ele, second_track_edge$, orbit, mat6, make_matrix)
call multipole_hard_edge_kick (ele, a_pole, b_pole, second_track_edge$, orbit, mat6, make_matrix)

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

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine sad_bend_linear_edge_kick (ele, g_bend, tilt, particle_at, orbit, mat6, make_matrix)
!
! Routine to track through the "linear" bend fringe field.
!
! Input:
!   ele         -- ele_struct: Element with fringe.
!   g_bend      -- real(rp): Bend bend strength.
!   tilt        -- real(rp): field rotation.
!   particle_at -- integer: Either first_track_edge$ or second_track_edge$.
!   orbit       -- coord_struct: Starting coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix up to the fringe.
!   make_matrix -- real(rp), optional: Make the transfer matrix? Default is False.
!
! Output:
!   orbit       -- coord_struct: Ending coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the fringe.
!-

subroutine sad_bend_linear_edge_kick (ele, g_bend, tilt, particle_at, orbit, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orbit

real(rp), optional :: mat6(6,6)
real(rp) g_bend, tilt
real(rp) g, px, y, y2, rel_p, p_zy, yg, kmat(6,6)

integer fringe_type, fringe_at, physical_end, particle_at
integer i, i_max

logical, optional :: make_matrix

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

if (logic_option(.false., make_matrix)) then
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

end subroutine sad_bend_linear_edge_kick

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!+
! Subroutine multipole_hard_edge_kick (ele, a_pole, b_pole, particle_at, orbit, mat6, make_matrix)
!
! Routine to track through the hard edge field of a multipole.
! The dipole component is ignored and only quadrupole and higher multipoles are included.
!
! Input:
!   ele         -- ele_struct: Element with fringe.
!   a_pole(0:)  -- real(rp): Multipole skew components.
!   b_pole(0:)  -- real(rp): Multipole normal components.
!   particle_at -- integer: Either first_track_edge$ or second_track_edge$.
!   orbit       -- coord_struct: Starting coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix up to the fringe.
!   make_matrix -- real(rp), optional: Make the transfer matrix? Default is False.
!
! Output:
!   orbit       -- coord_struct: Ending coordinates.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix including the fringe.
!-

subroutine multipole_hard_edge_kick (ele, a_pole, b_pole, particle_at, orbit, mat6, make_matrix)

implicit none

type (ele_struct) ele
type (coord_struct) orbit

real(rp), optional :: mat6(6,6)
real(rp) a_pole(0:), b_pole(0:)
real(rp) rel_p, cn, x, y, px, py, denom, ddenom_dx, ddenom_dy, ddenom_dpz, kmat(6,6)
real(rp) fx, dfx_dx, dfx_dy, d2fx_dxx, d2fx_dxy, d2fx_dyy
real(rp) fy, dfy_dx, dfy_dy, d2fy_dxx, d2fy_dxy, d2fy_dyy

complex(rp) poly, poly_n1, poly_n2, dpoly_dx, dpoly_dy, d2poly_dxx, d2poly_dxy, d2poly_dyy
complex(rp) xy, xny, dxny_dx, dxny_dy, cab

integer fringe_type, fringe_at, physical_end, particle_at
integer n, n_max

logical, optional :: make_matrix

! Fringe here?

fringe_type = nint(ele%value(fringe_type$))
if (fringe_type /= sad_full$ .and. fringe_type /= sad_nonlin_only$) return

fringe_at = nint(ele%value(fringe_at$))
physical_end = physical_ele_end (particle_at, orbit%direction, ele%orientation)
if (.not. at_this_ele_end(physical_end, fringe_at)) return

!

do n = ubound(a_pole, 1), 0, -1
  if (a_pole(n) /= 0 .or. b_pole(n) /= 0) exit
enddo
n_max = n

!

x = orbit%vec(1)
y = orbit%vec(3)
xy = cmplx(x, y)

poly_n1 = 1
poly = xy

rel_p = 1 + orbit%vec(6)

fx = 0
dfx_dx = 0
dfx_dy = 0
d2fx_dxx = 0
d2fx_dxy = 0
d2fx_dyy = 0

fy = 0
dfy_dx = 0
dfy_dy = 0
d2fy_dxx = 0
d2fy_dxy = 0
d2fy_dyy = 0

do n = 1, n_max

  poly_n2 = poly_n1
  poly_n1 = poly
  poly = poly * xy

  if (a_pole(n) == 0 .and. b_pole(n) == 0) cycle

  dpoly_dx = (n+1) * poly_n1
  dpoly_dy = i_imaginary * dpoly_dx

  d2poly_dxx = n * (n+1) * poly_n2
  d2poly_dxy = i_imaginary * d2poly_dxx
  d2poly_dyy = -d2poly_dxx

  cab = cmplx(b_pole(n), a_pole(n)) / (4 * (n + 2) * rel_p * ele%value(l$))
  if (particle_at == first_track_edge$) cab = -cab
  cn = real(n+3, rp) / (n+1) 

  xny = cmplx(x, -cn * y)
  dxny_dy = cmplx(0.0_rp, -cn)

  fx = fx + real(cab * poly * xny)
  dfx_dx = dfx_dx + real(cab * (dpoly_dx * xny + poly))
  dfx_dy = dfx_dy + real(cab * (dpoly_dy * xny + poly * dxny_dy))
  d2fx_dxx = d2fx_dxx + real(cab * (d2poly_dxx * xny + 2 * dpoly_dx))
  d2fx_dxy = d2fx_dxy + real(cab * (d2poly_dxy * xny + dpoly_dx * dxny_dy + dpoly_dy))
  d2fx_dyy = d2fx_dyy + real(cab * (d2poly_dyy * xny + 2 * dpoly_dy * dxny_dy))

  xny = cmplx(y, cn * x)
  dxny_dx = cmplx(0.0_rp, cn)

  fy = fy + real(cab * poly * xny)
  dfy_dx = dfy_dx + real(cab * (dpoly_dx * xny + poly * dxny_dx))
  dfy_dy = dfy_dy + real(cab * (dpoly_dy * xny + poly))
  d2fy_dxx = d2fy_dxx + real(cab * (d2poly_dxx * xny + 2 * dpoly_dx * dxny_dx))
  d2fy_dxy = d2fy_dxy + real(cab * (d2poly_dxy * xny + dpoly_dx + dpoly_dy * dxny_dx))
  d2fy_dyy = d2fy_dyy + real(cab * (d2poly_dyy * xny + 2 * dpoly_dy))
enddo

px = orbit%vec(2)
py = orbit%vec(4)
denom = (1 - dfx_dx) * (1 - dfy_dy) - dfx_dy * dfy_dx
ddenom_dx = -d2fx_dxx - d2fy_dxy + d2fx_dxx * dfy_dy + dfx_dx * d2fy_dxy - d2fx_dxy * dfy_dx - dfx_dy * d2fy_dxx 
ddenom_dy = -d2fx_dxy - d2fy_dyy + d2fx_dxy * dfy_dy + dfx_dx * d2fy_dyy - d2fx_dyy * dfy_dx - dfx_dy * d2fy_dxy 
ddenom_dpz = (dfx_dx + dfy_dy - 2 * dfx_dx * dfy_dy + 2 * dfx_dy * dfy_dx) / rel_p

orbit%vec(1) = orbit%vec(1) - fx
orbit%vec(2) = ((1 - dfy_dy) * px + dfy_dx * py) / denom
orbit%vec(3) = orbit%vec(3) - fy
orbit%vec(4) = (dfx_dy * px + (1 - dfx_dx) * py) / denom
orbit%vec(5) = orbit%vec(5) + (orbit%vec(2) * fx + orbit%vec(4) * fy ) / rel_p

if (logic_option(.false., make_matrix)) then
  kmat = 0
  kmat(1,1) = 1 - dfx_dx
  kmat(1,3) = -dfx_dy
  kmat(1,6) = fx / rel_p
  kmat(2,1) = (-d2fy_dxy * px + d2fy_dxx * py) / denom - orbit%vec(2) * ddenom_dx / denom
  kmat(2,2) = (1 - dfy_dy) / denom
  kmat(2,3) = (-d2fy_dyy * px + d2fy_dxy * py) / denom - orbit%vec(2) * ddenom_dy / denom
  kmat(2,4) = dfy_dx / denom
  kmat(2,6) = (dfy_dy * px - dfy_dx * py) / (denom * rel_p) - orbit%vec(2) * ddenom_dpz / denom
  kmat(3,1) = -dfy_dx
  kmat(3,3) = 1 - dfy_dy
  kmat(3,6) = fy / rel_p
  kmat(4,1) = (d2fx_dxy * px - d2fx_dxx * py) / denom - orbit%vec(4) * ddenom_dx / denom
  kmat(4,2) = dfx_dy / denom
  kmat(4,3) = (d2fx_dyy * px - d2fx_dxy * py) / denom - orbit%vec(4) * ddenom_dy / denom
  kmat(4,4) = (1 - dfx_dx) / denom
  kmat(4,6) = (-dfx_dy * px + dfx_dx * py) / (denom * rel_p) - orbit%vec(4) * ddenom_dpz / denom
  kmat(5,1) = (kmat(2,1) * fx + orbit%vec(2) * dfx_dx + kmat(4,1) * fy + orbit%vec(4) *dfy_dx) / rel_p
  kmat(5,2) = (kmat(2,2) * fx + kmat(4,2) * fy) / rel_p
  kmat(5,3) = (kmat(2,3) * fx + orbit%vec(2) * dfx_dy + kmat(4,3) * fy + orbit%vec(4) *dfy_dy) / rel_p
  kmat(5,4) = (kmat(2,4) * fx + kmat(4,4) * fy) / rel_p
  kmat(5,5) = 1
  kmat(5,6) = (kmat(2,6) * fx + kmat(4,6) * fy) / rel_p - 2 * (orbit%vec(2) * fx + orbit%vec(4) * fy) / rel_p**2
  kmat(6,6) = 1
  mat6 = matmul (kmat, mat6)
endif

end subroutine multipole_hard_edge_kick 

end module
