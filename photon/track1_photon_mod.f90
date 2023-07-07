module track1_photon_mod

use photon_utils_mod
use xraylib_interface

implicit none

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_lens (ele, param, orbit)
!
! Routine to track through a lens.
!
! Input:
!   ele      -- ele_struct: Element tracking through.
!   param    -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_lens (ele, param, orbit)

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

character(*), parameter :: r_name = 'track1_lens'

!

! Nothing implemented yet.

end subroutine track1_lens

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track_a_patch_photon (ele, orbit, drift_to_exit, use_z_pos)
!
! Routine to track through a patch element with a photon.
! The steps for tracking are:
!   1) Transform from entrance to exit coordinates.
!   2) Drift particle from the entrance to the exit coordinants.
!
! Input:
!   ele           -- ele_struct: patch element.
!   orbit         -- coord_struct: Starting phase space coords
!   drift_to_exit -- logical, optional: If False then do not drift the particle from
!                      start to ending faces. Default is True.
!   use_z_pos     -- loigical, optional: If present and True, use orbit%vec(5) as the true
!                      z-position relative to the start of the element instead of assuming 
!                      that the particle is at the patch edge.
!
! Output:
!   orbit         -- coord_struct: Coords after applying a patch transformation.
!-

Subroutine track_a_patch_photon (ele, orbit, drift_to_exit, use_z_pos)

type (coord_struct) orbit
type (ele_struct) ele

real(rp) w(3,3)

logical, optional :: drift_to_exit, use_z_pos

!

if (orbit%direction == 1) then
  ! Translate (x, y, z) to coordinate system with respect to downstream origin.
  orbit%vec(1) = orbit%vec(1) - ele%value(x_offset$)
  orbit%vec(3) = orbit%vec(3) - ele%value(y_offset$)
  if (logic_option(.false., use_z_pos)) then
    orbit%vec(5) = orbit%vec(5) - ele%value(z_offset$)   
  else
    orbit%vec(5) = -ele%value(z_offset$)   
  endif

  if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0) then
    call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat_inv = w)
    orbit%vec(2:6:2) = matmul(w, orbit%vec(2:6:2))
    orbit%vec(1:5:2) = matmul(w, orbit%vec(1:5:2))
  endif

  if (logic_option(.true., drift_to_exit)) then
    call track_a_drift_photon (orbit, -orbit%vec(5), .false.)
  endif

  orbit%vec(5) = orbit%vec(5) + ele%value(l$)
  orbit%s = ele%s_start + orbit%vec(5)

else
  
  ! Shift to z being with respect to exit end coords.
  if (logic_option(.false., use_z_pos)) then
    orbit%vec(5) = orbit%vec(5) - ele%value(l$)
  else
    orbit%vec(5) = 0   ! Assume particle starts at downstream face
  endif

  if (ele%value(x_pitch$) /= 0 .or. ele%value(y_pitch$) /= 0 .or. ele%value(tilt$) /= 0) then
    call floor_angles_to_w_mat (ele%value(x_pitch$), ele%value(y_pitch$), ele%value(tilt$), w_mat = w)
    orbit%vec(2:6:2) = matmul(w, orbit%vec(2:6:2))
    orbit%vec(1:5:2) = matmul(w, orbit%vec(1:5:2))
  endif

  orbit%vec(1) = orbit%vec(1) + ele%value(x_offset$)
  orbit%vec(3) = orbit%vec(3) + ele%value(y_offset$)
  orbit%vec(5) = orbit%vec(5) + ele%value(z_offset$)

  if (logic_option(.true., drift_to_exit)) then
    call track_a_drift_photon (orbit, -orbit%vec(5), .false.)
  endif

  orbit%s = ele%s_start + orbit%vec(5)

endif

end subroutine track_a_patch_photon

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_diffraction_plate_or_mask (ele, param, orbit)
!
! Routine to track through diffraction plate and mask elements.
!
! Input:
!   ele      -- ele_struct: Diffraction plate or mask element.
!   param    -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_diffraction_plate_or_mask (ele, param, orbit)

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param
type (wall3d_section_struct), pointer :: sec

real(rp) vz0
real(rp) wavelength, thickness, absorption, phase_shift

integer ix_sec

logical err_flag

character(60) material

! If the plate/mask is turned off then all photons are simply transmitted through.

if (.not. ele%is_on) return

! Photon is lost if in an opaque section

ix_sec = diffraction_plate_or_mask_hit_spot (ele, orbit)

if (ix_sec == 0) then
  if (ele%wall3d(1)%opaque_material == '') then
    orbit%state = lost$
    return
  endif
  material = ele%wall3d(1)%opaque_material
  thickness = ele%wall3d(1)%thickness

else
  material = ele%wall3d(1)%clear_material
  thickness = ele%wall3d(1)%thickness
  sec => ele%wall3d(1)%section(ix_sec)
  if (sec%material /= '') material = sec%material
  if (sec%thickness >= 0) thickness = sec%thickness
endif

! Transmit through material

if (material /= '') then
  call photon_absorption_and_phase_shift (material, orbit%p0c, absorption, phase_shift, err_flag)
  if (err_flag) then
    orbit%state = lost$
    return
  endif

  orbit%field = orbit%field * exp(-absorption * thickness)
  orbit%phase = orbit%phase - phase_shift * thickness
endif

! Choose outgoing direction

vz0 = orbit%vec(6)
if (ele%key == diffraction_plate$) call point_photon_emission (ele, param, orbit, +1, twopi)

! Rescale field

wavelength = c_light * h_planck / orbit%p0c
orbit%field = orbit%field * (vz0 + orbit%vec(6)) / (2 * wavelength)
orbit%phase = orbit%phase - pi / 2

if (ele%value(field_scale_factor$) /= 0) then
  orbit%field = orbit%field * ele%value(field_scale_factor$)
endif

end subroutine track1_diffraction_plate_or_mask

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_sample (ele, param, orbit)
!
! Routine to track reflection from a sample element.
!
! Input:
!   ele      -- ele_struct: Element tracking through.
!   param    -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_sample (ele, param, orbit)

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

real(rp) w_surface(3,3), absorption, phase_shift

logical err_flag

character(*), parameter :: r_name = 'track1_sample'

!

call track_to_surface (ele, orbit, param, w_surface)
if (orbit%state /= alive$) return

! Reflect 

select case (nint(ele%value(mode$)))
case (reflection$)

  call point_photon_emission (ele, param, orbit, -1, fourpi, w_surface)

case (transmission$)
  call photon_absorption_and_phase_shift (ele%component_name, orbit%p0c, absorption, phase_shift, err_flag)
  if (err_flag) then
    orbit%state = lost$
    return
  endif

  orbit%field = orbit%field * exp(-absorption * ele%value(l$))
  orbit%phase = orbit%phase - phase_shift * ele%value(l$)

case default
  call out_io (s_error$, r_name, 'MODE NOT SET.')
  orbit%state = lost$
end select

! Rotate back to uncurved element coords

if (has_curvature(ele%photon)) call rotate_for_curved_surface (ele, orbit, unset$, w_surface)

end subroutine track1_sample

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine point_photon_emission (ele, param, orbit, direction, max_target_area, w_to_surface)
!
! Routine to emit a photon from a point that may be on a surface.
! If there is a downstream target, the emission calc will take this into account.
!
! Input:
!   ele                -- ele_struct: Emitting element.
!   param              -- lat_param_struct: lattice parameters.
!   orbit              -- Coord_struct: phase-space coords of photon.
!                      --   Will be in curved surface coords if there is a curved surface.
!   direction          -- Integer: +1 -> Emit in forward +z direction, -1 -> emit backwards.
!   max_target_area    -- real(rp): Area of the solid angle photons may be emitted over.
!                           max_target_area is used for normalizing the photon field.
!                           generally will be equal to twopi or fourpi.
!   w_to_surface(3,3)  -- real(rp), optional: Rotation matrix for curved surface.
!
! Output:
!   orbit    -- Coord_struct: Final phase-space coords
!-

subroutine point_photon_emission (ele, param, orbit, direction, max_target_area, w_to_surface)

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (ele_struct), pointer :: det_ele
type (lat_param_struct) :: param
type (photon_target_struct), pointer :: target
type (target_point_struct) corner(8)
type (surface_grid_struct), pointer :: gr

real(rp), optional :: w_to_surface(3,3)
real(rp) zran(2), r_particle(3), w_to_target(3,3), w_to_ele(3,3)
real(rp) phi_min, phi_max, y_min, y_max, y, phi, rho, r(3), max_target_area, cos2_dphi
real(rp) lb(2), ub(2), r_len, dr_x(3), dr_y(3), area

integer direction
integer n, i, ix, iy

!

target => ele%photon%target

select case (target%type)
case (rectangular$)
  r_particle = orbit%vec(1:5:2)
  r = target%center%r - r_particle
  if (ele%photon%curvature%has_curvature) r = matmul(w_to_surface, r)
  call target_rot_mats (r, w_to_target, w_to_ele)

  do i = 1, target%n_corner
    r = target%corner(i)%r - r_particle
    if (direction == 1) then
      if (r(3) < 0) r(3) = 0   ! photon cannot be emitted backward
    elseif (direction == -1) then
      if (r(3) > 0) r(3) = 0   ! photon cannot be emitted forward
    else
      call err_exit
    endif
    if (ele%photon%curvature%has_curvature) r = matmul(w_to_surface, r)
    r = matmul(w_to_target, r)
    corner(i)%r = r / norm2(r)
  enddo

  call target_min_max_calc (corner(4)%r, corner(1)%r, y_min, y_max, phi_min, phi_max, .true.)
  call target_min_max_calc (corner(1)%r, corner(2)%r, y_min, y_max, phi_min, phi_max)
  call target_min_max_calc (corner(2)%r, corner(3)%r, y_min, y_max, phi_min, phi_max)
  call target_min_max_calc (corner(3)%r, corner(4)%r, y_min, y_max, phi_min, phi_max)

  if (target%n_corner == 8) then
    call target_min_max_calc (corner(8)%r, corner(5)%r, y_min, y_max, phi_min, phi_max)
    call target_min_max_calc (corner(5)%r, corner(6)%r, y_min, y_max, phi_min, phi_max)
    call target_min_max_calc (corner(6)%r, corner(7)%r, y_min, y_max, phi_min, phi_max)
    call target_min_max_calc (corner(7)%r, corner(8)%r, y_min, y_max, phi_min, phi_max)
  endif

  if (y_min >= y_max .or. phi_min >= phi_max) then
    orbit%state = lost$
    return
  endif

  ! Correction for bulge in line projected onto (y, phi) sphere.

  cos2_dphi = cos((phi_max - phi_min)/2)**2

  if (y_max > 0) y_max = y_max / sqrt((1 - y_max**2) * cos2_dphi + y_max**2)
  if (y_min < 0) y_min = y_min / sqrt((1 - y_min**2) * cos2_dphi + y_min**2)

  !

  call ran_uniform(zran)
  y = y_min + (y_max-y_min) * zran(1)
  phi = phi_min + (phi_max-phi_min) * zran(2)
  rho = sqrt(1 - y*y)
  orbit%vec(2:6:2) = [rho * sin(phi), y, rho * cos(phi)]
  orbit%vec(2:6:2) = matmul(w_to_ele, orbit%vec(2:6:2))

  ! Field scaling

  if (photon_type(ele) == coherent$) then
    orbit%field = orbit%field * (y_max - y_min) * (phi_max - phi_min) / max_target_area
    ! If dt_ref = 0 then assume that photon is being initialized so only normalize field if dt_ref /= 0
    if (orbit%dt_ref /= 0) then
      orbit%field = orbit%field * orbit%dt_ref
      orbit%dt_ref = 0
    endif
  else
    orbit%field = orbit%field * sqrt ((y_max - y_min) * (phi_max - phi_min) / max_target_area)
  endif

! No targeting

case (off$)
  call ran_uniform(zran)
  y = 2 * zran(1) - 1
  phi = pi * (zran(2) - 0.5_rp)
  rho = sqrt(1 - y**2)
  orbit%vec(2:6:2) = [rho * sin(phi), y, direction * rho * cos(phi)]
  ! Without targeting photons are emitted into twopi solid angle.
  if (photon_type(ele) == coherent$) then
    orbit%field = orbit%field * twopi / max_target_area
    ! If dt_ref = 0 then assume that photon is being initialized so only normalize field if dt_ref /= 0
    if (orbit%dt_ref /= 0) then
      orbit%field = orbit%field * orbit%dt_ref
      orbit%dt_ref = 0
    endif
  else
    orbit%field = orbit%field * sqrt(twopi / max_target_area)
  endif

case default
  call err_exit  ! Internal bookkeeping error
end select

end subroutine point_photon_emission 

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_mirror (ele, param, orbit)
!
! Routine to track reflection from a mirror.
!
! Input:
!   ele      -- ele_struct: Element tracking through.
!   param    -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_mirror (ele, param, orbit)

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

real(rp) wavelength, w_surface(3,3)
real(rp), pointer :: val(:)

character(*), parameter :: r_name = 'track1_mirror'

!

val => ele%value
wavelength = c_light * h_planck / orbit%p0c

call track_to_surface (ele, orbit, param, w_surface)
if (orbit%state /= alive$) return

! Reflect momentum vector

orbit%vec(6) = -orbit%vec(6)

! Rotate back to uncurved element coords

if (has_curvature(ele%photon)) call rotate_for_curved_surface (ele, orbit, unset$, w_surface)

end subroutine track1_mirror

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_multilayer_mirror (ele, param, orbit)
!
! Routine to track reflection from a multilayer_mirror.
! Basic equations are from Kohn, "On the Theory of Reflectivity of an X-Ray Multilayer Mirror".
!
! Input:
!   ele    -- ele_struct: Element tracking through.
!   param  -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_multilayer_mirror (ele, param, orbit)

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

real(rp) wavelength, kz_air, w_surface(3,3)
real(rp), pointer :: val(:)

complex(rp) zero, xi_1, xi_2, kz1, kz2, c1, c2

character(*), parameter :: r_name = 'track1_multilayer_mirror'

!

val => ele%value
wavelength = c_light * h_planck / orbit%p0c

call track_to_surface (ele, orbit, param, w_surface)
if (orbit%state /= alive$) return

! Note: Koln z-axis = Bmad x-axis.
! Note: f0_re and f0_im are both positive.

xi_1 = -conjg(ele%photon%material%f0_m1) * r_e * wavelength**2 / (pi * val(v1_unitcell$)) 
xi_2 = -conjg(ele%photon%material%f0_m2) * r_e * wavelength**2 / (pi * val(v2_unitcell$)) 

kz1 = twopi * sqrt(orbit%vec(6)**2 + xi_1) / wavelength
kz2 = twopi * sqrt(orbit%vec(6)**2 + xi_2) / wavelength
kz_air = twopi * orbit%vec(6) / wavelength

c1 = exp(I_imaginary * kz1 * val(d1_thickness$) / 2)
c2 = exp(I_imaginary * kz2 * val(d2_thickness$) / 2)

zero = cmplx(0.0_rp, 0.0_rp, rp)

call multilayer_track (xi_1, xi_2, orbit%field(1), orbit%phase(1))     ! pi polarization
call multilayer_track (zero, zero, orbit%field(2), orbit%phase(2))     ! sigma polarization

! Reflect momentum vector

orbit%vec(6) = -orbit%vec(6)

! Rotate back to uncurved element coords

if (has_curvature(ele%photon)) call rotate_for_curved_surface (ele, orbit, unset$, w_surface)

!-----------------------------------------------------------------------------------------------
contains

subroutine multilayer_track (xi1_coef, xi2_coef, e_field, e_phase)

real(rp) e_field, e_phase

complex(rp) xi1_coef, xi2_coef, r_11, r_22, tt, denom, k1, k2
complex(rp) a, v, f_minus, f_plus, r_ratio, f, r, ttbar, nu, exp_half, exp_n
complex(rp) Rc_n_top, R_tot, k_a, r_aa, Rc_1_bot, Rc_1_top

integer i, n1

! Rc_n_top is the field ratio for the top layer of cell n.
! Rc_n_bot is the field ratio for the bottom layer of cell n.
! Rc_1_bot for the bottom layer just above the substrate is assumed to be zero.
! Upgrade: If we knew the substrate material we would not have to make this assumption.

Rc_1_bot = 0

! Compute Rc_1_top.
! The top layer of a cell is labeled "2" and the bottom "1". See Kohn Eq 6.

k1 = (1 + xi2_coef) * kz1
k2 = (1 + xi1_coef) * kz2
denom = k1 + k2

r_11 = c1**2 * (k1 - k2) / denom
r_22 = c2**2 * (k2 - k1) / denom

tt = 4 * k1 * k2 * (c1 * c2 / denom)**2    ! = t_12 * t_21

Rc_1_top = r_22 + tt * Rc_1_bot / (1 - r_11 * Rc_1_bot)

! Now compute the single cell factors. See Kohn Eq. 12.
! Note: If you go through the math you will find r = r_bar.

f = tt / (1 - r_11**2)
r = r_22 + r_11 * f
ttbar = f**2

! Calc Rc_n_top. See Kohn Eq. 21. Note that there are n-1 cells in between. 

a = (1 - ttbar + r**2) / 2
nu = (1 + ttbar - r**2) / (2 * sqrt(ttbar))
n1 = nint(val(n_cell$)) - 1

exp_half = nu + I_imaginary * sqrt(1 - nu**2)
exp_n = exp_half ** (2 * n1)
f_plus  = a - I_imaginary * sqrt(ttbar) * sqrt(1 - nu**2)
f_minus = a + I_imaginary * sqrt(ttbar) * sqrt(1 - nu**2)
Rc_n_top = r * (r - Rc_1_top * f_minus - (r - Rc_1_top * f_plus) * exp_n) / &
               (f_plus * (r - Rc_1_top * f_minus) - f_minus * (r - Rc_1_top * f_plus) * exp_n)

! Now factor in the air interface

k_a = kz_air
denom = k_a + k2

tt = 4 * k_a * k2 * (c2 / denom)**2
r_aa = (k_a - k2) / denom
r_22 = c2**2 * (k2 - k_a) / denom

R_tot = r_aa + tt * Rc_n_top / (1 - r_22 * Rc_n_top)

e_field = e_field * abs(R_tot)
e_phase = e_phase + atan2(aimag(R_tot), real(R_tot))

end subroutine multilayer_track 

end subroutine track1_multilayer_mirror

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_mosaic_crystal (ele, param, orbit)
!
! Routine to track diffraction from a crystal.
!
! Input:
!   ele      -- ele_struct: Element tracking through.
!   param    -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_mosaic_crystal (ele, param, orbit)


type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param
type (crystal_param_struct) cp

type mosaic_layer_struct
  real(rp) :: theta_in = 0, theta_out = 0
  real(rp) :: branch_ratio = 0
  logical :: follow_diffracted = .false.
end type
type (mosaic_layer_struct), allocatable, target :: m_layer(:)
type (mosaic_layer_struct), pointer :: ml

real(rp) h_bar0(3), h_bar_here(3), e_tot, pc, p_factor, w_surface(3,3), vec_init(6)
real(rp) gamma_0, gamma_h, dr1(3), dr2(3), dr3(3), dr4(3), denom, f, field_init
real(rp) field1, field2, field3, field4, phase1, phase2, phase3, phase4, rr, rr2(2)
real(rp) thick, p_sum, prob
real(rp), allocatable :: follow_prob(:)

integer n, im, h_dir, n_layer, n_choose
integer, allocatable :: indx(:)

character(*), parameter :: r_name = 'track1_mosaic_crystal'

! A graze angle of zero means the wavelength of the reference photon was too large
! for the bragg condition. 

if (orbit_too_large(orbit, param)) return

if (ele%value(b_param$) < 0) then
  call out_io (s_fatal$, r_name, 'MOSAIC CRYSTAL TRACKING IN BRAGG ORIENTATION NOT YET IMPLEMENTED!')
  orbit%state = lost$
  if (global_com%exit_on_error) call err_exit
  return
endif


if (ele%value(bragg_angle_in$) == 0) then
  call out_io (s_fatal$, r_name, 'REFERENCE ENERGY TOO SMALL TO SATISFY BRAGG CONDITION! FOR ELEMENT: ' // ele%name)
  orbit%state = lost_pz_aperture$
  if (global_com%exit_on_error) call err_exit
  return
endif

if (ele%value(thickness$) == 0) then 
  call out_io (s_error$, r_name, 'LAUE CRYSTAL WITH ZERO THICKNESS WILL HAVE NO DIFFRACTION. FOR ELEMENT: ' // ele%name)
  return
endif

if (ele%value(mosaic_thickness$) == 0) then 
  call out_io (s_error$, r_name, 'A MOSAIC CRYSTAL MUST HAVE A FINITE MOSAIC_THICKNESS PARAMETER. FOR ELEMENT: ' // ele%name)
  orbit%state = lost_pz_aperture$
  if (global_com%exit_on_error) call err_exit
  return
endif

if (ele%value(mosaic_angle_rms_out_plane$) < 0) then 
  call out_io (s_error$, r_name, 'MOSAIC CRYSTAL WITH MOSAIC_ANGLE_RMS_OUT_PLANE PARAMETER NOT', &
                                 'SET (OR SET NEGATIVE). FOR ELEMENT:' // ele%name)
  orbit%state = lost_pz_aperture$
  if (global_com%exit_on_error) call err_exit
  return
endif

n_layer = max(1, nint(ele%value(thickness$) / ele%value(mosaic_thickness$)))
thick = ele%value(thickness$) / n_layer

! (px, py, pz) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

call track_to_surface (ele, orbit, param, w_surface)
if (orbit%state /= alive$) return

! Some init

field_init = sum(orbit%field)
cp%wavelength = c_light * h_planck / orbit%p0c
cp%cap_gamma = r_e * cp%wavelength**2 / (pi * ele%value(v_unitcell$)) 
vec_init = orbit%vec
h_bar0 = ele%photon%material%h_norm * cp%wavelength / ele%value(d_spacing$)
p_factor = cos(ele%value(bragg_angle_in$) + ele%value(bragg_angle_out$))

! Decide how many diffractions there are

n_choose = (n_layer+1)/2
allocate (m_layer(n_layer), indx(n_layer), follow_prob(n_choose))

do im = 1, n_layer
  ml => m_layer(im)
  call ran_gauss(rr2)
  ml%theta_in  = rr2(1) * ele%value(mosaic_angle_rms_in_plane$)
  ml%theta_out = rr2(2) * ele%value(mosaic_angle_rms_out_plane$)

  orbit%vec = vec_init
  call rotate_vec(orbit%vec(2:6:2), y_axis$, ml%theta_in)
  call rotate_vec(orbit%vec(2:6:2), z_axis$, ml%theta_out)
  if (ele%photon%grid%type == h_misalign$) call crystal_h_misalign (ele, orbit, h_bar0) 

  cp%old_vvec = orbit%vec(2:6:2)
  cp%new_vvec = orbit%vec(2:6:2) + h_bar0
  cp%new_vvec(3) = sqrt(1 - cp%new_vvec(1)**2 - cp%new_vvec(2)**2)

  ! Calculate some parameters

  gamma_0 = cp%old_vvec(3)
  gamma_h = cp%new_vvec(3)

  cp%b_eff = gamma_0 / gamma_h
  cp%dtheta_sin_2theta = -dot_product(h_bar0 + 2 * cp%old_vvec, h_bar0) / 2

  ! Diffraction calc...

  call crystal_diffraction_field_calc (cp, ele, thick, param, p_factor, .true.,  field1, phase1, orbit%state, dr1)
  call crystal_diffraction_field_calc (cp, ele, thick, param, 1.0_rp,   .false., field2, phase2, orbit%state, dr2)   ! Sigma polarity

  ! Undiffrated calc...

  call crystal_diffraction_field_calc (cp, ele, thick, param, p_factor, .true.,  field3, phase3, orbit%state, dr3, follow_undiffracted = .true.)
  call crystal_diffraction_field_calc (cp, ele, thick, param, 1.0_rp,   .false., field4, phase4, orbit%state, dr4, follow_undiffracted = .true.)

  ml%branch_ratio = (field1**2 + field2**2) / (field1**2 + field2**2 + field3**2 + field4**2)
enddo

call indexer(m_layer%branch_ratio, indx)
indx = indx(n_layer:1:-1)
follow_prob(1) = 1
p_sum = 1
prob = 1

do im = 2, n_choose
  prob = prob * m_layer(indx(2*im-2))%branch_ratio * m_layer(indx(2*im-1))%branch_ratio
  p_sum = p_sum + prob
  follow_prob(im) = prob
enddo
follow_prob = follow_prob / p_sum

!

n = nint(ele%value(mosaic_diffraction_num$))
call ran_uniform(rr)
m_layer(indx(1))%follow_diffracted = .true.
p_sum = follow_prob(1)
do im = 2, n_choose
  if (n == 0 .and. p_sum > rr) exit
  if (n /= 0 .and. 2*im-2 > n) exit
  m_layer(indx(2*im-2))%follow_diffracted = .true.
  m_layer(indx(2*im-1))%follow_diffracted = .true.
  p_sum = p_sum + follow_prob(im)  
enddo

! And track through the mosaic layers

h_dir = 1
orbit%vec = vec_init

do im = 1, n_layer
  ml => m_layer(im)
  h_bar_here = h_dir * h_bar0

  ! Rotate to mosaic coords

  call rotate_vec(orbit%vec(2:6:2), y_axis$, ml%theta_in)
  call rotate_vec(orbit%vec(2:6:2), z_axis$, ml%theta_out)

  ! cp%new_vvec is the normalized outgoing wavevector outside the crystal for diffracted phtons.

  cp%old_vvec = orbit%vec(2:6:2)
  cp%new_vvec = orbit%vec(2:6:2) + h_bar_here
  cp%new_vvec(3) = sqrt(1 - cp%new_vvec(1)**2 - cp%new_vvec(2)**2)

  ! Calculate some parameters

  gamma_0 = cp%old_vvec(3)
  gamma_h = cp%new_vvec(3)

  cp%b_eff = gamma_0 / gamma_h
  cp%dtheta_sin_2theta = -dot_product(h_bar_here + 2 * cp%old_vvec, h_bar_here) / 2


  if (ml%follow_diffracted) then  ! diffracted
    call crystal_diffraction_field_calc (cp, ele, thick, param, p_factor, .true.,  field1, phase1, orbit%state, dr1)
    call crystal_diffraction_field_calc (cp, ele, thick, param, 1.0_rp,   .false., field2, phase2, orbit%state, dr2)   ! Sigma polarity
    h_dir = -h_dir
    orbit%vec(1:5:2) = orbit%vec(1:5:2) + (dr1 * field1**2 + dr2 * field2**2) / (field1**2 + field2**2)
    orbit%vec(2:6:2) = cp%new_vvec
    
  else
    call crystal_diffraction_field_calc (cp, ele, thick, param, p_factor, .true.,  field1, phase1, orbit%state, dr1, follow_undiffracted = .true.)
    call crystal_diffraction_field_calc (cp, ele, thick, param, 1.0_rp,   .false., field2, phase2, orbit%state, dr2, follow_undiffracted = .true.)
  endif

  orbit%field(1) = orbit%field(1) * field1
  orbit%field(2) = orbit%field(2) * field2
  orbit%phase    = orbit%phase + [phase1, phase2]

  ! Rotate back out from mosaic coords

  call rotate_vec(orbit%vec(2:6:2), z_axis$, -ml%theta_out)
  call rotate_vec(orbit%vec(2:6:2), y_axis$, -ml%theta_in)

enddo

if (has_curvature(ele%photon)) call rotate_for_curved_surface (ele, orbit, unset$, w_surface)

end subroutine track1_mosaic_crystal

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track1_crystal (ele, param, orbit)
!
! Routine to track diffraction from a crystal.
!
! Input:
!   ele      -- ele_struct: Element tracking through.
!   param    -- lat_param_struct: lattice parameters.
!   orbit    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   orbit    -- Coord_struct: final phase-space coords
!-

subroutine track1_crystal (ele, param, orbit)

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param
type (crystal_param_struct) cp
type (photon_reflect_table_struct), pointer :: rt

real(rp) h_norm(3), h_bar(3), e_tot, pc, p_factor, field(2), w_surface(3,3)
real(rp) gamma_0, gamma_h, dr1(3), dr2(3), dphase(2), theta, angle(4), energy(4), prob(4), coef(3)
real(rp) bragg_ang, p0, lambda, beta

integer ie, ia

character(*), parameter :: r_name = 'track1_cyrstal'

! A graze angle of zero means the wavelength of the reference photon was too large
! for the bragg condition. 

if (orbit_too_large(orbit, param)) return

if (ele%value(bragg_angle_in$) == 0) then
  call out_io (s_fatal$, r_name, 'REFERENCE ENERGY TOO SMALL TO SATISFY BRAGG CONDITION!')
  orbit%state = lost_pz_aperture$
  if (global_com%exit_on_error) call err_exit
  return
endif

if (ele%value(b_param$) > 0 .and. ele%value(thickness$) == 0) then 
  call out_io (s_error$, r_name, 'LAUE CRYSTAL WITH ZERO THICKNESS WILL HAVE NO DIFFRACTION: ' // ele%name)
endif

!

cp%wavelength = c_light * h_planck / orbit%p0c
cp%cap_gamma = r_e * cp%wavelength**2 / (pi * ele%value(v_unitcell$)) 

! (px, py, pz) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

call track_to_surface (ele, orbit, param, w_surface)
if (orbit%state /= alive$) return

! Construct h_bar = H * wavelength.

h_norm = ele%photon%material%h_norm
h_bar = h_norm
if (ele%photon%grid%type == h_misalign$) call crystal_h_misalign (ele, orbit, h_bar) 
h_bar = h_bar * cp%wavelength / ele%value(d_spacing$)

! cp%new_vvec is the normalized outgoing wavevector outside the crystal

cp%old_vvec = orbit%vec(2:6:2)
cp%new_vvec = orbit%vec(2:6:2) + h_bar

if (ele%value(b_param$) < 0) then ! Bragg
  cp%new_vvec(3) = -sqrt(1 - cp%new_vvec(1)**2 - cp%new_vvec(2)**2)
else
  cp%new_vvec(3) = sqrt(1 - cp%new_vvec(1)**2 - cp%new_vvec(2)**2)
endif

! Calculate some parameters

gamma_0 = cp%old_vvec(3)
gamma_h = cp%new_vvec(3)

cp%b_eff = gamma_0 / gamma_h
cp%dtheta_sin_2theta = -dot_product(h_bar + 2 * cp%old_vvec, h_bar) / 2

! E field calc

if (is_true(ele%value(use_reflectivity_table$))) then
  rt => ele%photon%reflectivity_table_sigma
  theta = atan2(cp%old_vvec(3), norm2(cp%old_vvec(1:2)))

  lambda = c_light * h_planck / orbit%p0c
  beta = lambda / (2 * ele%value(d_spacing$))
  if (ele%value(b_param$) < 0) then ! Bragg
    bragg_ang = asin((-beta * h_norm(3) - h_norm(1) * sqrt(h_norm(1)**2 + h_norm(3)**2 - beta**2)) / (h_norm(1)**2 + h_norm(3)**2))
  else
    bragg_ang = asin((-beta * h_norm(1) + h_norm(3) * sqrt(h_norm(1)**2 + h_norm(3)**2 - beta**2)) / (h_norm(1)**2 + h_norm(3)**2))
  endif

  ie = bracket_index(orbit%p0c, rt%energy, 1)
  ia = bracket_index(theta-bragg_ang, rt%angle, 1)

  if (ie < 1 .or. ie >= size(rt%energy) .or. ia < 1 .or. ia >= size(rt%angle)) then
    orbit%field = 0

  else
    energy = [rt%energy(ie), rt%energy(ie), rt%energy(ie+1), rt%energy(ie+1)]
    angle = [rt%angle(ia), rt%angle(ia+1), rt%angle(ia), rt%angle(ia+1)]
    prob = [rt%p_reflect(ia, ie), rt%p_reflect(ia+1, ie), rt%p_reflect(ia, ie+1), rt%p_reflect(ia+1, ie+1)]
    call linear_fit_2D(energy, angle, prob, coef)
    p0 = coef(1) * orbit%p0c + coef(2) * (theta-bragg_ang) + coef(3)
    orbit%field = orbit%field * sqrt(p0)
  endif

else
  p_factor = cos(ele%value(bragg_angle_in$) + ele%value(bragg_angle_out$))
  call crystal_diffraction_field_calc (cp, ele, ele%value(thickness$), param, p_factor, .true.,  field(1), dphase(1), orbit%state, dr1)
  call crystal_diffraction_field_calc (cp, ele, ele%value(thickness$), param, 1.0_rp,   .false., field(2), dphase(2), orbit%state, dr2)   ! Sigma pol

  orbit%field = orbit%field * field
  orbit%phase = orbit%phase + dphase

  ! For Laue: Average trajectories for the two polarizations weighted by the fields.
  ! This approximation is valid as long as the two trajectories are "close" enough.

  if (ele%value(b_param$) > 0 .and. (field(1) /= 0 .or. field(2) /= 0)) then ! Laue
    orbit%vec(1:5:2) = orbit%vec(1:5:2) + (dr1 * field(1) + dr2 * field(2)) / (field(1) + field(2))
  endif
endif

! Rotate back from curved body coords to element coords

if (ele%value(b_param$) < 0) then ! Bragg
  orbit%vec(2:6:2) = cp%new_vvec
else
  ! forward_diffracted and undiffracted beams do not have an orientation change.
  if (nint(ele%value(ref_orbit_follows$)) == bragg_diffracted$) orbit%vec(2:6:2) = cp%new_vvec
endif

if (has_curvature(ele%photon)) call rotate_for_curved_surface (ele, orbit, unset$, w_surface)

end subroutine track1_crystal

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine crystal_h_misalign (ele, orbit, h_vec)
!
! Routine reorient the crystal H vector due to local imperfections in the crystal lattice.
!
! Input:
!   ele      -- ele_struct: Crystal element
!   orbit    -- coord_struct: Photon position at crystal surface.
!   h_vec(3) -- real(rp): H vector before misalignment.
!
! Output:
!   h_vec(3) -- real(rp): H vector after misalignment.
!-

subroutine crystal_h_misalign (ele, orbit, h_vec)

type (ele_struct), target :: ele
type (coord_struct) orbit
type (photon_element_struct), pointer :: ph
type (surface_orientation_struct), pointer :: orient

real(rp) h_vec(3), r(2)
integer ij(2)
character(*), parameter :: r_name = 'crystal_h_misalign'

!

ph => ele%photon
if (.not. ph%grid%active) return

ij = nint((orbit%vec(1:3:2) - ph%grid%r0) / ph%grid%dr)

if (any(ij < lbound(ph%grid%pt)) .or. any(ij > ubound(ph%grid%pt))) then
  call out_io (s_error$, r_name, &
              'Photon position on crystal surface outside of grid bounds for element: ' // ele%name)
  return
endif

! Make small angle approximation

orient => ph%grid%pt(ij(1), ij(2))%orientation

h_vec(1:2) = h_vec(1:2) + [orient%dz_dx, orient%dz_dy]
if (orient%dz_dx_rms /= 0 .or. orient%dz_dy /= 0) then
  call ran_gauss (r)
  h_vec(1:2) = h_vec(1:2) + [orient%dz_dx_rms, orient%dz_dy_rms] * r
endif

h_vec(3) = sqrt(1 - h_vec(1)**2 - h_vec(2)**2)

end subroutine crystal_h_misalign

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine target_rot_mats (r_center, w_to_target, w_to_ele)
!
! Routine to calculate the rotation matrices between ele coords and "target" coords.
! By definition, in target coords r_center = [0, 0, 1].
!
! Input:
!   r_center(3)   -- real(rp): In lab coords: Center of target relative to phton emission point.
!
! Output:
!   w_to_target(3,3) -- real(rp): Rotation matrix from ele to target coords.
!   w_to_ele(3,3)    -- real(rp): Rotation matrix from target to ele coords.
!-

subroutine target_rot_mats (r_center, w_to_target, w_to_ele)

real(rp) r_center(3), w_to_target(3,3), w_to_ele(3,3)
real(rp) r(3), cos_theta, sin_theta, cos_phi, sin_phi

!

r = r_center / norm2(r_center)
sin_phi = r(2)
cos_phi = sqrt(r(1)**2 + r(3)**2)
if (cos_phi == 0) then
  sin_theta = 0  ! Arbitrary
  cos_theta = 1
else
  sin_theta = r(1) / cos_phi
  cos_theta = r(3) / cos_phi
endif

w_to_ele(1,:) = [ cos_theta, -sin_theta * sin_phi, sin_theta * cos_phi]
w_to_ele(2,:) = [ 0.0_rp,     cos_phi,             sin_phi]
w_to_ele(3,:) = [-sin_theta, -cos_theta * sin_phi, cos_theta * cos_phi]

w_to_target(1,:) = [ cos_theta,           0.0_rp,  -sin_theta]
w_to_target(2,:) = [-sin_theta * sin_phi, cos_phi, -cos_theta * sin_phi]
w_to_target(3,:) = [ sin_theta * cos_phi, sin_phi,  cos_theta * cos_phi]

end subroutine target_rot_mats

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine target_min_max_calc (r_corner1, r_corner2, y_min, y_max, phi_min, phi_max, initial)
!
! Routine to calculate the min/max values for (y, phi).
! min/max values are cumulative.
!
! Input:
!   r_corner1(3)     -- real(rp): In target coords: A corner of the target. Must be normalized to 1.
!   r_corner2(3)     -- real(rp): In target coords: Adjacent corner of the target. Must be normalized to 1.
!   y_min, y_max     -- real(rp): min/max values. Only needed if initial = False.
!   phi_min, phi_max -- real(rp): min/max values. Only needed if initial = False.
!   initial          -- logical, optional: If present and True then this is the first edge for computation.
!
! Output:
!   y_min, y_max     -- real(rp): min/max values. 
!   phi_min, phi_max -- real(rp): min/max values. 
!-

subroutine target_min_max_calc (r_corner1, r_corner2, y_min, y_max, phi_min, phi_max, initial)

real(rp)  r_corner1(3), r_corner2(3), y_min, y_max, phi_min, phi_max
real(rp) phi1, phi2, k, t, alpha, beta, y

logical, optional :: initial

!

phi1 = atan2(r_corner1(1), r_corner1(3))
phi2 = atan2(r_corner2(1), r_corner2(3))

if (logic_option(.false., initial)) then
  y_max = max(r_corner1(2), r_corner2(2))
  y_min = min(r_corner1(2), r_corner2(2))
  phi_max = max(phi1, phi2)
  phi_min = min(phi1, phi2)
else
  y_max = max(y_max, r_corner1(2), r_corner2(2))
  y_min = min(y_min, r_corner1(2), r_corner2(2))
  phi_max = max(phi_max, phi1, phi2)
  phi_min = min(phi_min, phi1, phi2)
endif

k = dot_product(r_corner1, r_corner2)
t = abs(k * r_corner1(2) - r_corner2(2)) 
alpha = t / sqrt((1-k**2)*(1-k**2 + t**2))

if (alpha < 1) then
  beta = sqrt((alpha * k)**2 + 1 - alpha**2) - alpha * k
  y = r_corner1(2) * alpha + r_corner2(2) * beta
  y_max = max(y_max, y)
  y_min = min(y_min, y)
endif

end subroutine target_min_max_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_bend_photon (orb, ele, length)
!
! Routine to track a photon through a dipole bend.
! The photon is traveling in a straight line but the reference frame
! is curved in a circular shape.
!
! Input:
!   orb        -- Coord_struct: Starting position.
!   ele        -- Ele_struct: Bend element.
!   length     -- real(rp): length to track.
!
! Output:
!   orb         -- Coord_struct: End position.
!-

subroutine track_a_bend_photon (orb, ele, length)

type (coord_struct) orb
type (ele_struct) ele

real(rp) length
real(rp) g, radius, theta, tan_t, dl, st, ct, denom, sin_t, cos_t
real(rp) v_x, v_s

!

g = ele%value(g$)

! g = 0 case

if (g == 0) then
  call track_a_drift_photon (orb, length, .true.)
  return
endif

! Normal case

if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords_photon (ele%value(ref_tilt_tot$), orb%vec)

radius = 1 / g
theta = length * g
tan_t = tan(theta)
dl = tan_t * (radius + orb%vec(1)) / (orb%vec(6) - tan_t * orb%vec(2))

! Move to the stop point. 
! Need to remember that radius can be negative.

st = dl * orb%vec(6)
ct = radius + orb%vec(1) + dl * orb%vec(2)
if (abs(st) < 1d-3 * ct) then
  denom = sign (ct * (1 + (st/ct)**2/2 + (st/ct)**4/8), radius)
else
  denom = sign (sqrt((radius + orb%vec(1) + dl * orb%vec(2))**2 + (dl * orb%vec(6))**2), radius)
endif
sin_t = st / denom
cos_t = ct / denom
v_x = orb%vec(2); v_s = orb%vec(6)

orb%vec(1) = denom - radius
orb%vec(2) = v_s * sin_t + v_x * cos_t
orb%vec(3) = orb%vec(3) + dl * orb%vec(4)
orb%vec(5) = orb%vec(5) + length
orb%vec(6) = v_s * cos_t - v_x * sin_t
orb%s      = orb%s + length
orb%t      = orb%t + length * orb%vec(6) / c_light

if (ele%value(ref_tilt_tot$) /= 0) call tilt_coords_photon (-ele%value(ref_tilt_tot$), orb%vec)

end subroutine track_a_bend_photon

end module
