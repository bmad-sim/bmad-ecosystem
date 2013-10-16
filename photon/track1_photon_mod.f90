module track1_photon_mod

use track1_mod
use wall3d_mod

! This is for passing info into the field_calc_routine

type, private :: crystal_param_struct
  real(rp) cap_gamma, dtheta_sin_2theta, b_eff
  complex(rp) f0, fh
end type

private e_field_calc

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine diffraction_plate_hit_spot (ele, orbit, is_clear, ix_section)
!
! Routine to determine where a photon hits on a diffraction_plate element.
!
! Note: It is assumed that orbit is in the frame of reference of the element.
! That is, offset_photon needs to be called before this routine.
!
! Input:
!   ele     -- ele_struct: diffraction_plate element.
!   orbit   -- coord_struct: photon position.
!
! Output:
!   is_clear   -- logical, optional: True if hitting clear section. False otherwise.
!   ix_section -- integer, optional: Set to index of section hit. Set to zero if
!                   photon is outside all sections.
!-

subroutine diffraction_plate_hit_spot (ele, orbit, is_clear, ix_section)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (wall3d_section_struct), pointer :: sec

integer, optional :: ix_section
logical, optional :: is_clear

integer i, ix_sec
logical opening_found, in_this_aperture, cleared

! Logic: A particle is in a clear section if it is inside the section and outside
! all subsiquent opaque sections up to the next clear section.

in_this_aperture = .false.
cleared = .false.
ix_sec = 0

do i = 1, size(ele%wall3d%section)
  sec => ele%wall3d%section(i)

  if (sec%type == clear$) then
    if (cleared) exit
    if (.not. in_section(sec)) then
      in_this_aperture = .false.
      cycle
    endif
    in_this_aperture = .true.
    ix_sec = i

  else    ! opaque
    if (.not. in_this_aperture) cycle
    if (.not. in_section(sec)) cycle
    in_this_aperture = .false.
    cleared = .false.
    ix_sec = i

  endif

enddo
  
if (present(ix_section)) ix_section = ix_sec
if (present(is_clear)) is_clear = cleared

!------------------------------------------------------------
contains

function in_section(sec) result (is_in)

type (wall3d_section_struct) sec
logical is_in

real(rp) x, y, norm, r_wall, dr_dtheta

!

x = orbit%vec(1) - sec%x0;  y = orbit%vec(3) - sec%y0

if (x == 0 .and. y == 0) then
  is_in = .true.
  return
endif

norm = norm2([x, y])
call calc_wall_radius (sec%v, x/norm, y/norm, r_wall, dr_dtheta)
is_in = (norm <= r_wall)

end function in_section

end subroutine diffraction_plate_hit_spot

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

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param
type (photon_target_struct), pointer :: target
type (target_point_struct) corner(4)

real(rp) wavelength, ran(2), r_particle(3), w_to_target(3,3), w_to_ele(3,3), w_to_surface(3,3)
real(rp) phi_min, phi_max, y_min, y_max, y, phi, rho, r(3)
real(rp), pointer :: val(:)

integer n, i, ix

character(*), parameter :: r_name = 'track1_sample'

!

val => ele%value
wavelength = c_light * h_planck / orbit%p0c

call track_to_surface (ele, orbit, curved_surface_rot = .false.)
if (orbit%state /= alive$) return

if (ele%photon%surface%has_curvature) then
  call rotate_for_curved_surface (ele, orbit, set$, w_to_surface)
endif

! Check aperture

if (ele%aperture_at == surface$) then
  call check_aperture_limit (orbit, ele, surface$, param)
  if (orbit%state /= alive$) return
endif

! Reflect 

select case (ele%photon%surface%type)
case (isotropic_emission$)

  call ran_uniform(ran)

  target => ele%photon%target
  if (target%enabled) then
    r_particle = orbit%vec(1:5:2)
    r = target%center%r - r_particle
    if (ele%photon%surface%has_curvature) r = matmul(w_to_surface, r)
    call target_rot_mats (r, w_to_target, w_to_ele)
    do i = 1, 4
      corner(i)%r = matmul(w_to_target, target%corner(i)%r - r_particle)
      if (ele%photon%surface%has_curvature) corner(i)%r = matmul(w_to_surface, corner(i)%r)
      corner(i)%r = corner(i)%r / norm2(corner(i)%r)
    enddo
    call target_min_max_calc (corner(4)%r, corner(1)%r, y_min, y_max, phi_min, phi_max, .true.)
    call target_min_max_calc (corner(1)%r, corner(2)%r, y_min, y_max, phi_min, phi_max)
    call target_min_max_calc (corner(2)%r, corner(3)%r, y_min, y_max, phi_min, phi_max)
    call target_min_max_calc (corner(3)%r, corner(4)%r, y_min, y_max, phi_min, phi_max)
    y = y_min + (y_max-y_min) * ran(1)
    phi = phi_min + (phi_max-phi_min) * ran(2)
    rho = sqrt(1 - y*y)
    orbit%vec(2:6:2) = [rho * sin(phi), y, rho * cos(phi)]
    orbit%vec(2:6:2) = matmul(w_to_ele, orbit%vec(2:6:2))
    if (param%tracking_type == coherent$) then
      orbit%field = orbit%field * (y_max - y_min) * (phi_max - phi_min) / fourpi
    else
      orbit%field = orbit%field * sqrt ((y_max - y_min) * (phi_max - phi_min) / fourpi)
    endif

  else
    y = 2 * ran(1) - 1
    phi = pi * (ran(2) - 0.5_rp)
    rho = sqrt(1 - y**2)
    orbit%vec(2:6:2) = [rho * sin(phi), y, -rho * cos(phi)]
    if (param%tracking_type == coherent$) then
      orbit%field = orbit%field / 2       ! Half the photons get lost by being emitted into the bulk.
    else
      orbit%field = orbit%field / sqrt_2  ! Half the photons get lost by being emitted into the bulk.
    endif
  endif

case default
  call out_io (s_error$, r_name, 'SURFACE TYPE NOT SET.')
  orbit%state = lost$
end select

! Rotate back to uncurved element coords

if (ele%photon%surface%has_curvature) then
  call rotate_for_curved_surface (ele, orbit, unset$)
endif

end subroutine track1_sample

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

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

real(rp) wavelength
real(rp), pointer :: val(:)

character(*), parameter :: r_name = 'track1_mirror'

!

val => ele%value
wavelength = c_light * h_planck / orbit%p0c

call track_to_surface (ele, orbit)
if (orbit%state /= alive$) return

! Check aperture

if (ele%aperture_at == surface$) then
  call check_aperture_limit (orbit, ele, surface$, param)
  if (orbit%state /= alive$) return
endif

! Reflect momentum vector

orbit%vec(6) = -orbit%vec(6)

! Rotate back to uncurved element coords

if (ele%photon%surface%has_curvature) then
  call rotate_for_curved_surface (ele, orbit, unset$)
endif

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

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param

real(rp) wavelength, kz_air
real(rp), pointer :: val(:)

complex(rp) zero, xi_1, xi_2, kz1, kz2, c1, c2

character(*), parameter :: r_name = 'track1_multilayer_mirror'

!

val => ele%value
wavelength = c_light * h_planck / orbit%p0c

call track_to_surface (ele, orbit)
if (orbit%state /= alive$) return

! Check aperture

if (ele%aperture_at == surface$) then
  call check_aperture_limit (orbit, ele, surface$, param)
  if (orbit%state /= alive$) return
endif

! Note: Koln z-axis = Bmad x-axis.
! Note: f0_re and f0_im are both positive.

xi_1 = cmplx(-val(f0_re1$), val(f0_im1$)) * r_e * wavelength**2 / (pi * val(v1_unitcell$)) 
xi_2 = cmplx(-val(f0_re2$), val(f0_im2$)) * r_e * wavelength**2 / (pi * val(v2_unitcell$)) 

kz1 = twopi * sqrt(orbit%vec(6)**2 + xi_1) / wavelength
kz2 = twopi * sqrt(orbit%vec(6)**2 + xi_2) / wavelength
kz_air = twopi * orbit%vec(6) / wavelength

c1 = exp(I_imaginary * kz1 * val(d1_thickness$) / 2)
c2 = exp(I_imaginary * kz2 * val(d2_thickness$) / 2)

zero = cmplx(0.0_rp, 0.0_rp)

call multilayer_track (xi_1, xi_2, orbit%field(1), orbit%phase(1))     ! pi polarization
call multilayer_track (zero, zero, orbit%field(2), orbit%phase(2))     ! sigma polarization

! Reflect momentum vector

orbit%vec(6) = -orbit%vec(6)

! Rotate back to uncurved element coords

if (ele%photon%surface%has_curvature) then
  call rotate_for_curved_surface (ele, orbit, unset$)
endif

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
! Subroutine track1_crystal (ele, param, orbit)
!
! Routine to track reflection from a crystal.
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

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: orbit
type (lat_param_struct) :: param
type (crystal_param_struct) c_param

real(rp) h_bar(3), e_tot, pc, p_factor, wavelength
real(rp) gamma_0, gamma_h, old_vec(6)

character(*), parameter :: r_name = 'track1_cyrstal'

! A graze angle of zero means the wavelength of the reference photon was too large
! for the bragg condition. 

if (ele%value(bragg_angle_in$) == 0) then
  call out_io (s_fatal$, r_name, 'REFERENCE ENERGY TOO SMALL TO SATISFY BRAGG CONDITION!')
  orbit%state = lost_z_aperture$
  if (global_com%exit_on_error) call err_exit
  return
endif

!

wavelength = c_light * h_planck / orbit%p0c

! (px, py, pz) coords are with respect to laboratory reference trajectory.
! Convert this vector to k0_outside_norm which are coords with respect to crystal surface.
! k0_outside_norm is normalized to 1.

call track_to_surface (ele, orbit)
if (orbit%state /= alive$) return

old_vec = orbit%vec

! Check aperture

if (ele%aperture_at == surface$) then
  call check_aperture_limit (orbit, ele, surface$, param)
  if (orbit%state /= alive$) return
endif

! Construct h_bar = H * wavelength.


h_bar = [ele%value(h_x_norm$), ele%value(h_y_norm$), ele%value(h_z_norm$)] 
if (ele%photon%surface%grid%type == h_misalign$) call crystal_h_misalign (ele, orbit, h_bar) 
h_bar = h_bar * wavelength / ele%value(d_spacing$)

! kh_outside_norm is the normalized outgoing wavevector outside the crystal

c_param%cap_gamma = r_e * wavelength**2 / (pi * ele%value(v_unitcell$)) 

orbit%vec(2:6:2) = orbit%vec(2:6:2) + h_bar

if (ele%value(b_param$) < 0) then ! Bragg
  orbit%vec(6) = -sqrt(1 - orbit%vec(2)**2 - orbit%vec(4)**2)
else
  orbit%vec(6) = sqrt(1 - orbit%vec(2)**2 - orbit%vec(4)**2)
endif

gamma_0 = old_vec(6)
gamma_h = orbit%vec(6)

!-------------------------------------
! Calculate phase and intensity

c_param%b_eff             = gamma_0 / gamma_h
c_param%dtheta_sin_2theta = -dot_product(h_bar + 2 * old_vec(2:6:2), h_bar) / 2
c_param%f0                = cmplx(ele%value(f0_re$), ele%value(f0_im$)) 
c_param%fh                = cmplx(ele%value(fh_re$), ele%value(fh_im$))

p_factor = cos(ele%value(bragg_angle_in$) + ele%value(bragg_angle_out$))
call e_field_calc (c_param, ele, p_factor, orbit%field(1), orbit%phase(1))
call e_field_calc (c_param, ele, 1.0_rp,   orbit%field(2), orbit%phase(2))   ! Sigma polarization

! Rotate back from curved body coords to element coords

if (ele%photon%surface%has_curvature) then
  call rotate_for_curved_surface (ele, orbit, unset$)
endif

end subroutine track1_crystal

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine e_field_calc (cp, ele, p_factor, e_field, e_phase)
!
! Routine to compute position where crystal reflects in crystal coordinates.
!
! Input:
!   cp       -- crystal_param_struct: Crystal parameters.
!   ele      -- ele_struct: Crystal element.
!   p_factor -- Real(rp)
!
! Output:
!   e_field -- Real(rp)
!   e_phase -- Real(rp)
!-

subroutine e_field_calc (cp, ele, p_factor, e_field, e_phase)

type (crystal_param_struct) cp
type (ele_struct) ele

real(rp) p_factor, e_field, e_phase, sqrt_b

complex(rp) e_rel, e_rel_a, e_rel_b, eta, eta1, f_cmp, xi_0k_a, xi_hk_a, xi_0k_b, xi_hk_b

! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

sqrt_b = sqrt(abs(cp%b_eff))

eta = (cp%b_eff * cp%dtheta_sin_2theta + cp%f0 * cp%cap_gamma * (1.0_rp - cp%b_eff)/2) / &
                                              (cp%cap_gamma * abs(p_factor) * sqrt_b * cp%fh) 
eta1 = sqrt(eta**2 + sign(1.0_rp, cp%b_eff))
f_cmp = abs(p_factor) * sqrt_b * cp%cap_gamma * cp%fh / 2

xi_0k_b = f_cmp * (eta - eta1)
xi_hk_b = f_cmp / (abs(cp%b_eff) * (eta - eta1))

xi_0k_a = f_cmp * (eta + eta1)
xi_hk_a = f_cmp / (abs(cp%b_eff) * (eta + eta1))

! Bragg

if (ele%value(b_param$) < 0) then 
  if (abs(eta+eta1) > abs(eta-eta1)) then
    e_rel = -2.0_rp * xi_0k_b / (p_factor * cp%cap_gamma * cp%fh)
  else
    e_rel = -2.0_rp * xi_0k_a / (p_factor * cp%cap_gamma * cp%fh)
  endif

  ! Factor of sqrt_b comes from geometrical change in the transverse width of the photon beam

  e_field = e_field * abs(e_rel) / sqrt_b
  e_phase = atan2(aimag(e_rel), real(e_rel)) + e_phase

!---------------
! Laue calc

else 

  e_rel_a = -2.0_rp * xi_0k_a / (p_factor * cp%cap_gamma * cp%fh)
  e_rel_b = -2.0_rp * xi_0k_b / (p_factor * cp%cap_gamma * cp%fh)

endif

end subroutine e_field_calc

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine track_to_surface (ele, orbit, curved_surface_rot)
!
! Routine to track a photon to the surface of the element.
!
! If the surface is curved, the photon's velocity coordinates are rotated so that 
! orbit%vec(6) is maintained to be normal to the local surface and pointed inward.
!
! Notice that, to save time, the photon's position is not rotated when there is a curved surface.
!
! Input:
!   ele                 -- ele_struct: Element
!   orbit               -- coord_struct: Coordinates in the element coordinate frame
!   curved_surface_rot  -- Logical, optional, If present and False then do not rotate velocity coords.
!
! Output:
!   orbit      -- coord_struct: At surface in local surface coordinate frame
!   err        -- logical: Set true if surface intersection cannot be found. 
!-

subroutine track_to_surface (ele, orbit, curved_surface_rot)

use nr, only: zbrent

implicit none

type (ele_struct) ele
type (coord_struct) orbit
type (segmented_surface_struct), pointer :: segment

real(rp) :: s_len, s1, s2, s_center, x0, y0
character(*), parameter :: r_name = 'track_to_surface'
logical, optional :: curved_surface_rot

! If there is curvature, compute the reflection point which is where 
! the photon intersects the surface.

if (ele%photon%surface%has_curvature) then

  ele%photon%surface%segment%ix = int_garbage$; ele%photon%surface%segment%iy = int_garbage$

  ! Assume flat crystal, compute s required to hit the intersection

  s_center = orbit%vec(5) / orbit%vec(6)

  s1 = s_center
  s2 = s_center
  if (photon_depth_in_crystal(s_center) > 0) then
    do
      s1 = s1 - 0.1
      if (photon_depth_in_crystal(s1) < 0) exit
      if (s1 < -10) then
        !! call out_io (s_warn$, r_name, &
        !!      'PHOTON INTERSECTION WITH SURFACE NOT FOUND FOR ELEMENT: ' // ele%name)
        orbit%state = lost$
        return
      endif
    enddo
  else
    do
      s2 = s2 + 0.1
      if (photon_depth_in_crystal(s2) > 0) exit
      if (s1 > 10) then
        !! call out_io (s_warn$, r_name, &
        !!      'PHOTON INTERSECTION WITH SURFACE NOT FOUND FOR ELEMENT: ' // ele%name)
        orbit%state = lost$
        return
      endif
    enddo
  endif

  s_len = zbrent (photon_depth_in_crystal, s1, s2, 1d-10)

  ! Compute the intersection point

  orbit%vec(1:5:2) = s_len * orbit%vec(2:6:2) + orbit%vec(1:5:2)
  orbit%t = orbit%t + s_len / c_light

  if (logic_option(.true., curved_surface_rot)) call rotate_for_curved_surface (ele, orbit, set$)

else
  s_len = -orbit%vec(5) / orbit%vec(6)
  orbit%vec(1:5:2) = orbit%vec(1:5:2) + s_len * orbit%vec(2:6:2) ! Surface is at z = 0
  orbit%t = orbit%t + s_len / c_light
endif

contains

!-----------------------------------------------------------------------------------------------
!+
! Function photon_depth_in_crystal (s_len) result (delta_h)
! 
! Private routine to be used as an argument in zbrent. Propagates
! photon forward by a distance s_len. Returns delta_h = z-z0
! where z0 is the height of the crystal surface. 
! Since positive z points inward, positive delta_h => inside crystal.
!
! Input:
!   s_len   -- Real(rp): Place to position the photon.
!
! Output:
!   delta_h -- Real(rp): Depth of photon below surface in crystal coordinates.
!-

function photon_depth_in_crystal (s_len) result (delta_h)

implicit none

type (photon_surface_struct), pointer :: surf

real(rp), intent(in) :: s_len
real(rp) :: delta_h
real(rp) :: point(3), x, y
integer ix, iy, i_pt

!

point = s_len * orbit%vec(2:6:2) + orbit%vec(1:5:2)

surf => ele%photon%surface
x = point(1)
y = point(2)


if (surf%grid%type == segmented$) then
  call init_surface_segment (x, y, ele)

  delta_h = point(3) - surf%segment%z0 + (x - surf%segment%x0) * surf%segment%slope_x + &
                                         (y - surf%segment%y0) * surf%segment%slope_y

else
  delta_h = point(3)
  do ix = 0, ubound(surf%curvature_xy, 1)
  do iy = 0, ubound(surf%curvature_xy, 2) - ix
    if (ele%photon%surface%curvature_xy(ix, iy) == 0) cycle
    delta_h = delta_h + surf%curvature_xy(ix, iy) * x**ix * y**iy
  enddo
  enddo
endif

end function photon_depth_in_crystal

end subroutine track_to_surface

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine rotate_for_curved_surface (ele, orbit, set, rot_mat)
!
! Routine to rotate just the velocity coords between element body coords and effective 
! body coords ("curved body coords") with respect to the surface at the point of photon impact.
!
! Input:
!   ele      -- ele_struct: reflecting element
!   orbit    -- coord_struct: Photon position.
!   set      -- Logical: True -> Transform body to curved body. 
!                        False -> Transform curved body to body.
!
! Output:
!   orbit        -- coord_struct: Photon position.
!   rot_mat(3,3) -- real(rp), optional: rotation matrix applied to 
!-

subroutine rotate_for_curved_surface (ele, orbit, set, rot_mat)

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (photon_surface_struct), pointer :: s

real(rp), optional :: rot_mat(3,3)
real(rp) curve_rot(3,3), angle
real(rp) slope_y, slope_x, x, y
integer ix, iy

logical set

! Compute the slope of the crystal at that the point of impact.
! curve_rot transforms from standard body element coords to body element coords at point of impact.

s => ele%photon%surface
x = orbit%vec(1)
y = orbit%vec(3)

if (s%grid%type == segmented$) then
  call init_surface_segment (x, y, ele)
  slope_x = s%segment%slope_x; slope_y = s%segment%slope_y

else

  slope_x = 0
  slope_y = 0

  do ix = 0, ubound(s%curvature_xy, 1)
  do iy = 0, ubound(s%curvature_xy, 2) - ix
    if (s%curvature_xy(ix, iy) == 0) cycle
    if (ix > 0) slope_x = slope_x - ix * s%curvature_xy(ix, iy) * x**(ix-1) * y**iy
    if (iy > 0) slope_y = slope_y - iy * s%curvature_xy(ix, iy) * x**ix * y**(iy-1)
  enddo
  enddo

endif

if (slope_x == 0 .and. slope_y == 0) return

! Compute rotation matrix and goto body element coords at point of photon impact

angle = atan2(sqrt(slope_x**2 + slope_y**2), 1.0_rp)
if (set) angle = -angle
call axis_angle_to_w_mat ([slope_y, -slope_x, 0.0_rp], angle, curve_rot)

orbit%vec(2:6:2) = matmul(curve_rot, orbit%vec(2:6:2))
if (present(rot_mat)) rot_mat = curve_rot

end subroutine rotate_for_curved_surface

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine init_surface_segment (x, y, ele)
!
! Routine to init the componentes in ele%photon%surface%segment for use with segmented surface calculations.
! The segment used is determined by the (x,y) photon coordinates
!
! Input:
!   x, y   -- Real(rp): Coordinates of the photon.
!   ele    -- ele_struct: Elment containing a surface.
!
! Output:
!   ele    -- ele_struct: Element with ele%photon%surface%segment initialized.
!-

subroutine init_surface_segment (x, y, ele)

type (ele_struct), target :: ele
type (photon_surface_struct), pointer :: s
type (segmented_surface_struct), pointer :: seg

real(rp) x, y, x0, y0, dx, dy, coef_xx, coef_xy, coef_yy, coef_diag
integer ix, iy

! Only redo the cacluation if needed

s => ele%photon%surface
seg => s%segment

ix = nint(x / s%grid%dr(1))
iy = nint(y / s%grid%dr(2))

if (ix == seg%ix .and. iy == seg%iy) return

!

x0 = ix * s%grid%dr(1)
y0 = iy * s%grid%dr(2)

seg%ix = ix
seg%iy = iy

seg%x0 = x0
seg%y0 = y0
seg%z0 = 0

seg%slope_x = 0
seg%slope_y = 0
coef_xx = 0; coef_xy = 0; coef_yy = 0

do ix = 0, ubound(s%curvature_xy, 1)
do iy = 0, ubound(s%curvature_xy, 2) - ix
  if (s%curvature_xy(ix, iy) == 0) cycle
  seg%z0 = seg%z0 - s%curvature_xy(ix, iy) * x0**ix * y0**iy
  if (ix > 0) seg%slope_x = seg%slope_x - ix * s%curvature_xy(ix, iy) * x0**(ix-1) * y0**iy
  if (iy > 0) seg%slope_y = seg%slope_y - iy * s%curvature_xy(ix, iy) * x0**ix * y0**(iy-1)
  if (ix > 1) coef_xx = coef_xx - ix * (ix-1) * s%curvature_xy(ix, iy) * x0**(ix-2) * y0**iy / 2
  if (iy > 1) coef_yy = coef_yy - iy * (iy-1) * s%curvature_xy(ix, iy) * x0**ix * y0**(iy-2) / 2
  if (ix > 0 .and. iy > 0) coef_xy = coef_xy - ix * iy * s%curvature_xy(ix, iy) * x0**(ix-1) * y0**(iy-1)
enddo
enddo

! Correct for fact that segment is supported at the corners of the segment
! This correction only affects z0 and not the slopes

dx = s%grid%dr(1) / 2
dy = s%grid%dr(2) / 2
coef_xx = coef_xx * dx**2
coef_xy = coef_xy * dx * dy
coef_yy = coef_yy * dy**2
coef_diag = coef_xx + coef_yy - abs(coef_xy)

if (coef_diag < coef_xx .and. coef_diag < coef_yy) then
  seg%z0 = seg%z0 + coef_diag
else if (coef_xx < coef_yy) then
  seg%z0 = seg%z0 + coef_xx
else
  seg%z0 = seg%z0 + coef_yy
endif

end Subroutine init_surface_segment 

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

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (photon_surface_struct), pointer :: s
type (surface_grid_pt_struct), pointer :: pt

real(rp) h_vec(3), r(2)
integer ij(2)
character(*), parameter :: r_name = 'crystal_h_misalign'

!

s => ele%photon%surface

ij = nint((orbit%vec(1:3:2) + s%grid%r0) / s%grid%dr)

if (any(ij < lbound(s%grid%pt)) .or. any(ij > ubound(s%grid%pt))) then
  call out_io (s_error$, r_name, &
              'Photon position on crystal surface outside of grid bounds for element: ' // ele%name)
  return
endif

! Make small angle approximation

pt => s%grid%pt(ij(1), ij(2))
if (pt%x_pitch == 0 .and. pt%y_pitch == 0 .and. pt%x_pitch_rms == 0 .and. pt%y_pitch_rms == 0) return

h_vec(1:2) = h_vec(1:2) + [pt%x_pitch, pt%y_pitch]
if (pt%x_pitch_rms /= 0 .or. pt%y_pitch /= 0) then
  call ran_gauss (r)
  h_vec(1:2) = h_vec(1:2) + [pt%x_pitch_rms, pt%y_pitch_rms] * r
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
!   w_to_target(3,3) -- real(rp): Rotation matrix from lab to target coords.
!   w_to_ele(3,3)    -- real(rp): Rotation matrix from target to ele coords.
!-

subroutine target_rot_mats (r_center, w_to_target, w_to_ele)

implicit none

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

w_to_target(1,:) = [ cos_theta,           0,       -sin_theta]
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

implicit none

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

end module
