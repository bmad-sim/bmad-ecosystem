module photon_utils_mod

use pointer_to_ele_mod
use pointer_to_branch_mod

implicit none

! This is for passing info into the crystal_diffraction_field_calc routine

type :: crystal_param_struct
  real(rp) cap_gamma, dtheta_sin_2theta, b_eff, wavelength
  real(rp) old_vvec(3), new_vvec(3)
end type

contains

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function photon_type (ele) result (e_type)
!
! Routine to return the type of photon to be tracked: coherent$ or incoherent$.
!
! Input:
!   ele -- ele_struct: Element being tracked through.
!
! Output:
!   e_type -- integer: coherent$ or incoherent$
!-

function photon_type (ele) result (e_type)

type (ele_struct) ele
type (branch_struct), pointer :: branch
integer e_type

! Use

e_type = incoherent$   ! Default

if (associated(ele%branch)) then
  branch => pointer_to_branch(ele)
  e_type = branch%lat%photon_type
endif

end function photon_type

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function z_at_surface (ele, x, y, status) result (z)
!
! Routine return the height (z) of the surface for a particular (x,y) position. 
! Remember: +z points into the element.
!
! Input:
!   ele         -- ele_struct: Element
!   x, y        -- real(rp): coordinates on surface.
!
! Output:
!   z           -- real(rp): z coordinate.
!   status      -- integer: 0 -> Everythin OK.
!                           1 -> Cannot compute z due to point being outside of ellipseoid bounds.
!-

function z_at_surface (ele, x, y, status) result (z)

type (ele_struct), target :: ele
type (photon_surface_struct), pointer :: surf

real(rp) x, y, z, g(3), gs, f
integer status, ix, iy

!

surf => ele%photon%surface
status = 0

if (surf%grid%type == segmented$) then
  call init_surface_segment (x, y, ele)

  z = surf%segment%z0 - (x - surf%segment%x0) * surf%segment%slope_x - &
                        (y - surf%segment%y0) * surf%segment%slope_y

else
  z = 0
  do ix = 0, ubound(surf%curvature_xy, 1)
  do iy = 0, ubound(surf%curvature_xy, 2) - ix
    if (ele%photon%surface%curvature_xy(ix, iy) == 0) cycle
    z = z - surf%curvature_xy(ix, iy) * x**ix * y**iy
  enddo
  enddo


  g = surf%elliptical_curvature
  if (g(3) /= 0) then
    f = -((g(1) * x)**2 + (g(2) * y)**2)
    if (f < -1) then
      status = 1
      return
    endif
    z = z + sqrt_one(f) / g(3)
  endif

  gs = surf%spherical_curvature
  if (gs /= 0) then
    f = -((gs * x)**2 + (gs * y)**2)
    if (f < -1) then
      status = 1
      return
    endif
    z = z + sqrt_one(f) / gs
  endif

endif

end function z_at_surface

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

real(rp) x, y, zt, x0, y0, dx, dy, coef_xx, coef_xy, coef_yy, coef_diag, g(3), gs
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

g = s%elliptical_curvature
if (g(3) /= 0) then
  zt = sqrt(1 - (x0 * g(1))**2 - (y0 * g(2))**2)
  seg%z0 = seg%z0 + sqrt_one(-(g(1) * x)**2 - (g(2) * y)**2) / g(3)
  seg%slope_x = seg%slope_x - x0 * g(1)**2 / (g(3) * zt)
  seg%slope_y = seg%slope_y - y0 * g(2)**2 / (g(3) * zt)
  coef_xx = coef_xx - (g(1)**2 / zt - (x0 * g(1)**2)**2 / zt**3) / (2 * g(3))
  coef_yy = coef_yy - (g(2)**2 / zt - (y0 * g(2)**2)**2 / zt**3) / (2 * g(3))
  coef_xy = coef_xy - (x0 * y0 * (g(1) * g(2))**2 / zt**3) / (g(3))
endif

gs = s%spherical_curvature
if (gs /= 0) then
  zt = sqrt(1 - (x0 * gs)**2 - (y0 * gs)**2)
  seg%z0 = seg%z0 + sqrt_one(-(gs * x)**2 - (gs * y)**2) / gs
  seg%slope_x = seg%slope_x - x0 * gs**2 / (gs * zt)
  seg%slope_y = seg%slope_y - y0 * gs**2 / (gs * zt)
  coef_xx = coef_xx - (gs**2 / zt - (x0 * gs**2)**2 / zt**3) / (2 * gs)
  coef_yy = coef_yy - (gs**2 / zt - (y0 * gs**2)**2 / zt**3) / (2 * gs)
  coef_xy = coef_xy - (x0 * y0 * (gs * gs)**2 / zt**3) / (gs)
endif

! Correct for fact that segment is supported at the corners of the segment and the segment is flat.
! This correction only affects z0 and not the slopes

dx = s%grid%dr(1) / 2
dy = s%grid%dr(2) / 2
coef_xx = coef_xx * dx**2
coef_xy = coef_xy * dx * dy
coef_yy = coef_yy * dy**2
coef_diag = coef_xx + coef_yy - abs(coef_xy)

if (abs(coef_diag) > abs(coef_xx) .and. abs(coef_diag) > abs(coef_yy)) then
  seg%z0 = seg%z0 + coef_diag
else if (abs(coef_xx) > abs(coef_yy)) then
  seg%z0 = seg%z0 + coef_xx
else
  seg%z0 = seg%z0 + coef_yy
endif

end subroutine init_surface_segment 

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine crystal_diffraction_field_calc (cp, ele, thickness, param, p_factor, do_branch_calc, e_field, e_phase, orbit_state, dr, follow_undiffracted)
!
! Routine to compute the photon field after reflection.
!
! Input:
!   cp                  -- crystal_param_struct: Crystal parameters.
!   ele                 -- ele_struct: Crystal element.
!   thickness           -- real(rp): Crystal thickness
!   param               -- lat_param_struct: 
!   p_factor            -- real(rp): Polarization factor.
!   do_branch_calc      -- logical: Calculate probability of branching to alpha or beta branches?
!   follow_undiffracted -- logical, optional: Used with mosaic crystals to calcuate undefracted channel. Default is False.
!
! Output:
!   e_field        -- real(rp): Output field amplitude assuming initial field is 1.
!   e_phase        -- real(rp): Field phase advance.
!   orbit_state    -- integer: Set to lost$ if crystal is to thick to transmit a photon.
!   dr(3)          -- real(rp): (x,y,z) orbit change.
!-

subroutine crystal_diffraction_field_calc (cp, ele, thickness, param, p_factor, do_branch_calc, e_field, e_phase, orbit_state, dr, follow_undiffracted)

implicit none

type (crystal_param_struct) cp
type (ele_struct), target :: ele
type (lat_param_struct) param
type (photon_material_struct), pointer :: pms

real(rp) p_factor, e_field, e_phase, sqrt_b, delta1_0_a, delta1_H_a, delta1_0_b, delta1_H_b
real(rp) s_alpha(3), s_beta(3), dr_alpha(3), dr_beta(3), k_0_a(3), k_h_a(3), k_0_b(3), k_h_b(3), dr(3)
real(rp) kr, k0_im, denom, thickness
real(rp), save :: r_ran

complex(rp) e_rel, e_rel_a, e_rel_b, eta, eta1, f_cmp, xi_0k_a, xi_hk_a, xi_0k_b, xi_hk_b
complex(rp) E_hat_alpha, E_hat_beta, E_hat, exp_factor_a, exp_factor_b
complex(rp), save :: E_hat_alpha_saved, E_hat_beta_saved

integer orbit_state
logical, optional :: follow_undiffracted
logical to_alpha_branch, do_branch_calc, follow_diffracted

! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

e_phase = 0
follow_diffracted = (nint(ele%value(ref_orbit_follows$)) == bragg_diffracted$ .and. .not. logic_option(.false., follow_undiffracted))

pms => ele%photon%material
sqrt_b = sqrt(abs(cp%b_eff))

eta = (cp%b_eff * cp%dtheta_sin_2theta + pms%f_0 * cp%cap_gamma * (1.0_rp - cp%b_eff)/2) / &
                                              (cp%cap_gamma * abs(p_factor) * sqrt_b * pms%f_hkl) 
eta1 = sqrt(eta**2 + sign_of(cp%b_eff))
f_cmp = abs(p_factor) * sqrt_b * cp%cap_gamma * pms%f_hkl / 2

! It can happen that eta - eta1 or eta + eta1 is near zero.

if (abs(eta + eta1) > abs(eta - eta1)) then
  xi_0k_b = -f_cmp * sign_of(cp%b_eff) / (eta + eta1)     ! = f_cmp * (eta - eta1)  beta branch xi
  xi_hk_b = -f_cmp * (eta + eta1) / cp%b_eff              ! = f_cmp / (abs(cp%b_eff) * (eta - eta1))

  xi_0k_a = f_cmp * (eta + eta1)                          ! alpha branch xi
  xi_hk_a = f_cmp / (abs(cp%b_eff) * (eta + eta1))
else
  xi_0k_b = f_cmp * (eta - eta1)
  xi_hk_b = f_cmp / (abs(cp%b_eff) * (eta - eta1))

  xi_0k_a = -f_cmp * sign_of(cp%b_eff) / (eta - eta1)                          ! alpha branch xi
  xi_hk_a = -f_cmp * (eta - eta1) / cp%b_eff
endif

!---------------
! Bragg calc

if (ele%value(b_param$) < 0) then 
  if (abs(eta+eta1) > abs(eta-eta1)) then  ! beta branch excited
    e_rel = -2.0_rp * xi_0k_b / (p_factor * cp%cap_gamma * pms%f_hbar)
  else                                     ! alpha branch excited
    e_rel = -2.0_rp * xi_0k_a / (p_factor * cp%cap_gamma * pms%f_hbar)
  endif

  ! Factor of sqrt_b comes from geometrical change in the transverse width of the photon beam

  if (photon_type(ele) == coherent$) then
    e_field = abs(e_rel) / abs(cp%b_eff)
    e_phase = atan2(aimag(e_rel), real(e_rel))
  else
    e_field = abs(e_rel) / sqrt_b
    e_phase = atan2(aimag(e_rel), real(e_rel))
  endif

!---------------
! Laue calc

else
  e_rel_a = -2.0_rp * xi_0k_a / (p_factor * cp%cap_gamma * pms%f_hbar)
  e_rel_b = -2.0_rp * xi_0k_b / (p_factor * cp%cap_gamma * pms%f_hbar)

  ! Alpha branch

  delta1_0_a = real(xi_0k_a) + (1 - cp%cap_gamma * real(pms%f_0) / 2)
  k_0_a = [cp%old_vvec(1), cp%old_vvec(2), sqrt(delta1_0_a**2 - cp%old_vvec(1)**2 - cp%old_vvec(2)**2)] / cp%wavelength

  delta1_H_a = real(xi_hk_a) + (1 - cp%cap_gamma * real(pms%f_0) / 2)
  k_H_a = [cp%new_vvec(1), cp%new_vvec(2), sqrt(delta1_H_a**2 - cp%new_vvec(1)**2 - cp%new_vvec(2)**2)] / cp%wavelength

  if (follow_diffracted) then
    k0_im = (delta1_H_a / cp%wavelength)**2 * (cp%cap_gamma * imag(pms%f_0) / 2 - aimag(xi_0k_a)) / k_0_a(3)
    exp_factor_a = exp(-twopi * k0_im * thickness) * e_rel_a * e_rel_b / (e_rel_b - e_rel_a)
    s_alpha = k_0_a + abs(e_rel_a)**2 * k_H_a
    dr_alpha = thickness * s_alpha / s_alpha(3)

  else
    k0_im = (delta1_0_a / cp%wavelength)**2 * (cp%cap_gamma * imag(pms%f_0) / 2 - imag(xi_0k_a)) / k_0_a(3)
    exp_factor_a = exp(-twopi * k0_im * thickness) * e_rel_b / (e_rel_b - e_rel_a)
    dr_alpha = thickness * k_0_a / k_0_a(3)
  endif

  ! Beta branch

  delta1_0_b = real(xi_0k_b) + (1 - cp%cap_gamma * real(pms%f_0) / 2)
  k_0_b = [cp%old_vvec(1), cp%old_vvec(2), sqrt(delta1_0_b**2 - cp%old_vvec(1)**2 - cp%old_vvec(2)**2)] / cp%wavelength

  delta1_H_b = real(xi_hk_b) + (1 - cp%cap_gamma * real(pms%f_0) / 2)
  k_H_b = [cp%new_vvec(1), cp%new_vvec(2), sqrt(delta1_H_b**2 - cp%new_vvec(1)**2 - cp%new_vvec(2)**2)] / cp%wavelength

  if (follow_diffracted) then
    k0_im = (delta1_H_b / cp%wavelength)**2 * (cp%cap_gamma * imag(pms%f_0) / 2 - imag(xi_0k_b)) / k_0_b(3)
    exp_factor_b = exp(-twopi * k0_im * thickness) * e_rel_a * e_rel_b / (e_rel_b - e_rel_a)
    s_beta = k_0_b + abs(e_rel_b)**2 * k_H_b
    dr_beta = thickness * s_beta / s_beta(3)
  else
    k0_im = (delta1_H_b / cp%wavelength)**2 * (cp%cap_gamma * imag(pms%f_0) / 2 - imag(xi_0k_b)) / k_0_b(3)
    exp_factor_b = exp(-twopi * k0_im * thickness) * e_rel_a / (e_rel_b - e_rel_a)
    dr_beta = thickness * k_0_b / k_0_b(3)
  endif

  ! If crystal is too thick then mark particle as lost

  if (exp_factor_a == 0 .and. exp_factor_b == 0) then
    orbit_state = lost$
    e_field = 0
    return
  endif

  !--------------------------
  ! Calculate phase and field

  if (photon_type(ele) == coherent$) then
    if (follow_diffracted) then
      kr = -twopi * dot_product(k_h_a, dr_alpha)
      E_hat_alpha = cmplx(cos(kr), sin(kr), rp) * exp_factor_a 

      kr = -twopi * dot_product(k_h_b, dr_beta)
      E_hat_beta  = -cmplx(cos(kr), sin(kr), rp) * exp_factor_b

    else
      kr = -twopi * dot_product(k_0_a, dr_alpha)
      E_hat_alpha = cmplx(cos(kr), sin(kr), rp) * exp_factor_a

      kr = -twopi * dot_product(k_0_b, dr_beta)
      E_hat_beta  = cmplx(cos(kr), sin(kr), rp) * exp_factor_b
    endif

    ! Calculate branching numbers (first pass only)

    if (do_branch_calc) then
      call ran_uniform(r_ran)
      E_hat_alpha_saved = E_hat_alpha
      E_hat_beta_saved = E_hat_beta
    endif

    ! And branch

    denom = abs(E_hat_beta_saved) + abs(E_hat_alpha_saved)
    if (r_ran < abs(E_hat_alpha_saved) / denom) then ! alpha branch
      e_field = denom * abs(E_hat_alpha) / abs(E_hat_alpha_saved)
      e_phase = atan2(aimag(E_hat_alpha), real(E_hat_alpha))
      dr = dr_alpha
    else
      e_field = denom * abs(E_hat_beta) / abs(E_hat_beta_saved)
      e_phase = atan2(aimag(E_hat_beta), real(E_hat_beta))
      dr = dr_beta
    endif

  ! Incoherent

  else

    ! Take dr as average

    dr = (dr_alpha * abs(exp_factor_a) + dr_beta * abs(exp_factor_b)) / (abs(exp_factor_a) + abs(exp_factor_b))

    ! Alpha branch 

    if (follow_diffracted) then
      kr = -twopi * dot_product(k_h_a, dr)
      E_hat_alpha = cmplx(cos(kr), sin(kr), rp) * exp_factor_a 
    else
      kr = -twopi * dot_product(k_0_a, dr)
      E_hat_alpha = cmplx(cos(kr), sin(kr), rp) * exp_factor_a 
    endif

    ! Beta branch

    if (follow_diffracted) then
      kr = -twopi * dot_product(k_h_B, dr)
      E_hat_beta  = -cmplx(cos(kr), sin(kr), rp) * exp_factor_b 
    else
      kr = -twopi * dot_product(k_0_b, dr)
      E_hat_beta  = -cmplx(cos(kr), sin(kr), rp) * exp_factor_b 
    endif

    E_hat = E_hat_alpha + E_hat_beta
    e_field = abs(E_hat)
    e_phase = atan2(aimag(E_hat), real(E_hat))
    dr = (dr_alpha * abs(E_hat_alpha) + dr_beta * abs(E_hat_beta)) / (abs(E_hat_alpha) + abs(E_hat_beta))
  endif
endif

end subroutine crystal_diffraction_field_calc



end module
