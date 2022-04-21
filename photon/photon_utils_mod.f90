module photon_utils_mod

use pointer_to_ele_mod

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
! Function has_curvature (phot_ele) result (curved)
!
! Routine to determine if a surface is potentially curved or is flat.
!
! Input:
!   phot_ele    -- photon_element_struct: From ele%photon
!
! Output:
!   curved      -- logical: Set True if phot_eleace is curved.
!-

function has_curvature (phot_ele) result (curved)

type (photon_element_struct) phot_ele
logical curved

!

curved = (phot_ele%curvature%has_curvature .or. (phot_ele%grid%type == displacement$ .and. phot_ele%grid%active))

end function has_curvature

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
! Function z_at_surface (ele, x, y, err_flag, extend_grid, dz_dxy) result (z)
!
! Routine return the height (z) of the surface for a particular (x,y) position. 
! Remember: +z points into the element.
!
! Input:
!   ele         -- ele_struct: Element
!   x, y        -- real(rp): Photon coordinates on surface.
!   extend_grid -- logical, optional: If a grid is involved and (x, y) is outside of the grid, and 
!                   extend_grid = True: Pretend (x, y) is at edge. Default is False.
!
! Output:
!   z           -- real(rp): z coordinate.
!   err_flag    -- logical: Set True if cannot compute z due to, say, point 
!                             being outside of ellipseoid or grid bounds.
!   dz_dxy(2)   -- real(rp), optional: Surface slope at (x, y).
!-

function z_at_surface (ele, x, y, err_flag, extend_grid, dz_dxy) result (z)

type (ele_struct), target :: ele
type (photon_element_struct), pointer :: ph
type (surface_grid_pt_struct), pointer :: pt

real(rp) x, y, z, g(3), gs, f, dz_dx, dz_dy, xx, yy
real(rp), optional :: dz_dxy(2)

integer ix, iy
logical err_flag
logical, optional :: extend_grid

!

ph => ele%photon
err_flag = .true.
z = 0

if (ph%grid%type == segmented$  .and. ph%grid%active) then
  pt => pointer_to_surface_grid_pt(ele, .true., x, y, ix, iy, extend_grid, xx, yy)
  if (.not. associated(pt)) return
  z = pt%z0 - (xx - pt%x0) * pt%dz_dx - (yy - pt%y0) * pt%dz_dy
  if (present(dz_dxy)) dz_dxy = [pt%dz_dx, pt%dz_dy]

else
  if (ph%grid%type == displacement$) then
    call surface_grid_displacement (ele, x, y, err_flag, z, dz_dxy, extend_grid); if (err_flag) return
  else
    if (present(dz_dxy)) dz_dxy = 0
  endif

  do ix = 0, ubound(ph%curvature%xy, 1)
    do iy = 0, ubound(ph%curvature%xy, 2) - ix
      if (ph%curvature%xy(ix, iy) == 0) cycle
      z = z - ph%curvature%xy(ix, iy) * x**ix * y**iy
      if (present(dz_dxy)) then
        if (ix /= 0) dz_dxy(1) = dz_dxy(1) + ph%curvature%xy(ix, iy) * ix * x**(ix-1) * y**iy
        if (iy /= 0) dz_dxy(2) = dz_dxy(2) + ph%curvature%xy(ix, iy) * iy * x**ix * y**(iy-1)
      endif
    enddo
  enddo

  g = ph%curvature%elliptical
  if (g(3) /= 0) then
    f = -((g(1) * x)**2 + (g(2) * y)**2)
    if (f < -1) return
    z = z + sqrt_one(f) / g(3)
    if (present(dz_dxy)) dz_dxy = dz_dxy - (2 * sqrt_one(f,1) / g(3)) * g(1:2)**2 * [x, y]
  endif

  gs = ph%curvature%spherical
  if (gs /= 0) then
    f = -((gs * x)**2 + (gs * y)**2)
    if (f < -1) return
    z = z + sqrt_one(f) / gs
    if (present(dz_dxy)) dz_dxy = dz_dxy - (2 * sqrt_one(f,1) / gs) * gs**2 * [x, y]
  endif
endif

err_flag = .false.

end function z_at_surface

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Function pointer_to_surface_grid_pt (ele, nearest, x, y, ix, iy, extend_grid, xx, yy) result (pt)
!
! Routine to point to the grid point struct associated with point (x,y).
!
! Note: If nearest = True, the grid boundary is a length dr/2 from the boundary grid points.
!
! Input:
!   ele         -- ele_struct: Element containing the grid
!   nearest     -- logical: If True, return pointer to nearest grid point. 
!                           If False, return pointer to the grid point lower and left of (x,y).
!   x, y        -- real(rp): Photon position.
!   extend_grid -- logical, optional: If (x,y) past grid pretend (x,y) is at grid boundary.
!                   Default is False.
!
! Output:
!   ix, iy      -- integer, optional: Grid point index.
!   pt          -- grid_point_struct: Pointer to grid point. 
!                   Will not be associated if (x,y) outside the grid.
!   xx, yy      -- real(rp), optional: Set equal to (x, y) except if (x,y) is outside of the grid.
!                   In this case, (xx, yy) will be set to be on the nearest grid boundary point.
!-

function pointer_to_surface_grid_pt (ele, nearest, x, y, ix, iy, extend_grid, xx, yy) result (pt)

type (ele_struct), target :: ele
type (surface_grid_struct), pointer :: grid
type (surface_grid_pt_struct), pointer :: pt

real(rp) x, y, xh, yh, ff
real(rp), optional :: xx, yy

integer, optional :: ix, iy
integer kx, ky, nx0, ny0, nx1, ny1

logical, optional :: extend_grid
logical nearest, outside

character(*), parameter :: r_name = 'pointer_to_surface_grid_pt'

!

grid => ele%photon%grid

if (nearest) then
  ff = 0.5_rp
else
  ff = 0
endif

nx0 = lbound(grid%pt, 1);  nx1 = ubound(grid%pt, 1)
ny0 = lbound(grid%pt, 2);  ny1 = ubound(grid%pt, 2)

if (x < grid%pt(nx0,ny0)%x0 - ff * grid%dr(1)) then
  xh = grid%pt(nx0,ny0)%x0 - ff * grid%dr(1)
elseif (x > grid%pt(nx1,ny1)%x0 + ff * grid%dr(1)) then
  xh = grid%pt(nx1,ny1)%x0 + ff * grid%dr(1)
else
  xh = x
endif

if (y < grid%pt(nx0,ny0)%y0 - ff * grid%dr(2)) then
  yh = grid%pt(nx0,ny0)%y0 - ff * grid%dr(2)
elseif (y > grid%pt(nx1,ny1)%y0 + ff * grid%dr(2)) then
  yh = grid%pt(nx1,ny1)%y0 + ff * grid%dr(2)
else
  yh = y
endif

if (present(xx)) xx = xh
if (present(yy)) yy = yh

outside = (x /= xh .or. y /= yh)

if (.not. logic_option(.false., extend_grid) .and. outside) then
  call out_io (s_info$, r_name, 'Photon position: (\2f12.8\) is outside of grid for: ' // ele%name, r_array = [x, y])
  pt => null()
  return
endif

if (nearest) then
  kx = nint((xh - grid%r0(1)) / grid%dr(1))
  ky = nint((yh - grid%r0(2)) / grid%dr(2))
  kx = min(max(kx, nx0), nx1)   ! Can happen due to roundoff
  ky = min(max(ky, ny0), ny1)   ! Can happen due to roundoff
else
  kx = floor((xh - grid%r0(1)) / grid%dr(1))
  ky = floor((yh - grid%r0(2)) / grid%dr(2))
  kx = min(max(kx, nx0), nx1-1)   ! Can happen due to roundoff
  ky = min(max(ky, ny0), ny1-1)   ! Can happen due to roundoff
endif

pt => grid%pt(kx, ky)

if (present(ix)) ix = kx
if (present(iy)) iy = ky

end function pointer_to_surface_grid_pt

!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------
!+
! Subroutine surface_grid_displacement (ele, x, y, err_flag, z, dz_dxy, extend_grid)
!
! Routine to add in the z displacement defined by the grid 
!
! Input:
!   ele           -- ele_struct: Element containing the grid
!   x, y          -- real(rp): Photon coords at surface.
!   extend_grid   -- logical, optional: If (x,y) past grid pretend (x,y) is at grid boundary.
!                     Default is False.
!
! Output
!   err_flag      -- logical: Set True if there is a problem.
!   z             -- real(rp): surface height at (x, y).
!   dz_dxy(2)     -- real(rp), optional: Surface slope at (x, y).
!-

subroutine surface_grid_displacement (ele, x, y, err_flag, z, dz_dxy, extend_grid)

use super_recipes_mod

type (ele_struct), target :: ele
type (photon_element_struct), pointer :: ph
type (surface_grid_pt_struct), pointer :: pt00, pt01, pt10, pt11

real(rp) x, y
real(rp), optional :: dz_dxy(2)
real(rp) z, xx, yy, dz_dx, dz_dy
integer ix, iy
logical err_flag
logical, optional :: extend_grid

!

ph => ele%photon

if (.not. ph%grid%active) then
  z = 0
  if (present(dz_dxy)) dz_dxy = 0
  err_flag = .false.
  return
endif

!

err_flag = .true.
pt00 => pointer_to_surface_grid_pt(ele, .false., x, y, ix, iy, extend_grid, xx, yy)
if (.not. associated(pt00)) return

pt01 => ph%grid%pt(ix,iy+1)
pt10 => ph%grid%pt(ix+1,iy)
pt11 => ph%grid%pt(ix+1,iy+1)

call super_bicubic_interpolation([pt00%z0, pt10%z0, pt11%z0, pt01%z0], [pt00%dz_dx, pt10%dz_dx, pt11%dz_dx, pt01%dz_dx], &
        [pt00%dz_dy, pt10%dz_dy, pt11%dz_dy, pt01%dz_dy], [pt00%d2z_dxdy, pt10%d2z_dxdy, pt11%d2z_dxdy, pt01%d2z_dxdy], &
        pt00%x0, pt11%x0, pt00%y0, pt11%y0, xx, yy, z, dz_dx, dz_dy)

if (present(dz_dxy)) dz_dxy = [dz_dx, dz_dy]

err_flag = .false.

end subroutine surface_grid_displacement

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
