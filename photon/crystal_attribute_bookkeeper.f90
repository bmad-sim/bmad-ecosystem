!+
! Subroutine crystal_attribute_bookkeeper (ele)
! 
! Routine to orient the crystal element.
!
! Input:
!   ele   -- Ele_struct: Crystal element.
!     %value(bragg_angle$)     -- Bragg angle.
!     %value(e_tot$)           -- Photon reference energy.
!     %value(alpha_angle$)     -- H-vector orientation
!     %value(d_spacing$)       -- Crystal plane spacing.
!     ... etc ...
!
! Output
!   ele       -- Ele_struct: Crystal element.
!     %value(bragg_angle_in$)
!     %value(bragg_angle_out$)
!     %value(tilt_corr$)
!     ... etc.
!-

subroutine crystal_attribute_bookkeeper (ele)

use photon_utils_mod, dummy => crystal_attribute_bookkeeper

implicit none

type (ele_struct), target :: ele
type (photon_material_struct), pointer :: pms

real(rp) lambda, gamma, delta1, lambda_in, d, alpha, psi, theta0
real(rp) bragg_angle_in, ang_tot, k0_x_norm, k0_y_norm, k0_z_norm
real(rp) kh_x_norm, kh_y_norm, kh_z_norm, h_x_norm, h_y_norm, h_z_norm
real(rp) ent_kh_x_norm, ent_kh_y_norm, ent_kh_z_norm, b_param, beta
real(rp) cos_graze_in, sin_graze_in, s_vec(3), k0_norm(3), h_vec(3), dtheta_sin_2theta
real(rp) p_factor, f

complex(rp) eta, eta1, f_cmp, xi_0k, xi_hk

character(*), parameter :: r_name = 'crystal_attribute_bookkeeper'

! If the photon energy or the bragg angle has not been set then cannot do the calc yet.

if (ele%value(e_tot$) == 0) return
if (ele%value(bragg_angle$) == 0) return

pms => ele%photon%material
b_param = ele%value(b_param$)

if (b_param == 0) then
  call out_io (s_error$, r_name, 'B_PARAM NOT SET FOR CRYSTAL ELEMENT: ' // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

lambda = ele%value(ref_wavelength$)
gamma = lambda**2 * r_e / (pi * ele%value(v_unitcell$))
delta1 = 1 - gamma * real(pms%f_0) / 2
lambda_in = lambda / delta1
d = ele%value(d_spacing$)

alpha = ele%value(alpha_angle$)
psi   = ele%value(psi_angle$)

h_x_norm = -sin(alpha) * cos(psi)
h_y_norm =  sin(alpha) * sin(psi) 
h_z_norm = -cos(alpha) 

pms%h_norm = [h_x_norm, h_y_norm, h_z_norm]

! Compute bragg_angle_in

beta = lambda_in / (2 * d)

if (b_param < 0) then ! Bragg
  theta0 = asin((-beta * h_z_norm - h_x_norm * sqrt(h_x_norm**2 + h_z_norm**2 - beta**2)) / (h_x_norm**2 + h_z_norm**2))
  bragg_angle_in = acos(delta1 * cos(theta0))
  k0_x_norm = -cos(bragg_angle_in)
  k0_z_norm =  sin(bragg_angle_in)

else                              ! Laue
  theta0 = asin((-beta * h_x_norm + h_z_norm * sqrt(h_x_norm**2 + h_z_norm**2 - beta**2)) / (h_x_norm**2 + h_z_norm**2))
  bragg_angle_in = asin(delta1 * sin(theta0))
  k0_x_norm =  sin(bragg_angle_in)
  k0_z_norm =  cos(bragg_angle_in)
endif

ele%value(bragg_angle_in$) = bragg_angle_in

! kh_norm is the normalized k_H vector.
! Obtained from: k_H = H + K_0 + q*n_surface

kh_x_norm = lambda * h_x_norm/d + k0_x_norm
kh_y_norm = lambda * h_y_norm/d
f = 1 - kh_x_norm**2 - kh_y_norm**2
if (f < 0) then   ! Can happen due to roundoff
  kh_z_norm = 0
else
  kh_z_norm = sqrt(f)
endif
if (ele%value(b_param$) < 0) kh_z_norm = -kh_z_norm   ! Bragg

! ent_kh_norm is kh_norm in entrance coordinates 

ent_kh_x_norm = kh_z_norm * k0_x_norm - kh_x_norm * k0_z_norm
ent_kh_y_norm = kh_y_norm
ent_kh_z_norm = kh_z_norm * k0_z_norm + kh_x_norm * k0_x_norm

! There is no tilt correction if we are following the non-diffracted beam.

if (nint(ele%value(ref_orbit_follows$)) == bragg_diffracted$) then
  ele%value(tilt_corr$) = atan2(ent_kh_y_norm, ent_kh_x_norm)
else
  ele%value(tilt_corr$) = 0
endif

ang_tot = atan2(sqrt(ent_kh_x_norm**2 + ent_kh_y_norm**2), ent_kh_z_norm)
ele%value(bragg_angle_out$) = ang_tot - bragg_angle_in

! Displacement due to finite crystal thickness for Laue diffraction.

if (ele%value(b_param$) > 0) then
  k0_norm = [sin(bragg_angle_in), 0.0_rp, cos(bragg_angle_in)]
  h_vec = lambda * pms%h_norm / d
  dtheta_sin_2theta = -dot_product(h_vec + 2 * k0_norm, h_vec) / 2

  p_factor = 1
  eta = (b_param * dtheta_sin_2theta + pms%f_0 * gamma * (1.0_rp - b_param)/2) / (gamma * abs(p_factor) * sqrt(b_param) * pms%f_hkl) 
  eta1 = sqrt(eta**2 + 1)
  f_cmp = abs(p_factor) * sqrt(b_param) * gamma * pms%f_hkl / 2
  ele%value(pendellosung_period_sigma$) = cos(bragg_angle_in) * lambda_in / (2 * real(f_cmp * eta1)) 

  p_factor = cos(ang_tot)
  eta = (b_param * dtheta_sin_2theta + pms%f_0 * gamma * (1.0_rp - b_param)/2) / (gamma * abs(p_factor) * sqrt(b_param) * pms%f_hkl) 
  eta1 = sqrt(eta**2 + 1)
  f_cmp = abs(p_factor) * sqrt(b_param) * gamma * pms%f_hkl / 2
  ele%value(pendellosung_period_pi$) = cos(bragg_angle_in) * lambda_in / (2 * real(f_cmp * eta1)) 

  select case (nint(ele%value(ref_orbit_follows$)))
  case (undiffracted$)
    ! reference orbit direction is same as k_0 (outside)
    pms%l_ref = [tan(bragg_angle_in), 0.0_rp, 1.0_rp] * ele%value(thickness$) 
  case (forward_diffracted$, bragg_diffracted$)
    p_factor = (1 + cos(ang_tot)) / 2   ! Average of sigma and pi polarization factors

    eta = (b_param * dtheta_sin_2theta + pms%f_0 * gamma * (1.0_rp - b_param)/2) / &
                                                  (gamma * abs(p_factor) * sqrt(b_param) * pms%f_hkl) 
    eta1 = sqrt(eta**2 + 1)
    f_cmp = abs(p_factor) * sqrt(b_param) * gamma * pms%f_hkl / 2

    xi_0k = f_cmp * (eta + eta1)                         ! alpha branch xi
    xi_hk = f_cmp / (b_param * (eta + eta1))

    s_vec = real(xi_hk) * k0_norm + real(xi_0k) * (k0_norm + h_vec)
    pms%l_ref = s_vec * ele%value(thickness$) / s_vec(3) 
  end select

else              ! Bragg
  pms%l_ref = 0
  ele%value(pendellosung_period_sigma$) = 0
  ele%value(pendellosung_period_pi$) = 0
endif

ele%value(l$) = norm2(pms%l_ref)

!

ele%value(ref_cap_gamma$) = gamma

if (ang_tot == 0) then
  ele%value(darwin_width_sigma$) = 0
  ele%value(darwin_width_pi$) = 0
  ele%value(dbragg_angle_de$) = 0
else
  ele%value(darwin_width_sigma$) = 2 * gamma * real(pms%f_hkl) / (abs(sin(ang_tot)) * sqrt(abs(ele%value(b_param$))))
  ele%value(darwin_width_pi$) = ele%value(darwin_width_sigma$) * abs(cos(ang_tot))
  ele%value(dbragg_angle_de$) = -lambda / (2 * d * ele%value(e_tot$) * cos(abs(ang_tot)/2))
endif

! Reflectivity table

call setup_reflect_table(ele%photon%reflectivity_table_sigma)
call setup_reflect_table(ele%photon%reflectivity_table_pi)

!------------------------------------------------------------------------
contains

subroutine setup_reflect_table(rt)

type (photon_reflect_table_struct) :: rt
real(rp) lambda, beta
integer i, ne

!

if (allocated(rt%p_reflect)) then
  ne = size(rt%energy)
  call re_allocate (rt%bragg_angle, ne)
  do i = 1, ne
    lambda = c_light * h_planck / rt%energy(i)
    beta = lambda / (2 * ele%value(d_spacing$))
    if (b_param < 0) then ! Bragg
      rt%bragg_angle(i) = asin((-beta * h_z_norm - h_x_norm * sqrt(h_x_norm**2 + h_z_norm**2 - beta**2)) / (h_x_norm**2 + h_z_norm**2))
    else
      rt%bragg_angle(i) = asin((-beta * h_x_norm + h_z_norm * sqrt(h_x_norm**2 + h_z_norm**2 - beta**2)) / (h_x_norm**2 + h_z_norm**2))
    endif
  enddo
endif

end subroutine setup_reflect_table

end subroutine
