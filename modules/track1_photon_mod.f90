module track1_photon_mod

use bmad_struct
use bmad_interface

! This is for passing info to the photon_hit_func routine used by zbrent

type (ele_struct), private, pointer, save :: ele_com
type (coord_struct), private, pointer, save :: start_com

private photon_hit_crystal

contains

!---------------------------------------------------------------------------
!+
! Subroutine track_a_photon (ele,param,end)
!
! Routine to track reflection from a crystal.
!
! Input:
!   
!
! Output:
!   
!-

subroutine track_a_photon (ele, param, end)

implicit none

type (ele_struct) ele
type (coord_struct) end
type (lat_param_struct) :: param

real(rp) f, p_factor, sin_alpha, cos_alpha, sin_psi, cos_psi, wavelength
real(rp) cos_g, sin_g, cos_tc, sin_tc
real(rp) k_in_norm(3), h_norm(3), k_out_norm(3), e_tot, pc
real(rp) cap_gamma, gamma_0, gamma_h, b_err, dtheta_sin_2theta, b_eff

real(rp) m_in(3,3) , m_out(3,3), y_out(3), x_out(3), k_out(3)
real(rp) test, nn, mm, temp_vec(3)

real(rp) refpoint(3),slope,c2,c3,c4,curveangle,cos_c,sin_c,curverot(3,3)

complex(rp) f0, fh, f0_g, eta, eta1, f_cmp, xi_0k, xi_hk, e_rel, e_rel2

wavelength = ele%value(ref_wavelength$) / (1 + end%vec(6))

call reflection_point(ele,end,refpoint)
c2 = ele%value(c2_curve$)
c3 = ele%value(c3_curve$)
c4 = ele%value(c4_curve$)
slope = -2*c2*refpoint(3) - 3*c3*refpoint(3)**2 - 4*c4*refpoint(3)**4
curveangle = atan(slope)
cos_c=cos(curveangle)
sin_c=sin(curveangle)
curverot = reshape([cos_c, 0.0_rp, -sin_c, 0.0_rp, 1.0_rp, 0.0_rp, sin_c, 0.0_rp, cos_c], [3,3])

! (px, py, sart(1-px^2+py^2)) coords are with respect to the incoming reference trajectory.
! Convert this vector to k_in_norm which are coords with respect to crystal surface.
! k_in_norm is incoming wavevector * wavelength so has unit length.

sin_g = sin(ele%value(graze_angle_in$))
cos_g = cos(ele%value(graze_angle_in$))
f = sqrt (1 - end%vec(2)**2 - end%vec(4)**2)

k_in_norm(1) =  cos_g * end%vec(2) + f * sin_g
k_in_norm(2) = end%vec(4)
k_in_norm(3) = -sin_g * end%vec(2) + f * cos_g

! Construct h_norm = H vector * wavelength

sin_alpha = sin(ele%value(alpha_angle$))
cos_alpha = cos(ele%value(alpha_angle$))
sin_psi = sin(ele%value(psi_angle$))
cos_psi = cos(ele%value(psi_angle$))

h_norm = [-cos_alpha, sin_alpha * sin_psi, sin_alpha * cos_psi] * wavelength / ele%value(d_spacing$)
h_norm = matmul(curverot,h_norm)

! Construct m_in and m_out matrices
! m_in = [x_in y_in k_in]
! m_out = [x_out y_out k_out]

m_in = reshape([cos_g, 0.0_rp, -sin_g, 0.0_rp, 1.0_rp, 0.0_rp, sin_g, 0.0_rp, cos_g], [3,3])

sin_tc = sin(ele%value(tilt_corr$))
cos_tc = cos(ele%value(tilt_corr$))

! y_out = inverse(m_in) . m_tiltcorr . m_in . (0,1,0)

y_out = matmul(m_in, [0.0_rp, 1.0_rp, 0.0_rp])
y_out = matmul(reshape([cos_tc, sin_tc, 0.0_rp, -sin_tc, cos_tc, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp], [3, 3]),  y_out)
y_out = matmul(transpose(m_in), y_out)

! x_out = vector orthogonal to y and z

x_out(1) = y_out(2)*ele%value(nz_out$)-y_out(3)*ele%value(ny_out$)
x_out(2) = -y_out(1)*ele%value(nz_out$)+y_out(3)*ele%value(nx_out$)
x_out(3) = y_out(1)*ele%value(ny_out$)-y_out(2)*ele%value(nx_out$)
  
k_out(1) = ele%value(nx_out$)
k_out(2) = ele%value(ny_out$)
k_out(3) = ele%value(nz_out$)
test = dot_product( k_out, y_out ) ! assert 0
m_out = reshape( [x_out, y_out, k_out], [3, 3])

  
! k_out_norm is the outgoing wavevector outside the crystal

k_out_norm = k_in_norm + h_norm
k_out_norm = matmul(transpose(curverot),k_out_norm)
k_out_norm(1) = - sqrt( 1 - k_out_norm(2)**2 - k_out_norm(3)**2)
k_out_norm = matmul(curverot,k_out_norm)

temp_vec = matmul(transpose(m_out), k_out_norm)

end%vec(2) = temp_vec(1)
end%vec(4) = temp_vec(2)

!======= Position in Phase Space
  
temp_vec = matmul(m_in, [end%vec(1), end%vec(3), 0.0_rp])
nn = -dot_product( [1.0_rp, 0.0_rp, 0.0_rp], temp_vec )
nn = nn/dot_product( [1.0_rp, 0.0_rp, 0.0_rp], k_in_norm )

temp_vec = temp_vec + nn * k_in_norm
  
mm = -dot_product(k_out, temp_vec)
mm = nn/dot_product(k_out, k_out_norm)
  
temp_vec = temp_vec + mm * k_out_norm

temp_vec = matmul(transpose(m_out), temp_vec)

end%vec(1) = temp_vec(1)
end%vec(3) = temp_vec(2)
! %vec(5) doesn't include phase change due to wave nature of radiation
end%vec(5) = nn + mm 


!======== Calculate phase and intensity
cap_gamma = r_e * wavelength**2 / (pi * ele%value(v_unitcell$)) 
gamma_0 = k_in_norm(1)
gamma_h = k_out_norm(1)

b_eff = gamma_0 / gamma_h
dtheta_sin_2theta = dot_product(h_norm + 2 * k_in_norm, h_norm) / 2
f0 = cmplx(ele%value(f0_re$), ele%value(f0_im$)) 
fh = cmplx(ele%value(fh_re$), ele%value(fh_im$))
f0_g = cap_gamma * f0 / 2

! For the x direction
! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

p_factor = cos(2*ele%value(graze_angle_in$))
eta = (-b_eff * dtheta_sin_2theta + f0_g * (1 - b_eff)) / &
          (cap_gamma * abs(p_factor) * sqrt(abs(b_eff)) * fh) 
eta1 = sqrt(eta**2 + sign(1.0_rp, b_eff))
f_cmp = abs(p_factor) * sqrt(abs(b_eff)) * cap_gamma * fh / 2
if (abs(eta+eta1) > abs(eta-eta1)) then
  xi_0k = f_cmp * (eta - eta1)
  xi_hk = f_cmp / (abs(b_eff) * (eta - eta1))
else        
  xi_0k = f_cmp * (eta + eta1)
  xi_hk = f_cmp / (abs(b_eff) * (eta + eta1))
endif

! relative electric field, or reflectivity calculated in 2 equivalent ways

e_rel = -2 * xi_0k / (p_factor * cap_gamma * fh)
e_rel2 = sqrt(xi_0k/xi_hk) ! assert = e_rel

end%e_field_x = end%e_field_x * abs(e_rel)
end%phase_x = atan2(aimag(e_rel),real(e_rel))+end%phase_x

! For the y direction
! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

p_factor = 1
eta = (-b_eff * dtheta_sin_2theta + f0_g * (1 - b_eff)) / &
          (cap_gamma * abs(p_factor) * sqrt(abs(b_eff)) * fh) 
eta1 = sqrt(eta**2 + sign(1.0_rp, b_eff))

f_cmp = abs(p_factor) * sqrt(abs(b_eff)) * cap_gamma * fh / 2
if (abs(eta+eta1) > abs(eta-eta1)) then
  xi_0k = f_cmp * (eta - eta1)
  xi_hk = f_cmp / (abs(b_eff) * (eta - eta1))
else        
  xi_0k = f_cmp * (eta + eta1)
  xi_hk = f_cmp / (abs(b_eff) * (eta + eta1))
endif

! relative electric field, or reflectivity calculated in 2 equivalent ways
e_rel = -2 * xi_0k / (p_factor * cap_gamma * fh)
e_rel2 = sqrt(xi_0k/xi_hk)! assert = e_rel

end%e_field_y = end%e_field_y * abs(e_rel)
end%phase_y = atan2(aimag(e_rel),real(e_rel))+end%phase_y

end subroutine track_a_photon

!---------------------------------------------------------------------------
!+
! Subroutine reflection_point (ele,end,point)
!
! Routine to compute position where crystal reflects in crystal coordinates.
!
! Input:
!   
!
! Output:
!   
!-

subroutine reflection_point (ele, end, point)

use nr, only: zbrent
implicit none
type (ele_struct), target:: ele
type (coord_struct), target:: end

real (rp),target :: point(3)
real(rp) :: m_in(3,3), k_in_norm(3),cos_g, sin_g,f
real (rp) ::x1,x2,z_len

ele_com=>ele
start_com=>end

x1=-100.0
x2=100.0
z_len = zbrent (photon_hit_crystal, x1, x2, 1d-10)

! Construct m_in matrices
! m_in = [x_in y_in k_in]
sin_g = sin(ele_com%value(graze_angle_in$))
cos_g = cos(ele_com%value(graze_angle_in$))
m_in = reshape([cos_g, 0.0_rp, -sin_g, 0.0_rp, 1.0_rp, 0.0_rp, sin_g, 0.0_rp, cos_g], [3,3])
! Compute the incoming vector
f = sqrt(1 - start_com%vec(2)**2 - start_com%vec(4)**2)
k_in_norm(1) =  cos_g * start_com%vec(2) + f * sin_g
k_in_norm(2) = start_com%vec(4)
k_in_norm(3) = -sin_g * start_com%vec(2) + f * cos_g
! Compute the intersection point
point(1) = k_in_norm(1)*z_len + m_in(1,1)*start_com%vec(1)+m_in(1,3)*start_com%vec(3)
point(2) = k_in_norm(2)*z_len + m_in(2,1)*start_com%vec(1)+m_in(2,3)*start_com%vec(3)
point(3) = k_in_norm(3)*z_len + m_in(3,1)*start_com%vec(1)+m_in(3,3)*start_com%vec(3)


end subroutine reflection_point


!---------------------------------------------------------------------------
!+
! Function photon_hit_crystal (z_len) result (d_radius)
! 
! Routine to be used as an argument in zbrent
!
! Input:
!   z_len -- Real(rp): Place to position the photon.
!
! Output:
!-

function photon_hit_crystal (z_len) result (delta_h)
real(rp), intent(in) :: z_len
real(rp) :: delta_h
real(rp) :: c2, c3, c4, cos_g, sin_g
real(rp) :: m_in(3,3), k_in_norm(3), h, r, f

c2 = ele_com%value(c2_curve$)
c3 = ele_com%value(c3_curve$)
c4 = ele_com%value(c4_curve$)

! Construct m_in matrices
! m_in = [x_in y_in k_in]
sin_g = sin(ele_com%value(graze_angle_in$))
cos_g = cos(ele_com%value(graze_angle_in$))
m_in = reshape([cos_g, 0.0_rp, -sin_g, 0.0_rp, 1.0_rp, 0.0_rp, sin_g, 0.0_rp, cos_g], [3,3])


f = sqrt(1 - start_com%vec(2)**2 - start_com%vec(4)**2)
k_in_norm(1) =  cos_g * start_com%vec(2) + f * sin_g
k_in_norm(2) = start_com%vec(4)
k_in_norm(3) = -sin_g * start_com%vec(2) + f * cos_g

h = k_in_norm(1)*z_len + m_in(1,1)*start_com%vec(1)+m_in(1,3)*start_com%vec(3)
r = k_in_norm(3)*z_len + m_in(3,1)*start_com%vec(1)+m_in(3,3)*start_com%vec(3)

delta_h = h-c2*r*r-c3*r*r*r-c4*r*r*r*r

end function photon_hit_crystal

end module
