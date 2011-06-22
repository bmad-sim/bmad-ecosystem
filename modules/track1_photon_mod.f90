module track1_photon_mod

use bmad_struct
use bmad_interface

! This is for passing info to the photon_hit_crystal routine used by zbrent and reflection_point

type (ele_struct), private, pointer, save :: ele_com
type (coord_struct), private, pointer, save :: start_com
real (rp), private, pointer, save :: m_in_com(:,:)
real (rp), private, pointer, save :: k_in_norm_com(:)

private photon_hit_crystal

contains

!---------------------------------------------------------------------------
!+
! Subroutine track_a_photon (ele,param,end)
!
! Routine to track reflection from a crystal.
!
! Input:
!   ele    -- ele_struct: Element tracking through.
!   param  -- lat_param_struct: lattice parameters.
!   end    -- Coord_struct: phase-space coords to be transformed
!
! Output:
!   end    -- Coord_struct: final phase-space coords
!-

subroutine track_a_photon (ele, param, end)

implicit none

type (ele_struct), target:: ele
type (coord_struct), target:: end
type (lat_param_struct) :: param

real(rp), target::m_in(3,3)
real(rp), target::k_in_norm(3)
real(rp) f, p_factor, sin_alpha, cos_alpha, sin_psi, cos_psi, wavelength
real(rp) cos_g, sin_g, cos_tc, sin_tc
real(rp) h_norm(3), k_out_norm(3), e_tot, pc
real(rp) cap_gamma, gamma_0, gamma_h, b_err, dtheta_sin_2theta, b_eff

real(rp) m_out(3,3), y_out(3), x_out(3), k_out(3)
real(rp) test, temp_vec(3), direction(3),z_len

real(rp) refpoint(3),slope,c2,c3,c4,curveangle,cos_c,sin_c,curverot(3,3)

complex(rp) f0, fh, f0_g, eta, eta1, f_cmp, xi_0k, xi_hk, e_rel, e_rel2

wavelength = ele%value(ref_wavelength$) / (1 + end%vec(6))

! (px, py, sqrt(1-px^2+py^2)) coords are with respect to laboratory reference trajectory.
! Convert this vector to k_in_norm which are coords with respect to crystal surface.
! k_in_norm is normalized to 1.

sin_g = sin(ele%value(graze_angle_in$))
cos_g = cos(ele%value(graze_angle_in$))
f = sqrt (1 - end%vec(2)**2 - end%vec(4)**2)

k_in_norm(1) =  cos_g * end%vec(2) + f * sin_g
k_in_norm(2) = end%vec(4)
k_in_norm(3) = -sin_g * end%vec(2) + f * cos_g

! Construct m_in and m_out matrices of basis vectors
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


! To account for curvature, compute the reflection point

ele_com => ele
start_com => end
m_in_com => m_in
k_in_norm_com => k_in_norm
call reflection_point(ele,end,refpoint,z_len)

! Compute the slope of the crystal at that point

c2 = ele%value(c2_curve_tot$)
c3 = ele%value(c3_curve_tot$)
c4 = ele%value(c4_curve_tot$)
slope = -2.0_rp*c2*refpoint(3) - 3.0_rp*c3*refpoint(3)**2 - 4.0_rp*c4*refpoint(3)**4
curveangle = atan(slope)
cos_c=cos(curveangle)
sin_c=sin(curveangle)

! Form a rotation matrix

curverot = reshape([cos_c, 0.0_rp, -sin_c, 0.0_rp, 1.0_rp, 0.0_rp, sin_c, 0.0_rp, cos_c], [3,3])

! Construct h_norm = H vector * wavelength, including a rotation due to curvature
sin_alpha = sin(ele%value(alpha_angle$))
cos_alpha = cos(ele%value(alpha_angle$))
sin_psi = sin(ele%value(psi_angle$))
cos_psi = cos(ele%value(psi_angle$))
h_norm = [-cos_alpha, sin_alpha * sin_psi, sin_alpha * cos_psi] * wavelength / ele%value(d_spacing$)
h_norm = matmul(curverot,h_norm)
  
! k_out_norm is the outgoing wavevector outside the crystal
k_out_norm = k_in_norm + h_norm
k_out_norm = matmul(transpose(curverot),k_out_norm)
k_out_norm(1) = - sqrt( 1 - k_out_norm(2)**2 - k_out_norm(3)**2)
k_out_norm = matmul(curverot,k_out_norm)

!======= (x,px,y,py,x,pz) - Phase Space Calculations
!Translate to outgoing basis
direction = matmul(transpose(m_out), k_out_norm)
end%vec(2) = direction(1)
end%vec(4) = direction(2)

!Compute position in phase space, backpropagating the ray
temp_vec = matmul(transpose(m_out),refpoint)
temp_vec = temp_vec - direction *temp_vec(3)

end%vec(1) = temp_vec(1)
end%vec(3) = temp_vec(2)
! %vec(5) doesn't include phase change due to wave nature of radiation
end%vec(5) = temp_vec(3)+z_len


!======== Calculate phase and intensity
cap_gamma = r_e * wavelength**2 / (pi * ele%value(v_unitcell$)) 
gamma_0 = k_in_norm(1)
gamma_h = k_out_norm(1)

b_eff = gamma_0 / gamma_h
dtheta_sin_2theta = dot_product(h_norm + 2.0_rp * k_in_norm, h_norm) *0.5_rp
f0 = cmplx(ele%value(f0_re$), ele%value(f0_im$)) 
fh = cmplx(ele%value(fh_re$), ele%value(fh_im$))
f0_g = cap_gamma * f0 *0.5_rp

! For the x direction
! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

p_factor = cos(2.0_rp*ele%value(graze_angle_in$))
eta = (-b_eff * dtheta_sin_2theta + f0_g * (1.0_rp - b_eff)) / &
          (cap_gamma * abs(p_factor) * sqrt(abs(b_eff)) * fh) 
eta1 = sqrt(eta**2 + sign(1.0_rp, b_eff))
f_cmp = abs(p_factor) * sqrt(abs(b_eff)) * cap_gamma * fh * 0.5_rp
if (abs(eta+eta1) > abs(eta-eta1)) then
  xi_0k = f_cmp * (eta - eta1)
  xi_hk = f_cmp / (abs(b_eff) * (eta - eta1))
else        
  xi_0k = f_cmp * (eta + eta1)
  xi_hk = f_cmp / (abs(b_eff) * (eta + eta1))
endif

! relative electric field, or reflectivity calculated in 2 equivalent ways

e_rel = -2.0_rp * xi_0k / (p_factor * cap_gamma * fh)
! e_rel2 = sqrt(xi_0k/xi_hk) ! assert = e_rel

end%e_field_x = end%e_field_x * abs(e_rel)
end%phase_x = atan2(aimag(e_rel),real(e_rel))+end%phase_x

! For the y direction
! Construct xi_0k = xi_0 / k and xi_hk = xi_h / k

p_factor = 1.0_rp
eta = (-b_eff * dtheta_sin_2theta + f0_g * (1.0_rp - b_eff)) / &
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
e_rel = -2.0_rp * xi_0k / (p_factor * cap_gamma * fh)
! e_rel2 = sqrt(xi_0k/xi_hk)! assert = e_rel

end%e_field_y = end%e_field_y * abs(e_rel)
end%phase_y = atan2(aimag(e_rel),real(e_rel))+end%phase_y

end subroutine track_a_photon

!---------------------------------------------------------------------------
!+
! Subroutine reflection_point (ele, end,refpoint, z_len)
!
! Routine to compute position where crystal reflects in crystal coordinates.
!
! Input:
!   
!
! Output:
!   refpoint(3) -- Real(rp): point of reflection in crystal coordinates
!   z_len       -- Real(rp): distance for photon to propagate to refpoint
!-

subroutine reflection_point (ele, end, refpoint, z_len)

use nr, only: zbrent
implicit none
type (ele_struct), target:: ele
type (coord_struct), target:: end
real (rp),target ::z_len
real (rp), target :: m_in(3,3)
real (rp), target :: refpoint(3)
real (rp) :: vec(3),cos_g, sin_g,f
real (rp) :: x1,x2,fa,fb

!Assume flat crystal, compute z required to hit the intersection
!Choose a Bracket of 1m around this point.

vec = matmul(m_in_com, [start_com%vec(1),0.0_rp,0.0_rp])
x2 = vec(1) / k_in_norm_com(1)

x1 = x2 - 0.5_rp
x2 = x2 + 0.5_rp

z_len = zbrent (photon_hit_crystal, x1, x2, 1d-10)

! Compute the intersection point
refpoint = z_len * k_in_norm_com + matmul( m_in_com , [start_com%vec(1), start_com%vec(3), 0.0_rp] )

end subroutine reflection_point


!---------------------------------------------------------------------------
!+
! Function photon_hit_crystal (z_len) result (delta_h)
! 
! Routine to be used as an argument in zbrent. Propagates
! photon forward by a distance z_len. Returns delta_h = x-x0
! where x0 is the height of the crystal surface. 
! Returns a positive number if photon is inside crystal, 
! negative if outside
!
! Input:
!   z_len   -- Real(rp): Place to position the photon.
!
! Output:
!   delta_h -- Real(rp): height of photon below surface in crystal coordinates
!-

function photon_hit_crystal (z_len) result (delta_h)
real(rp), intent(in) :: z_len
real(rp) :: delta_h
real(rp) :: c2, c3, c4
real(rp) :: point(3), h, r

c2 = ele_com%value(c2_curve_tot$)
c3 = ele_com%value(c3_curve_tot$)
c4 = ele_com%value(c4_curve_tot$)

point = k_in_norm_com * z_len + matmul( m_in_com, [start_com%vec(1), start_com%vec(3), 0.0_rp] )
h = point(1)
r = point(3)
delta_h = h - c2*r*r - c3*r*r*r - c4*r*r*r*r

end function photon_hit_crystal

end module
