module make_mat6_mod

use bmad_utils_mod

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_add_pitch (x_pitch_tot, y_pitch_tot, orientation, mat6)
!
! Subroutine to modify a first order transfer matrix to include the affect
! of an element pitch. Note that this routine does not correct the 0th order
! part of the map. It is assumed that on input the transfer map
! does not include the affect of any pitches.
!
! Modules needed:
!   use bmad
!
! Input:
!   x_pitch_tot -- Real(rp): Horizontal pitch
!   y_pitch_tot -- Real(rp): Vertical pitch
!   orientation -- integer: Element longitudinal orientation. +1 or -1.
!   mat6(6,6)   -- Real(rp): 1st order part of the transfer map (Jacobian).
!
! Output:
!   mat6(6,6) -- Real(rp): 1st order xfer map with pitches.
!-

subroutine mat6_add_pitch (x_pitch_tot, y_pitch_tot, orientation, mat6)

implicit none

real(rp) mat6(6,6), x_pitch_tot, y_pitch_tot
integer orientation

!

if (x_pitch_tot == 0 .and. y_pitch_tot == 0) return

! The equations below are performing matrix multiplication. The original matrix
! is being multiplied from left and right by matrices that correspond to the pitches. 
! The pitch matrices are obtained by differentiating the corresponding equations in   
! the offset_particle subroutine. The (i,j) numbers mentioned as comments refer to  
! the non-zero elements present in the pitch matrices. 

mat6(:,6) = mat6(:,6) - mat6(:,2) * orientation * x_pitch_tot ! (2,6)
mat6(:,1) = mat6(:,1) + mat6(:,5) * orientation * x_pitch_tot ! (5,1)

mat6(:,6) = mat6(:,6) - mat6(:,4) * orientation * y_pitch_tot ! (4,6)
mat6(:,3) = mat6(:,3) + mat6(:,5) * orientation * y_pitch_tot ! (5,3)

mat6(2,:) = mat6(2,:) + orientation * x_pitch_tot * mat6(6,:) ! (2,6)
mat6(5,:) = mat6(5,:) - orientation * x_pitch_tot * mat6(1,:) ! (5,1)

mat6(4,:) = mat6(4,:) + orientation * y_pitch_tot * mat6(6,:) ! (4,6)
mat6(5,:) = mat6(5,:) - orientation * y_pitch_tot * mat6(3,:) ! (5,3)

end subroutine mat6_add_pitch

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine quad_mat2_calc (k1, length, rel_p, mat2, z_coef, dz_dpz_coef)
!
! Subroutine to calculate the 2x2 transfer matrix for a quad for one plane. 
!
! Modules needed:
!   use bmad
!
! Input:
!   k1       -- Real(rp): Quad strength: k1 > 0 ==> defocus.
!   length   -- Real(rp): Quad length
!   rel_p    -- Real(rp), Relative momentum P/P0.      
!
! Output:
!   mat2(2,2)      -- Real(rp): Transfer matrix.
!   z_coef(3)      -- Real(rp), optional: Coefficients for calculating the
!                       the change in z position:
!                          z = Integral [-(px/(1+pz))^2/2 ds]
!                            = c(1) * x_0^2 + c(2) * x_0 * px_0 + c(3) * px_0^2 
!   dz_dpz_coef(3) -- Real(rp), optional: Coefficients for calculating the
!                       the mat6(5,6) Jacobian matrix element:
!                         dz_dpz = c(1) * x_0^2 + c(2) * x_0 * px_0 + c(3) * px_0^2 
!-

subroutine quad_mat2_calc (k1, length, rel_p, mat2, z_coef, dz_dpz_coef)

implicit none

real(rp) length, mat2(:,:), cx, sx
real(rp) k1, sqrt_k, sk_l, k_l2, zc(3), dsx, dcx, rel_p
real(rp), optional :: z_coef(3), dz_dpz_coef(3)

!

sqrt_k = sqrt(abs(k1))
sk_l = sqrt_k * length

if (abs(sk_l) < 1d-10) then
  k_l2 = k1 * length**2
  cx = 1 + k_l2 / 2
  sx = (1 + k_l2 / 6) * length
elseif (k1 < 0) then       ! focus
  cx = cos(sk_l)
  sx = sin(sk_l) / sqrt_k
else                       ! defocus
  cx = cosh(sk_l)
  sx = sinh(sk_l) / sqrt_k
endif

mat2(1,1) = cx
mat2(1,2) = sx / rel_p
mat2(2,1) = k1 * sx * rel_p
mat2(2,2) = cx

!

if (present(z_coef) .or. present(dz_dpz_coef)) then
  zc(1) = k1 * (-cx * sx + length) / 4
  zc(2) = -k1 * sx**2 / (2 * rel_p)
  zc(3) = -(cx * sx + length) / (4 * rel_p**2)
  if (present(z_coef)) z_coef = zc
endif

! dz_dpz_coef

if (present(dz_dpz_coef)) then

  if (abs(sk_l) < 1d-10) then
    dcx = -k_l2 / (2 * rel_p)
    dsx = -k_l2 * length / (6 * rel_p)
  else
    dcx = -k1 * sx * length / (2 * rel_p)
    dsx = (sx - length * cx) / (2 * rel_p)
  endif

  dz_dpz_coef(1) =   -zc(1)/rel_p - k1 * (cx * dsx + dcx * sx) / 4
  dz_dpz_coef(2) = -2*zc(2)/rel_p - k1 * sx * dsx/rel_p
  dz_dpz_coef(3) = -2*zc(3)/rel_p - (cx * dsx + dcx * sx) / (4 * rel_p**2)

endif

end subroutine quad_mat2_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine sol_quad_mat6_calc (ks, k1, length, orb, mat6, dz_coef)
!
! Subroutine to calculate the transfer matrix for a combination 
! solenoid/quadrupole element (without a tilt).
!
! Modules Needed:
!   use bmad
!
! Input:
!   ks      -- Real(rp): Solenoid strength.
!   k1      -- Real(rp): Quadrupole strength.
!   length  -- Real(rp): Sol_quad length.
!   orb(6)  -- Real(rp): Orbit at beginning of the sol_quad.
!
! Output:
!   mat6(6,6)    -- Real(rp): Transfer matrix across the sol_quad.
!   dz_coef(4,4) -- Real(rp), optional: coefs for dz calc
!                     dz = sum_ij dz_coef(i,j) * orb(i) * orb(j)
!-

subroutine sol_quad_mat6_calc (ks_in, k1_in, s_len, orb, m, dz_coef)

implicit none

real(rp) ks, k1, s_len
real(rp) m(6,6)
real(rp) orb(6)
real(rp), optional :: dz_coef(4,4)

real(rp) ks2, s, c, snh, csh, rel_p, ks_in, k1_in
real(rp) darg1, alpha, alpha2, beta, beta2, f, q, r, a, b
real(rp) df, dalpha2, dalpha, dbeta2, dbeta, darg
real(rp) dC, dCsh, dS, dSnh, dq, dr, da, db, k1_2
real(rp) ks3, fp, fm, dfm, dfp, df_f, ug, orbit(6)
real(rp) s1, s2, snh1, snh2, dsnh1, dsnh2, ds1, ds2
real(rp) coef1, coef2, dcoef1, dcoef2, ks4
real(rp) t4(4,4), ts(4,4), m0(6,6), xp_start, xp_end, yp_start, yp_end
real(rp) dt4(4,4), dts(4,4), dz_co(4,4), d_dz_co(4,4), dm(4,4)
real(rp) r_orb(4), d_orb(4), tsd(4,4), dtsd(4,4)
real(rp) d2_f, dug, d2_ug, d2_C, d2_fp, d2_fm, d2_Csh, d2_q, d2_s1, d2_s2 
real(rp) d2_r, d2_snh1, d2_snh2, d2_a, d2_b, d2_coef1, d2_coef2
real(rp) d2_alpha, d2_beta, d2_arg, d2_arg1, d2_S, d2_Snh, d2_f_f, rk

! Calculation is done in (x, x', y, y') coordinates and then converted
! to (x, p_x, y, p_y) coordinates.

rel_p = 1 + orb(6)

orbit = orb
orbit(2) = orb(2) / rel_p
orbit(4) = orb(4) / rel_p

k1 = k1_in / rel_p
ks = ks_in / rel_p

k1_2 = k1**2
ks2 = ks**2
ks3 = ks2 * ks 
ks4 = ks2*ks2
f = sqrt(ks4 + 4*k1_2)
ug = 1 / (4*f)
alpha2 = (f + ks2) / 2; alpha = sqrt(alpha2)

if (abs(k1) < 1d-2 * f) then
  rk = (k1 / ks2)**2
  beta2 = ks2 * (rk - rk**2 + 2 * rk**3) 
else
  beta2  = (f - ks2) / 2
endif
beta  = sqrt(beta2)

S = sin(alpha*s_len)                              
C = cos(alpha*s_len)
Snh = sinh(beta*s_len)
Csh = cosh(beta*s_len)
q = 2 * beta2  + 2*k1
r = 2 * alpha2 - 2*k1
a = 2 * alpha2 + 2*k1
b = 2 * beta2  - 2*k1 
fp = f + 2*k1
fm = f - 2*k1

S1 = S * alpha
S2 = S / alpha

Snh1 = Snh * beta

if (abs(beta) < 1d-10) then
  Snh2 = s_len
else
  Snh2 = Snh / beta
endif

coef1 = ks2*r + 4*k1*a
coef2 = ks2*q + 4*k1*b

! m0 is the transfer matrix in (x, x', y, y') space.

call mat_make_unit(m0)
             
m0(1,1) = 2*ug * (fp*C + fm*Csh)
m0(1,2) = (2*ug/k1) * (q*S1 - r*Snh1)
m0(1,3) = (ks*ug/k1) * (-b*S1 + a*Snh1)
m0(1,4) = 4*ug*ks * (-C + Csh)

m0(2,1) = -(ug/2) * (coef1*S2 + coef2*Snh2)
m0(2,2) = m0(1,1)             
m0(2,3) = ug*ks3 * (C - Csh)
m0(2,4) = ug*ks * (a*S2 + b*Snh2)

m0(3,1) = -m0(2,4)
m0(3,2) = -m0(1,4)
m0(3,3) = 2*ug * (fm*C + fp*Csh)  
m0(3,4) = 2*ug * (r*S2 + q*Snh2)

m0(4,1) = -m0(2,3)     
m0(4,2) = -m0(1,3)
m0(4,3) = (ug/(2*k1)) * (-coef2*S1 + coef1*Snh1)
m0(4,4) = m0(3,3)

!

df      = -2 * (ks4 + 2*k1_2) / f
dalpha2 = df/2 - ks2
dalpha  = (df/2 - ks2)/(2*alpha)
if (k1_2 < 1d-5 * ks4) then
  dbeta   = abs(k1**3/(ks3*ks2)) * (-1 + 3.5 * k1_2/ks4)
else
  dbeta   = (ks2 + df/2)/(2*beta)
endif
dbeta2  = 2 * beta * dbeta
darg    = s_len*dalpha
darg1   = s_len*dbeta         
dC      = -darg*S
dCsh    = darg1*Snh
dS      = darg*C
dSnh    = darg1*Csh
dq      = -2*k1 + 2*dbeta2
dr      =  2*k1 + 2*dalpha2
da      = -2*k1 + 2*dalpha2
db      =  2*k1 + 2*dbeta2
dfp = df - 2*k1
dfm = df + 2*k1
df_f =  -df/f

dS1 = dS * alpha + S * dalpha
dS2 = dS / alpha - S * dalpha / alpha2

dSnh1 = dSnh * beta + Snh * dbeta

if (k1_2 < 1d-5 * ks4) then
  dSnh2 = k1**4 * s_len**3 * (-1/3.0d0 + (40 - ks2 * s_len**2) * k1_2 / (30 * ks4)) / ks3**2
else
  dSnh2 = dSnh / beta - Snh * dbeta / beta2
endif

dcoef1 = -2*ks2*r + ks2*dr - 4*k1*a + 4*k1*da
dcoef2 = -2*ks2*q + ks2*dq - 4*k1*b + 4*k1*db                     

! t4(i,j) is dm(i,j)/dE at pz = 0 and is used to calculate the m0(x,6) terms.

t4(1,1) = m0(1,1)*df_f + 2*ug*(fp*dC + C*dfp + fm*dCsh + Csh*dfm)
t4(1,2) = m0(1,2)*df_f + (2*ug/k1) * (dq*S1 + q*dS1 - dr*Snh1 - r*dSnh1)
t4(1,3) = m0(1,3)*df_f + (ks*ug/k1)*(-db*S1 - b*dS1 + da*Snh1 + a*dSnh1)
t4(1,4) = m0(1,4)*(df_f - 2) + 4*ks*ug*(-dC + dCsh) 

t4(2,1) = m0(2,1)*(df_f + 1) - &
            (ug/2)*(dcoef1*S2 + coef1*dS2 + dcoef2*Snh2 + coef2*dSnh2)
t4(2,2) = t4(1,1)
t4(2,3) = m0(2,3)*(df_f - 2) + ks3*ug*(dC - dCsh) 
t4(2,4) = m0(2,4)*(df_f - 1) + ug*ks*(da*S2 + a*dS2 + db*Snh2 + b*dSnh2)

t4(3,1) = -t4(2,4)
t4(3,2) = -t4(1,4)
t4(3,3) = m0(3,3)*df_f + 2*ug*(fm*dC + C*dfm + fp*dCsh + Csh*dfp)
t4(3,4) = m0(3,4)*(df_f - 1) + 2*ug*(dr*S2 + r*dS2 + dq*Snh2 + q*dSnh2)

t4(4,1) = -t4(2,3)        
t4(4,2) = -t4(1,3)
t4(4,3) = m0(4,3)*(df_f + 2) + &
             (ug/(2*k1))*(-dcoef2*S1 - coef2*dS1 + dcoef1*Snh1 + coef1*dSnh1)
t4(4,4) = t4(3,3)

m0(1:4,6) = matmul(t4(1:4,1:4), orbit(1:4))

! m(5,6) calc

d2_f  = 2*(ks4+2*k1_2)*df/f**2 + 8*(ks4+k1_2)/f
dug  = -df/(4*f**2)
d2_ug = -d2_f/(4*f**2) + df**2/(2*f**3)
d2_f_f = (df/f)**2 - d2_f / f

d2_alpha  = -dalpha**2/alpha + (d2_f/2+2*ks2)/(2*alpha)
if (k1_2 < 1d-5 * ks4) then
  d2_beta = -abs(-3*k1**3/ks**5+5*k1**3/ks**5)
else
  d2_beta = -dbeta**2/beta + (d2_f/2-2*ks2)/(2*beta)
endif
d2_arg  = s_len*d2_alpha
d2_arg1 = s_len*d2_beta
d2_C    = -d2_arg*S - darg*dS
d2_Csh  = d2_arg1*Snh + darg1*dSnh
d2_S    = d2_arg*C + darg*dC
d2_Snh  = d2_arg1*Csh + darg1*dCsh
d2_q    =  2*k1 - 4*ks2 + d2_f
d2_r    = -2*k1 + 4*ks2 + d2_f
d2_a    =  2*k1 + 4*ks2 + d2_f
d2_b    = -2*k1 - 4*ks2 + d2_f
d2_fp   = d2_f + 2*k1
d2_fm   = d2_f - 2*k1

d2_S1 = S*d2_alpha+alpha*d2_S+2*dS*dalpha
d2_S2 = -dS*dalpha/alpha2+dalpha*S*dalpha2/alpha2**2-S*d2_alpha/alpha2+d2_S/alpha-dalpha*dS/alpha2

d2_Snh1 = Snh*d2_beta+beta*d2_Snh+2*dSnh*dbeta

if (k1_2 < 1d-5 * ks4) then
  d2_Snh2 = 4*k1**4*s_len**3/(3*ks**6)-2*k1**4*s_len**3/ks**6
else
  d2_Snh2 = -dSnh*dbeta/beta2+dbeta*Snh*dbeta2/beta2**2-Snh*d2_beta/beta2+d2_Snh/beta-dbeta*dSnh/beta2
endif

d2_coef1 = 4*ks2*r-4*ks2*dr+ks2*d2_r+4*k1*a-8*k1*da+4*k1*d2_a
d2_coef2 = 4*ks2*q-4*ks2*dq+ks2*d2_q+4*k1*b-8*k1*db+4*k1*d2_b

! First deal with m elements that are not affected by rel_p

dt4(1,1) = t4(1,1)*df_f + m0(1,1)*d2_f_f + 2*dug*(fp*dC + C*dfp + fm*dCsh + Csh*dfm) + &
            2*ug*(d2_fp*C + 2*dfp*dC + fp*d2_C + d2_Csh*fm + 2*dCsh*dfm + Csh*d2_fm)
dt4(1,2) = (t4(1,2) + m0(1,2))*df_f + m0(1,2)*d2_f_f + &
           2*(dug/k1 + ug/k1) * (dq*S1 + q*dS1 - dr*Snh1 - r*dSnh1) + &
           (2*ug/k1) * (d2_q*S1 + 2*dq*dS1 + q*d2_S1 - d2_r*Snh1 - 2*dr*dSnh1 - r*d2_Snh1)
dt4(1,3) = t4(1,3)*df_f + m0(1,3)*d2_f_f + (ks*dug/k1) *(-db*S1 - b*dS1 + da*Snh1 + a*dSnh1) + &
           (ks*ug/k1)*(-d2_b*S1 - 2*db*dS1 - b*d2_S1 + d2_a*Snh1 + 2*da*dSnh1 + a*d2_Snh1)
dt4(1,4) = (t4(1,4) + m0(1,4))*(df_f - 2) + m0(1,4)*d2_f_f + 4*(ks*dug - ks*ug)*(-dC + dCsh) + &
           4*ks*ug*(-d2_C + d2_Csh)

dt4(2,1) = (t4(2,1) - m0(2,1))*(df_f + 1) + m0(2,1)*d2_f_f - &
           (dug/2)*(dcoef1*S2 + coef1*dS2 + dcoef2*Snh2 + coef2*dSnh2) - &
           (ug/2)*(d2_coef1*S2 + 2*dcoef1*dS2 + coef1*d2_S2 + d2_coef2*Snh2 + 2*dcoef2*dSnh2 + coef2*d2_Snh2)
dt4(2,2) = dt4(1,1)
dt4(2,3) = (t4(2,3) - m0(2,3))*(df_f - 2) + m0(2,3)*d2_f_f + (ks3*dug - 3*ks3*ug)*(dC - dCsh) + &
           ks3*ug*(d2_C - d2_Csh)
dt4(2,4) = t4(2,4)*(df_f - 1) + m0(2,4)*d2_f_f + (dug*ks - ug*ks)*(da*S2 + a*dS2 + db*Snh2 + b*dSnh2) + &
           ug*ks*(d2_a*S2 + 2*da*dS2 + a*d2_S2 + d2_b*Snh2 + 2*db*dSnh2 + b*d2_Snh2)

dt4(3,1) = -dt4(2,4)
dt4(3,2) = -dt4(1,4)
dt4(3,3) = 4*(fm*dC+fp*dCSh+C*dfm+Csh*dfp)*dug+2*(C*fm+Csh*fp)*d2_ug + &
           2*ug*(2*dC*dfm+2*dCsh*dfp+fm*d2_C+fp*d2_Csh+C*d2_fm+Csh*d2_fp)
dt4(3,4) = (t4(3,4) + m0(3,4))*(df_f - 1) + m0(3,4)*d2_f_f + 2*dug*(dr*S2 + r*dS2 + dq*Snh2 + q*dSnh2) + &
           2*ug*(d2_r*S2 + 2*dr*dS2 + r*d2_S2 + d2_q*Snh2 + 2*dq*dSnh2 + q*d2_Snh2)

dt4(4,1) = -dt4(2,3)
dt4(4,2) = -dt4(1,3)
dt4(4,3) = (t4(4,3) - m0(4,3))*(df_f + 2) + m0(4,3)*d2_f_f + &
           ((dug + ug)/(2*k1))*(-dcoef2*S1 - coef2*dS1 + dcoef1*Snh1 + coef1*dSnh1) + &
           (ug/(2*k1))*(-d2_coef2*S1 - 2*dcoef2*dS1 - coef2*d2_S1 + d2_coef1*Snh1 + 2*dcoef1*dSnh1 + coef1*d2_Snh1)
dt4(4,4) = dt4(3,3)

! The m(5,6) term is computed 

ts(1:4,1) = -t4(2,1:4)
ts(1:4,2) =  t4(1,1:4)
ts(1:4,3) = -t4(4,1:4)
ts(1:4,4) =  t4(3,1:4)

dts(1:4,1) = -dt4(2,1:4)
dts(1:4,2) =  dt4(1,1:4)
dts(1:4,3) = -dt4(4,1:4)
dts(1:4,4) =  dt4(3,1:4)

tsd = ts
tsd(1:4,2) = ts(1:4,2) / rel_p
tsd(1:4,4) = ts(1:4,4) / rel_p

dtsd(1:4,1) = dts(1:4,1) / rel_p
dtsd(1:4,2) = (dts(1:4,2) - ts(1:4,2)) / rel_p**2
dtsd(1:4,3) = dts(1:4,3) / rel_p
dtsd(1:4,4) = (dts(1:4,4) - ts(1:4,4)) / rel_p**2

r_orb = [orbit(1), orbit(2)*rel_p, orbit(3), orbit(4)*rel_p]
d_orb = [0.0_rp, -orbit(2)/rel_p, 0.0_rp, -orbit(4)/rel_p]

! dz = Sum_ij dz_coef(i,j) * orb(i) * orb(j)

if (present(dz_coef)) then
  dz_coef = matmul (ts, m0(1:4,1:4)) / 2
  dz_coef(:,2) = dz_coef(:,2) / rel_p 
  dz_coef(:,4) = dz_coef(:,4) / rel_p 
  dz_coef(2,:) = dz_coef(2,:) / rel_p 
  dz_coef(4,:) = dz_coef(4,:) / rel_p 
endif

! energy corrections

if (all(orbit == 0) .and. .not. present(dz_coef)) then
  m = m0
  return
endif

m = m0

m(1,2) = m0(1,2) / rel_p
m(1,4) = m0(1,4) / rel_p

m(2,1) = m0(2,1) * rel_p
m(2,3) = m0(2,3) * rel_p

m(3,2) = m0(3,2) / rel_p
m(3,4) = m0(3,4) / rel_p

m(4,1) = m0(4,1) * rel_p
m(4,3) = m0(4,3) * rel_p

m(1,6) = m0(1,6) / rel_p
m(3,6) = m0(3,6) / rel_p

!

dm = t4 / rel_p
dm(2,1) = t4(2,1)
dm(2,3) = t4(2,3)
dm(4,1) = t4(4,1)
dm(4,3) = t4(4,3)
dm(1,2) = t4(1,2) / rel_p**2
dm(1,4) = t4(1,4) / rel_p**2
dm(3,2) = t4(3,2) / rel_p**2
dm(3,4) = t4(3,4) / rel_p**2

m(5,6) = (dot_product(matmul(matmul(d_orb, tsd), m(1:4,1:4)), r_orb) + &
          dot_product(matmul(matmul(orbit(1:4), dtsd), m(1:4,1:4)), r_orb) + &
          dot_product(matmul(matmul(orbit(1:4), tsd), dm(1:4,1:4)), r_orb)) / 2

! The m(5,x) terms follow from the symplectic condition.

m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(2,1) - m(4,6)*m(3,1) + m(3,6)*m(4,1)
m(5,2) = -m(2,6)*m(1,2) + m(1,6)*m(2,2) - m(4,6)*m(3,2) + m(3,6)*m(4,2)
m(5,3) = -m(2,6)*m(1,3) + m(1,6)*m(2,3) - m(4,6)*m(3,3) + m(3,6)*m(4,3)
m(5,4) = -m(2,6)*m(1,4) + m(1,6)*m(2,4) - m(4,6)*m(3,4) + m(3,6)*m(4,4)

end subroutine sol_quad_mat6_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine sbend_body_with_k1_map (ele, param, n_step, start, end, mat6)
!
! Subroutine to calculate for a single step the transfer matrix and/or 
! ending coordinates for a sbend with a finite k1 but without a tilt.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele      -- Ele_struct: Sbend element.
!   param    -- Lat_param_struct: Branch parameters.
!   n_step   -- Integer: Number of steps to divide the bend into.
!               Only one step is taken by this routine.
!   start    -- coord_struct: Orbit at beginning of the bend.
!
! Output:
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix across the sol_quad.
!   end        -- coord_struct, optional: Ending coordinates.
!-

subroutine sbend_body_with_k1_map (ele, param, n_step, start_orb, end_orb, mat6)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) start_orb
type (coord_struct), optional :: end_orb

real(rp) g, g_err, length
real(rp), optional :: mat6(6,6)
real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, s_y, c_y, r(6)
real(rp) z0, z1, z2, z11, z12, z22, z33, z34, z44
real(rp) dom_x, dom_xx, dx_c, dc_x, ds_x, dom_y, dom_yy, dc_y, ds_y, dcs_x, dcs_y
real(rp) g_tot, rel_p, rel_p2, charge_dir
real(rp) rel_pc, px, py, pxy2, pz

integer n_step

! Degenerate case

charge_dir = relative_tracking_charge(start_orb, param) * ele%orientation * start_orb%direction

k_1 = ele%value(k1$) * charge_dir
g = ele%value(g$)
g_tot = (g + ele%value(g_err$)) * charge_dir
g_err = g_tot - g
length = ele%value(l$) / n_step

!

g_tot = g + g_err
rel_p = (1 + start_orb%vec(6))
rel_p2 = rel_p**2


k_x = k_1 + g * g_tot
x_c = (g * rel_p - g_tot) / k_x

om_x = sqrt(abs(k_x) / rel_p)
om_y = sqrt(abs(k_1) / rel_p)

tau_x = -sign (1.0_rp, k_x)
tau_y =  sign (1.0_rp, k_1)

arg = om_x * length
if (arg < 1d-6) then
  s_x = (1 + tau_x * arg**2 / 6) * length
  c_x = 1 + tau_x * arg**2 / 2
  z2 = g * length**2 / (2 * rel_p)
elseif (k_x > 0) then
  s_x = sin(arg) / om_x
  c_x = cos(arg)
  z2 = tau_x * g * (1 - c_x) / (rel_p * om_x**2)
else
  s_x = sinh(arg) / om_x
  c_x = cosh(arg)
  z2 = tau_x * g * (1 - c_x) / (rel_p * om_x**2)
endif

arg = om_y * length
if (arg < 1d-6) then
  s_y = (1 + tau_y * arg**2 / 6) * length
  c_y = 1 + tau_y * arg**2 / 2
elseif (k_1 < 0) then
  s_y = sin(om_y * length) / om_y
  c_y = cos(om_y * length)
else
  s_y = sinh(om_y * length) / om_y
  c_y = cosh(om_y * length)
endif

r = start_orb%vec
r(1) = r(1) - x_c

!

z0  = -g * x_c * Length
z1  = -g * s_x
z11 = tau_x * om_x**2 * (length - c_x * s_x) / 4
z12 = -tau_x * om_x**2 * s_x**2 / (2 * rel_p) 
z22 = -(length + c_x * s_x) / (4 * rel_p2) 
z33 = tau_y * om_y**2 * (length - c_y * s_y) / 4
z34 = -tau_y * om_y**2 * s_y**2 / (2 * rel_p) 
z44 = -(length + c_y * s_y) / (4 * rel_p2)

! Jacobian matrix

if (present(mat6)) then

  dom_x = -om_x / (2 * rel_p)
  dom_xx = -1 / (2 * rel_p)  ! [d(om_x) / d(p_z)] / om_x
  dx_c = g / k_x
  dc_x = tau_x * s_x * om_x * dom_x * length
  ds_x = (c_x * length - s_x) * dom_xx 

  dom_y = -om_y / (2 * rel_p)
  dom_yy = -1 / (2 * rel_p)  ! [d(om_y) / d(p_z)] / om_y
  dc_y = tau_y * s_y * om_y * dom_y * length
  ds_y = (c_y * length - s_y) * dom_yy

  dcs_x = c_x * ds_x + dc_x * s_x
  dcs_y = c_y * ds_y + dc_y * s_y

  mat6 = 0

  mat6(1,1) = c_x
  mat6(1,2) = s_x / rel_p
  mat6(1,6) = dx_c * (1 - c_x) + dc_x * r(1) + &
            ds_x * r(2) / rel_p - s_x * r(2) / rel_p2
  mat6(2,1) = tau_x * om_x**2 * rel_p * s_x
  mat6(2,2) = c_x
  mat6(2,6) = tau_x * r(1) * 2 * om_x * dom_x * rel_p * s_x + &
              tau_x * r(1) * om_x**2 * s_x + &
              tau_x * r(1) * om_x**2 * rel_p * ds_x - &
              tau_x * dx_c * om_x**2 * rel_p * s_x + dc_x * r(2)

  mat6(3,3) = c_y
  mat6(3,4) = s_y / rel_p
  mat6(3,6) = dc_y * r(3) + ds_y * r(4) / rel_p - s_y * r(4) / rel_p2
  mat6(4,3) = tau_y * om_y**2 * rel_p * s_y
  mat6(4,4) = c_y
  mat6(4,6) = tau_y * r(3) * 2 * om_y * dom_y * rel_p * s_y + &
              tau_y * r(3) * om_y**2 * s_y + &
              tau_y * r(3) * om_y**2 * rel_p * ds_y + &
              dc_y * r(4)

  mat6(5,1) = z1 + 2 * z11 * r(1) +     z12 * r(2)  
  mat6(5,2) = z2 +     z12 * r(1) + 2 * z22 * r(2)
  mat6(5,3) =      2 * z33 * r(3) +     z34 * r(4)  
  mat6(5,4) =          z34 * r(3) + 2 * z44 * r(4)
  mat6(5,5) = 1
  mat6(5,6) = -dx_c * (z1 + 2 * z11 * r(1) + z12 * r(2)) - & 
              g * length * dx_c - g * ds_x * r(1) - &           ! dz0 & dz1
              tau_x * g * dc_x * r(2) / (om_x**2 * rel_p) - &   ! dz2
              (z11 / rel_p + tau_x * om_x**2 * dcs_x / 4) * r(1)**2 - &
              (2 * z12 / rel_p + tau_x * om_x**2 * s_x * ds_x / rel_p) * r(1) * r(2) - &
              (2 * z22 / rel_p + dcs_x / (4 * rel_p2)) * r(2)**2 - &
              (z33 / rel_p + tau_y * om_y**2 * dcs_y / 4) * r(3)**2 - &
              (2 * z34 / rel_p + tau_y * om_y**2 * s_y * ds_y / rel_p) * r(3) * r(4) - &
              (2 * z44 / rel_p + dcs_y / (4 * rel_p2)) * r(4)**2

  mat6(6,6) = 1

endif

! Ending coords

if (present(end_orb)) then
  end_orb%vec(1) = c_x * r(1) + s_x * r(2) / rel_p + x_c
  end_orb%vec(2) = tau_x * om_x**2 * rel_p * s_x * r(1) + c_x * r(2)
  end_orb%vec(3) = c_y * r(3) + s_y * r(4) / rel_p
  end_orb%vec(4) = tau_y * om_y**2 * rel_p * s_y * r(3) + c_y * r(4)
  end_orb%vec(5) = r(5) + z0 + z1 * r(1) + z2 * r(2) + &
                  z11 * r(1)**2 + z12 * r(1) * r(2) + z22 * r(2)**2 + &
                  z33 * r(3)**2 + z34 * r(3) * r(4) + z44 * r(4)**2 
endif

end subroutine sbend_body_with_k1_map

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine multipole_kick_mat (knl, tilt, orbit, factor, mat6)
!
! Subroutine to return the multipole kick components needed to
! construct the transfer matrix.
! This routine is not meant for general use.
!
! Input:
!   knl(0:)   -- Real(rp): Strength of multipoles
!   tilt(0:)  -- Real(rp): Tilt of multipoles
!   orbit     -- Coord_struct: coordinates of particle around which the
!                  multipole kick matrix is computed.
!   factor    -- real(rp): Factor to scale knl by.
!
! Output:
!   mat6(6,6) -- Real(rp): matrix with kick values at mat6(2:4:2, 1:3:2).
!                 The rest of the matrix is untouched.
!-

subroutine multipole_kick_mat (knl, tilt, orbit, factor, mat6)

implicit none

type (coord_struct) orbit

real(rp) mat6(6,6), kmat1(4,4), factor
real(rp) knl(0:), tilt(0:)

integer n

!                        

mat6(2:4:2, 1:3:2) = 0
if (orbit%vec(1) == 0 .and. orbit%vec(3) == 0 .and. knl(1) == 0) return

do n = 1, ubound(knl, 1)
  if (knl(n) /= 0) then
    call mat4_multipole (knl(n), tilt(n), n, orbit, kmat1)
    mat6(2:4:2, 1:3:2) = mat6(2:4:2, 1:3:2) + factor * kmat1(2:4:2, 1:3:2)
  endif
enddo

end subroutine multipole_kick_mat

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat4_multipole (knl, tilt, n, orbit, kick_mat)
!
! Subroutine to find the kick matrix (Jacobian) due to a multipole.
! This routine is not meant for general use.
!
! Input:
!   orbit   -- Coord_struct: coordinates of particle
!   knl     -- Real(rp): Strength of multipole
!   tilt    -- Real(rp): Tilt of multipole
!
! Output:
!   kick_mat(4,4) -- Real(rp): Kick matrix (Jacobian) at orbit.
!-


subroutine mat4_multipole (knl, tilt, n, orbit, kick_mat)

implicit none

type (coord_struct) orbit

real(rp) x_pos, y_pos, x, y, knl, tilt
real(rp) sin_ang, cos_ang, mat(2,2), rot(2,2)
real(rp) kick_mat(4,4)

integer m, n

! init

kick_mat = 0
forall (m = 1:4) kick_mat(m,m) = 1

x_pos = orbit%vec(1)
y_pos = orbit%vec(3)
         
! simple case

if (knl == 0 .or. (x_pos == 0 .and. y_pos == 0 .and. n /= 1)) then
  kick_mat(2:4:2, 1:3:2) = 0
  return
endif

! get position of particle in frame of multipole

if (tilt == 0) then
  x = x_pos
  y = y_pos
else
  sin_ang = sin(tilt)
  cos_ang = cos(tilt)
  x =  x_pos * cos_ang + y_pos * sin_ang
  y = -x_pos * sin_ang + y_pos * cos_ang
endif

! compute kick matrix

mat = 0

do m = 0, n, 2
  mat(1,1) = mat(1,1) +  &
                  knl * (n-m) * c_multi(n, m) * mexp(x, n-m-1) * mexp(y, m)
  mat(1,2) = mat(1,2) +  &
                  knl * m * c_multi(n, m) * mexp(x, n-m) * mexp (y, m-1)
enddo

do m = 1, n, 2
  mat(2,1) = mat(2,1) +  &
                  knl * (n-m) * c_multi(n, m) * mexp(x, n-m-1) * mexp(y, m)
  mat(2,2) = mat(2,2) +  &
                  knl * m * c_multi(n, m) * mexp(x, n-m) * mexp(y, m-1)
enddo

! transform back to lab frame

if (tilt /= 0) then
  rot(1,1) =  cos_ang
  rot(1,2) = -sin_ang
  rot(2,1) =  sin_ang
  rot(2,2) =  cos_ang
  mat = matmul(rot, mat)
  rot(1,2) =  sin_ang
  rot(2,1) = -sin_ang
  mat = matmul (mat, rot)
endif

kick_mat(2,1) = mat(1,1)
kick_mat(2,3) = mat(1,2)
kick_mat(4,1) = mat(2,1)
kick_mat(4,3) = mat(2,2)

end subroutine mat4_multipole

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine bbi_slice_calc (ele, n_slice, z_slice)
!
! Routine to compute the longitudinal positions of the slices of
! a beambeam element.
!
! Input:
!   ele        -- ele_struct: beambeam element
!   n_slice    -- integer: Number of slices
!
! Output:
!   z_slice(:) -- real(rp): Array of slice positions 
!-

subroutine bbi_slice_calc (ele, n_slice, z_slice)

implicit none

type (ele_struct) ele

integer :: i, n_slice
real(rp) z_slice(:), y
real(rp) :: z_norm

!

if (n_slice == 1) then
  z_slice(1) = ele%value(z_offset$)
elseif (n_slice > 1) then
  do i = 1, n_slice
    y = (i - 0.5) / n_slice - 0.5
    z_norm = inverse(probability_funct, y, -5.0_rp, 5.0_rp, 1.0e-5_rp)
    z_slice(i) = ele%value(sig_z$) * z_norm + ele%value(z_offset$)
  enddo
else
  print *, 'ERROR IN BBI_SLICE_CALC: N_SLICE IS NEGATIVE:', n_slice
  if (global_com%exit_on_error) call err_exit
endif

z_slice(n_slice+1) = 0

end subroutine bbi_slice_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+      
! Subroutine tilt_mat6 (mat6, tilt)
!
! Subroutine to transform a 6x6 transfer matrix to a new reference frame
! that is tilted in (x, Px, y, Py) with respect to the old reference frame.
!     mat6 -> tilt_mat * mat6 * tilt_mat_inverse
!
! Modules needed:
!   use bmad
!
! Input:
!   mat6(6,6) -- Real(rp): Untilted matrix.
!   tilt      -- Real(rp): Tilt angle.
!
! Output:
!   mat6(6,6) -- Real(rp): Tilted matrix.
!-

subroutine tilt_mat6 (mat6, tilt)

implicit none

real(rp) tilt, mat6(6,6), mm(6,6)
real(rp) c, s

!

if (tilt == 0) return

c = cos(tilt)
s = sin(tilt)

mm(1,:) = c * mat6(1,:) - s * mat6(3,:)
mm(2,:) = c * mat6(2,:) - s * mat6(4,:)
mm(3,:) = c * mat6(3,:) + s * mat6(1,:)
mm(4,:) = c * mat6(4,:) + s * mat6(2,:)
mm(5,:) =     mat6(5,:)
mm(6,:) =     mat6(6,:)

mat6(:,1) = mm(:,1) * c - mm(:,3) * s
mat6(:,2) = mm(:,2) * c - mm(:,4) * s
mat6(:,3) = mm(:,3) * c + mm(:,1) * s
mat6(:,4) = mm(:,4) * c + mm(:,2) * s
mat6(:,5) = mm(:,5)
mat6(:,6) = mm(:,6)

end subroutine tilt_mat6

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine drift_mat6_calc (mat6, length, ele, param, start, end)
!
! Subroutine to calculate a drift transfer matrix with a possible kick.
!
! Modules needed:
!   use bmad
!
! Input:
!  length   -- Real(rp): Drift length. Can be different from ele%value(l$).
!  ele      -- Real(rp): Element to drift through.
!  start    -- coord_struct: Starting coords
!  end      -- coord_struct, optional: Ending coords. Only needed if there is a kick.
!
! Output:
!   mat6(6,6) -- Real(rp): Transfer matrix
!-

subroutine drift_mat6_calc (mat6, length, ele, param, start, end)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) start
type (coord_struct), optional :: end

real(rp) ave(6), e_tot
real(rp) mat6(:,:), length, rel_pc, px, py, pxy2, pz, rel_len

!

call mat_make_unit(mat6)

if (length == 0) return

if (present(end)) then
  ave = (start%vec + end%vec) / 2
else
  ave = start%vec
endif

if (start%species == photon$) then
  rel_pc = start%p0c / ele%value(p0c$)
  px = ave(2)
  py = ave(4)
  pz = ave(6)
  pxy2 = px**2 + py**2
  e_tot = rel_pc
else
  rel_pc = 1 + ave(6)
  px = ave(2) / rel_pc
  py = ave(4) / rel_pc
  pxy2 = px**2 + py**2
  if (pxy2 > 1) return  ! prevent bomb
  pz = sqrt(1 - pxy2)
  e_tot = ele%value(p0c$) * (1 + ave(6)) / start%beta
endif

rel_len = length / (rel_pc * pz)

mat6(1,2) =  rel_len * (px**2 / pz**2 + 1)
mat6(3,4) =  rel_len * (py**2 / pz**2 + 1)
mat6(1,4) =  rel_len * px*py / pz**2
mat6(3,2) =  rel_len * px*py / pz**2
mat6(1,6) = -rel_len * px / pz**2
mat6(3,6) = -rel_len * py / pz**2
mat6(5,2) = -rel_len * px / pz**2 
mat6(5,4) = -rel_len * py / pz**2
mat6(5,6) =  rel_len * (px**2 + py**2) / pz**2 + &
                  length * mass_of(start%species)**2 * ele%value(e_tot$) / e_tot**3

end subroutine drift_mat6_calc

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_coord_transformation (mat6, orbit, param, dr, w_mat)
!
! Subroutine to calculate the matrix associated with a coordinate transformation.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele        -- real(rp): Element to drift through.
!   orbit      -- coord_struct: Starting coords
!   param      -- lat_param_struct:
!   dr(3)      -- real(rp): Transformation offset.
!   w_mat(3,3) -- real(rp): Transformation rotation matrix. 
!
! Output:
!   mat6(6,6) -- Real(rp): Transfer matrix
!-

subroutine mat6_coord_transformation (mat6, ele, param, orbit, dr, w_mat)

implicit none

type (ele_struct) ele
type (lat_param_struct) param
type (coord_struct) orbit


real(rp) dr(3), w_mat(3,3), mat6(6,6), m6_drift(6,6), rel_len
real(rp) rel_p, w_inv(3,3), beta_ref, vec(3), pvec(3), px, py, pz

!

w_inv = transpose(w_mat)

rel_p = 1 + orbit%vec(6)
pz = sqrt(rel_p**2 - orbit%vec(2)**2 - orbit%vec(4)**2)
beta_ref = ele%value(p0c$) / ele%value(e_tot$)

mat6(1,:) = [w_inv(1,1), 0.0_rp, w_inv(1,2), 0.0_rp, 0.0_rp, 0.0_rp]
mat6(3,:) = [w_inv(2,1), 0.0_rp, w_inv(2,2), 0.0_rp, 0.0_rp, 0.0_rp]

mat6(2,:) = [0.0_rp, (w_inv(1,1) - w_inv(1,3) * orbit%vec(2) / pz), &
             0.0_rp, (w_inv(1,2) - w_inv(1,3) * orbit%vec(4) / pz), 0.0_rp, w_inv(1,3) * rel_p / pz]
mat6(4,:) = [0.0_rp, (w_inv(2,1) - w_inv(2,3) * orbit%vec(2) / pz), &
             0.0_rp, (w_inv(2,2) - w_inv(2,3) * orbit%vec(4) / pz), 0.0_rp, w_inv(2,3) * rel_p / pz]

mat6(5,:) = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp, 0.0_rp]
mat6(6,:) = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]

! Drift

call mat_make_unit (m6_drift)

vec  = matmul (w_inv, [orbit%vec(1), orbit%vec(3), 0.0_rp] - dr)
pvec = matmul (w_inv, [orbit%vec(2), orbit%vec(4), pz])
px = pvec(1); py = pvec(2); pz = pvec(3)
rel_len = -vec(3) / pz

m6_drift(1,2) = rel_len * (px**2 / pz**2 + rel_p**2)
m6_drift(3,4) = rel_len * (py**2 / pz**2 + rel_p**2)
m6_drift(1,4) = rel_len * px*py / pz**2
m6_drift(3,2) = rel_len * px*py / pz**2
m6_drift(1,6) = - rel_len * rel_p * px / pz**2
m6_drift(3,6) = - rel_len * rel_p * py / pz**2
m6_drift(5,2) = - rel_len * rel_p * px / pz**2
m6_drift(5,4) = - rel_len * rel_p * py / pz**2

! Reference particle drift length is zero.
m6_drift(5,6) = -vec(3) * (px**2 + py**2) / pz**3 

! These matrix terms are due to the variation of vec(3) drift length
m6_drift(1,1) = m6_drift(1,1) + (w_inv(1,3) / w_inv(3,3)) * px / pz
m6_drift(1,3) = m6_drift(1,3) + (w_inv(2,3) / w_inv(3,3)) * px / pz

m6_drift(3,1) = m6_drift(3,1) + (w_inv(1,3) / w_inv(3,3)) * py / pz
m6_drift(3,3) = m6_drift(3,3) + (w_inv(2,3) / w_inv(3,3)) * py / pz

m6_drift(5,1) = m6_drift(5,1) - (w_inv(1,3) / w_inv(3,3)) * rel_p / pz
m6_drift(5,3) = m6_drift(5,3) - (w_inv(2,3) / w_inv(3,3)) * rel_p / pz

mat6 = matmul (m6_drift, mat6)

end subroutine mat6_coord_transformation

end module
