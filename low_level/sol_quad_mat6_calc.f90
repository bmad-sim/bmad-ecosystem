!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine sol_quad_mat6_calc (ks, k1, tilt, length, ele, orbit, mat6, make_matrix)
!
! Subroutine to calculate the transfer matrix for a combination 
! solenoid/quadrupole element (without a tilt).
!
! Input:
!   ks          -- real(rp): Solenoid strength.
!   k1          -- real(rp): Quadrupole strength.
!   tilt        -- real(rp): quadrupole tilt.
!   length      -- real(rp): Sol_quad length.
!   ele         -- ele_struct: Sol_quad element.
!   orbit       -- coord_struct: Orbit at beginning of the sol_quad.
!   mat6(6,6)   -- real(rp), optional: Transfer matrix up to the sol_quad.
!   make_matrix -- logical, optional: Extend the matrix?
!
! Output:
!   orbit       -- coord_struct: Orbit at beginning of the sol_quad.
!   mat6(6,6)   -- real(rp): Transfer matrix includeing the sol_quad.
!-

subroutine sol_quad_mat6_calc (ks_in, k1_in, tilt, length, ele, orbit, mat6, make_matrix)

use equal_mod, dummy => sol_quad_mat6_calc

implicit none

type (coord_struct) orbit
type (ele_struct) ele

real(rp) ks, k1, length
real(rp), optional :: mat6(6,6)

real(rp) ks2, s, c, snh, csh, rel_p, ks_in, k1_in, tilt, orb(6)
real(rp) darg1, alpha, alpha2, beta, beta2, f, q, r, a, b
real(rp) df, dalpha2, dalpha, dbeta2, dbeta, darg, kmat(6,6)
real(rp) dC, dCsh, dS, dSnh, dq, dr, da, db, k1_2
real(rp) ks3, fp, fm, dfm, dfp, df_f, ug, e_tot
real(rp) s1, s2, snh1, snh2, dsnh1, dsnh2, ds1, ds2
real(rp) coef1, coef2, dcoef1, dcoef2, ks4
real(rp) t4(4,4), ts(4,4), m0(6,6), xp_start, xp_end, yp_start, yp_end
real(rp) dt4(4,4), dts(4,4), dz_coef(4,4), d_dz_co(4,4), dm(4,4)
real(rp) r_orb(4), d_orb(4), tsd(4,4), dtsd(4,4)
real(rp) d2_f, dug, d2_ug, d2_C, d2_fp, d2_fm, d2_Csh, d2_q, d2_s1, d2_s2 
real(rp) d2_r, d2_snh1, d2_snh2, d2_a, d2_b, d2_coef1, d2_coef2
real(rp) d2_alpha, d2_beta, d2_arg, d2_arg1, d2_S, d2_Snh, d2_f_f, rk

logical, optional :: make_matrix

! Calculation is done in (x, x', y, y') coordinates and then converted
! to (x, p_x, y, p_y) coordinates.

call tilt_coords (tilt, orbit%vec, mat6, make_matrix)

rel_p = 1 + orbit%vec(6)

orb = orbit%vec
orb(2) = orbit%vec(2) / rel_p
orb(4) = orbit%vec(4) / rel_p

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

S = sin(alpha*length)                              
C = cos(alpha*length)
Snh = sinh(beta*length)
Csh = cosh(beta*length)
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
  Snh2 = length
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
darg    = length*dalpha
darg1   = length*dbeta         
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
  dSnh2 = k1**4 * length**3 * (-1/3.0d0 + (40 - ks2 * length**2) * k1_2 / (30 * ks4)) / ks3**2
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

m0(1:4,6) = matmul(t4(1:4,1:4), orb(1:4))

! kmat(5,6) calc

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
d2_arg  = length*d2_alpha
d2_arg1 = length*d2_beta
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
  d2_Snh2 = 4*k1**4*length**3/(3*ks**6)-2*k1**4*length**3/ks**6
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

! The kmat(5,6) term is computed 

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

r_orb = [orb(1), orb(2)*rel_p, orb(3), orb(4)*rel_p]
d_orb = [0.0_rp, -orb(2)/rel_p, 0.0_rp, -orb(4)/rel_p]

! dz = Sum_ij dz_coef(i,j) * orb(i) * orb(j)

dz_coef = matmul (ts, m0(1:4,1:4)) / 2
dz_coef(:,2) = dz_coef(:,2) / rel_p 
dz_coef(:,4) = dz_coef(:,4) / rel_p 
dz_coef(2,:) = dz_coef(2,:) / rel_p 
dz_coef(4,:) = dz_coef(4,:) / rel_p 

! Energy corrections

kmat = m0

kmat(1,2) = m0(1,2) / rel_p
kmat(1,4) = m0(1,4) / rel_p

kmat(2,1) = m0(2,1) * rel_p
kmat(2,3) = m0(2,3) * rel_p

kmat(3,2) = m0(3,2) / rel_p
kmat(3,4) = m0(3,4) / rel_p

kmat(4,1) = m0(4,1) * rel_p
kmat(4,3) = m0(4,3) * rel_p

kmat(1,6) = m0(1,6) / rel_p
kmat(3,6) = m0(3,6) / rel_p

! Propagate orbit

orbit%vec(5) = orbit%vec(5) + sum(orbit%vec(1:4) * matmul(dz_coef, orbit%vec(1:4))) + low_energy_z_correction (orbit, ele, length)
orbit%vec(1:4) = matmul (kmat(1:4,1:4), orbit%vec(1:4))

! mat6 calc

if (logic_option(.false., make_matrix)) then
  dm = t4 / rel_p
  dm(2,1) = t4(2,1)
  dm(2,3) = t4(2,3)
  dm(4,1) = t4(4,1)
  dm(4,3) = t4(4,3)
  dm(1,2) = t4(1,2) / rel_p**2
  dm(1,4) = t4(1,4) / rel_p**2
  dm(3,2) = t4(3,2) / rel_p**2
  dm(3,4) = t4(3,4) / rel_p**2

  kmat(5,6) = (dot_product(matmul(matmul(d_orb, tsd), kmat(1:4,1:4)), r_orb) + &
               dot_product(matmul(matmul(orb(1:4), dtsd), kmat(1:4,1:4)), r_orb) + &
               dot_product(matmul(matmul(orb(1:4), tsd), dm(1:4,1:4)), r_orb)) / 2

  ! The m(5,x) terms follow from the symplectic condition.

  kmat(5,1) = -kmat(2,6)*kmat(1,1) + kmat(1,6)*kmat(2,1) - kmat(4,6)*kmat(3,1) + kmat(3,6)*kmat(4,1)
  kmat(5,2) = -kmat(2,6)*kmat(1,2) + kmat(1,6)*kmat(2,2) - kmat(4,6)*kmat(3,2) + kmat(3,6)*kmat(4,2)
  kmat(5,3) = -kmat(2,6)*kmat(1,3) + kmat(1,6)*kmat(2,3) - kmat(4,6)*kmat(3,3) + kmat(3,6)*kmat(4,3)
  kmat(5,4) = -kmat(2,6)*kmat(1,4) + kmat(1,6)*kmat(2,4) - kmat(4,6)*kmat(3,4) + kmat(3,6)*kmat(4,4)

  ! 1/gamma^2 m56 correction
  e_tot = orbit%p0c * (1 + orbit%vec(6)) / orbit%beta
  kmat(5,6) = kmat(5,6) + length * (mass_of(orbit%species))**2 * ele%value(e_tot$) / e_tot**3

  mat6 = matmul(kmat, mat6)
endif

call tilt_coords (-tilt, orbit%vec, mat6, make_matrix)

end subroutine sol_quad_mat6_calc

