#include "CESR_platform.inc"

module make_mat6_mod

  use dcslib
  use physical_constants
  use bmad_basic_mod

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine quad_mat2_calc (k1, length, mat2, dz_coef)
!
! Subroutine to calculate the 2x2 transfer matrix for a quad for
! one plane. 
! Note: This calc does not include any energy corrections for an
! off-energy particle.
!
! Modules needed:
!   use bmad
!
! Input:
!   k1     -- Real(rp): Quad strength: k1 > 0 ==> defocus.
!   length -- Real(rp): Quad length
!
! Output:
!   mat2(2,2)  -- Real(rp): Transfer matrix.
!   dz_coef(3) -- Real(rp), optional: Coefficients for calculating the
!                 the change in z position:
!                   dz = Integral [-x'^2/2 ds]
!                      = c(1) * x_0^2 + c(2) * x_0 * x'_0 + c(3) * x'_0^2 
!-

subroutine quad_mat2_calc (k1, length, mat2, dz_coef)

  implicit none

  real(rp) length, mat2(2,2), cx, sx
  real(rp) k1, sqrt_k, arg, arg2
  real(rp), optional :: dz_coef(3)

!

  sqrt_k = sqrt(abs(k1))
  arg = sqrt_k * length

  if (arg < 1e-10) then
    arg2 = k1 * length**2
    cx = 1 - arg2 / 2
    sx = (1 - arg2 / 6) * length
  elseif (k1 < 0) then       ! focus
    cx = cos(arg)
    sx = sin(arg) / sqrt_k
  else                       ! defocus
    cx = cosh(arg)
    sx = sinh(arg) / sqrt_k
  endif

  mat2(1,1) = cx
  mat2(1,2) = sx
  mat2(2,1) = k1 * sx
  mat2(2,2) = cx

  if (present(dz_coef)) then
    dz_coef(1) = k1 * (-cx * sx + length) / 4
    dz_coef(2) = -k1 * sx**2 / 2
    dz_coef(3) = -(cx * sx + length) / 4
  endif

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine sol_quad_mat6_calc (ks, k1, length, mat6, orb, dz_coef)
!
! Subroutine to calculate the transfer matrix for a  combination 
! solenoid/quadrupole element (without a tilt).
!
! Note: This routine is not meant to be for general use.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ks      -- Real(rp): Solenoid strength.
!   k1      -- Real(rp): Quadrupole strength.
!   length  -- Real(rp): Sol_quad length.
!   orb(6)  -- Real(rp): Orbit at beginning of the sol_quad.
!              Note: using non-cononical coords here!
!
! Output:
!   mat6(6,6)    -- Real(rp): Transfer matrix across the sol_quad.
!   dz_coef(4,4) -- Real(rp), optional: coefs for dz calc
!                     dz = sum_ij dz_coef(i,j) * orb(i) * orb(j)
!-

subroutine sol_quad_mat6_calc (ks_in, k1_in, s_len, m, orb, dz_coef)

  implicit none

  real(rp) ks, k1, s_len
  real(rp) m(6,6)
  real(rp) orb(6)

  real(rp) ks2, s, c, snh, csh, rel_p, ks_in, k1_in
  real(rp) darg1, alpha, alpha2, beta, beta2, f, q, r, a, b
  real(rp) df, dalpha2, dalpha, dbeta2, dbeta, darg
  real(rp) dC, dCsh, dS, dSnh, dq, dr, da, db
  real(rp) ks3, fp, fm, dfm, dfp, df_f, ug
  real(rp) s1, s2, snh1, snh2, dsnh1, dsnh2, ds1, ds2
  real(rp) coef1, coef2, dcoef1, dcoef2, ks4
  real(rp) t4(4,4), ts(4,4), xp_start, xp_end, yp_start, yp_end
  real(rp), optional :: dz_coef(4,4)

! Calculation is done in (x, x', y, y') coordinates and then converted
! to (x, p_x, y, p_y) coordinates.

  rel_p = 1 + orb(6)
  k1 = k1_in / rel_p
  ks = ks_in / rel_p

  ks2 = ks**2
  ks3 = ks2 * ks 
  ks4 = ks2*ks2
  f = sqrt(ks4 + 4*k1**2)
  ug = 1 / (4*f)
  alpha2 = (f + ks2) / 2; alpha = sqrt(alpha2)
  beta2  = (f - ks2) / 2; beta  = sqrt(beta2)
  S = sin(alpha*s_len)                              
  C = cos(alpha*s_len)
  Snh = sinh(beta*s_len)
  Csh = cosh(beta*s_len)
  q = f + 2*k1 - ks2
  r = f - 2*k1 + ks2
  a = f + 2*k1 + ks2
  b = f - 2*k1 - ks2
  fp = f + 2*k1
  fm = f - 2*k1

  S1 = S * alpha
  S2 = S / alpha

  Snh1 = Snh * beta

  if (abs(beta) < 1e-10) then
    Snh2 = s_len
  else
    Snh2 = Snh / beta
  endif

  coef1 = ks2*r + 4*k1*a
  coef2 = ks2*q + 4*k1*b

  call mat_make_unit(m)
               
  m(1,1) = 2*ug * (fp*C + fm*Csh)
  m(1,2) = (2*ug/k1) * (q*S1 - r*Snh1)
  m(1,3) = (ks*ug/k1) * (-b*S1 + a*Snh1)
  m(1,4) = 4*ug*ks * (-C + Csh)

  m(2,1) = -(ug/2) * (coef1*S2 + coef2*Snh2)
  m(2,2) = m(1,1)             
  m(2,3) = ug*ks3 * (C - Csh)
  m(2,4) = ug*ks * (a*S2 + b*Snh2)

  m(3,1) = -m(2,4)
  m(3,2) = -m(1,4)
  m(3,3) = 2*ug * (fm*C + fp*Csh)  
  m(3,4) = 2*ug * (r*S2 + q*Snh2)

  m(4,1) = -m(2,3)     
  m(4,2) = -m(1,3)
  m(4,3) = (ug/(2*k1)) * (-coef2*S1 + coef1*Snh1)
  m(4,4) = m(3,3)

! 

  if (all(orb == 0) .and. .not. present(dz_coef)) return

  df      = -2 * (ks4 + 2*k1**2) / f
  dalpha2 = df/2 - ks2
  dalpha  = (df/2 - ks2)/(2*alpha)
  dbeta2  = ks2 + df/2
  if (beta < 1e-4) then
    dbeta   = -abs(k1**3/(ks3*ks2))
  else
    dbeta   = (ks2 + df/2)/(2*beta)
  endif
  darg    = s_len*dalpha
  darg1   = s_len*dbeta         
  dC      = -darg*S
  dCsh    = darg1*Snh
  dS      = darg*C
  dSnh    = darg1*Csh
  dq      = -2*k1 + 2*ks2 + df
  dr      =  2*k1 - 2*ks2 + df
  da      = -2*k1 - 2*ks2 + df
  db      =  2*k1 + 2*ks2 + df
  dfp = df - 2*k1
  dfm = df + 2*k1
  df_f =  -df/f

  dS1 = dS * alpha + S * dalpha
  dS2 = dS / alpha - S * dalpha / alpha2

  dSnh1 = dSnh * beta + Snh * dbeta

  if (beta < 1e-4) then
    dSnh2 = -k1**4 * s_len**3 / (3 * ks3**2)
  else
    dSnh2 = dSnh / beta - Snh * dbeta / beta2
  endif

  dcoef1 = -2*ks2*r + ks2*dr - 4*k1*a + 4*k1*da
  dcoef2 = -2*ks2*q + ks2*dq - 4*k1*b + 4*k1*db                     

! t4(i,j) is dm(i,j)/dE and give the m(x,6) terms.

  t4(1,1) = m(1,1)*df_f + 2*ug*(fp*dC + C*dfp + fm*dCsh + Csh*dfm)
  t4(1,2) = m(1,2)*df_f + (2*ug/k1) * (dq*S1 + q*dS1 - dr*Snh1 - r*dSnh1)
  t4(1,3) = m(1,3)*df_f + (ks*ug/k1)*(-db*S1 - b*dS1 + da*Snh1 + a*dSnh1)
  t4(1,4) = m(1,4)*(df_f - 2) + 4*ks*ug*(-dC + dCsh) 

  t4(2,1) = m(2,1)*(df_f + 1) - &
              (ug/2)*(dcoef1*S2 + coef1*dS2 + dcoef2*Snh2 + coef2*dSnh2)
  t4(2,2) = t4(1,1)
  t4(2,3) = m(2,3)*(df_f - 2) + ks3*ug*(dC - dCsh) 
  t4(2,4) = m(2,4)*(df_f - 1) + ug*ks*(da*S2 + a*dS2 + db*Snh2 + b*dSnh2)

  t4(3,1) = -t4(2,4)
  t4(3,2) = -t4(1,4)
  t4(3,3) = m(3,3)*df_f + 2*ug*(fm*dC + C*dfm + fp*dCsh + Csh*dfp)
  t4(3,4) = m(3,4)*(df_f - 1) + 2*ug*(dr*S2 + r*dS2 + dq*Snh2 + q*dSnh2)

  t4(4,1) = -t4(2,3)        
  t4(4,2) = -t4(1,3)
  t4(4,3) = m(4,3)*(df_f + 2) + &
               (ug/(2*k1))*(-dcoef2*S1 - coef2*dS1 + dcoef1*Snh1 + coef1*dSnh1)
  t4(4,4) = t4(3,3)

  m(1:4,6) = matmul(t4(1:4,1:4), orb(1:4))

! dz = Sum_ij dz_coef(i,j) * orb(i) * orb(j)

  if (present(dz_coef)) then
    ts(1:4,1) = -t4(2,1:4)
    ts(1:4,2) =  t4(1,1:4)
    ts(1:4,3) = -t4(4,1:4)
    ts(1:4,4) =  t4(3,1:4)
    dz_coef = matmul (ts, m(1:4,1:4)) / 2
  endif

! energy corrections

  m(1,2) = m(1,2) / rel_p
  m(1,4) = m(1,4) / rel_p

  m(2,1) = m(2,1) * rel_p
  m(2,3) = m(2,3) * rel_p

  m(3,2) = m(3,2) / rel_p
  m(3,4) = m(3,4) / rel_p

  m(4,1) = m(4,1) * rel_p
  m(4,3) = m(4,3) * rel_p

  m(1,6) = m(1,6) / rel_p
  m(3,6) = m(3,6) / rel_p

! The m(5,x) terms follow from the symplectic condition.

  m(5,1) = -m(2,6)*m(1,1) + m(1,6)*m(2,1) - m(4,6)*m(3,1) + m(3,6)*m(4,1)
  m(5,2) = -m(2,6)*m(1,2) + m(1,6)*m(2,2) - m(4,6)*m(3,2) + m(3,6)*m(4,2)
  m(5,3) = -m(2,6)*m(1,3) + m(1,6)*m(2,3) - m(4,6)*m(3,3) + m(3,6)*m(4,3)
  m(5,4) = -m(2,6)*m(1,4) + m(1,6)*m(2,4) - m(4,6)*m(3,4) + m(3,6)*m(4,4)

! The m(5,6) term is computed 

  xp_start = orb(2) + ks_in * orb(3) / 2
  yp_start = orb(4) - ks_in * orb(1) / 2
  xp_end = sum(orb(1:4) * m(2,1:4)) + ks_in * orb(3) / 2
  yp_end = sum(orb(1:4) * m(4,1:4)) - ks_in * orb(1) / 2
  m(5,6) = s_len * (xp_start**2 + yp_start**2 +xp_end**2 + yp_end**2) / (2 * rel_p**3)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine sbend_body_with_k1_map (g, g_err, length, k_1, start, end, mat6)
!
! Subroutine to calculate the transfer matrix and/or ending coordinates for a 
! sbend with a finite k1 but without a tilt.
!
! Modules Needed:
!   use bmad
!
! Input:
!   g        -- Real(rp): 1/rho. rho = radius of curvature of reference orbit.
!   g_err    -- Real(rp): Actual bend field proportional to (g + g_err)
!   length   -- Real(rp): Length of the bend.
!   k_1      -- Real(rp): focusing strength.
!   start(6) -- Real(rp): Orbit at beginning of the bend.
!
! Output:
!   mat6(6,6)  -- Real(rp), optional: Transfer matrix across the sol_quad.
!   end(6)     -- Real(rp), optional: Ending coordinates.
!-

subroutine sbend_body_with_k1_map (g, g_err, length, k_1, start, end, mat6)

  implicit none

  real(rp) g, g_err, length, start(6)
  real(rp), optional :: end(6), mat6(6,6)
  real(rp) k_1, k_x, x_c, om_x, om_y, tau_x, tau_y, arg, s_x, c_x, s_y, c_y, r(6)
  real(rp) z0, z1, z2, z11, z12, z22, z33, z34, z44
  real(rp) dom_x, dom_xx, dx_c, dc_x, ds_x, dom_y, dom_yy, dc_y, ds_y, dcs_x, dcs_y
  real(rp) g_tot, rel_p, rel_p2

!

  g_tot = g + g_err
  rel_p = (1 + start(6))
  rel_p2 = rel_p**2
  

  k_x = k_1 + g * g_tot
  x_c = (g * rel_p - g_tot) / k_x

  om_x = sqrt(abs(k_x) / rel_p)
  om_y = sqrt(abs(k_1) / rel_p)

  tau_x = -sign (1.0_rp, k_x)
  tau_y =  sign (1.0_rp, k_1)

  arg = om_x * length
  if (arg < 1e-6) then
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
  if (arg < 1e-6) then
    s_y = (1 + tau_y * arg**2 / 6) * length
    c_y = 1 + tau_y * arg**2 / 2
  elseif (k_1 < 0) then
    s_y = sin(om_y * length) / om_y
    c_y = cos(om_y * length)
  else
    s_y = sinh(om_y * length) / om_y
    c_y = cosh(om_y * length)
  endif

  r = start
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

  if (present(end)) then
    end(1) = c_x * r(1) + s_x * r(2) / rel_p + x_c
    end(2) = tau_x * om_x**2 * rel_p * s_x * r(1) + c_x * r(2)
    end(3) = c_y * r(3) + s_y * r(4) / rel_p
    end(4) = tau_y * om_y**2 * rel_p * s_y * r(3) + c_y * r(4)
    end(5) = r(5) + z0 + z1 * r(1) + z2 * r(2) + &
                    z11 * r(1)**2 + z12 * r(1) * r(2) + z22 * r(2)**2 + &
                    z33 * r(3)**2 + z34 * r(3) * r(4) + z44 * r(4)**2 
  endif




end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat6_multipole (knl, tilt, c00, factor, mat6)
!
! Subroutine to add to a transfer matrix the effect of a multipole kick.
! This routine is not meant for general use.
!
! Input:
!   knl  -- Real(rp): Strength of multipole
!   tilt -- Real(rp): Tilt of multipole
!   c00  -- Coord_struct: coordinates of particle around which the
!             multipole kick matrix is computed.
!   mat6(6,6) -- Real(rp): Transfer matrix.
!
! Output:
!   mat6(6,6) -- Real(rp): Transfer matrix with effect of a multiplole added.
!-

subroutine mat6_multipole (knl, tilt, c00, factor, mat6)

  implicit none

  real(rp) c00(6)
  real(rp) mat6(6,6), kmat1(4,4), factor
  real(rp) knl(0:), tilt(0:)

  integer n

!                        

  if (c00(1) == 0 .and. c00(3) == 0 .and. knl(1) == 0) return

  do n = 1, ubound(knl, 1)
    if (knl(n) /= 0) then
      call mat4_multipole (knl(n), tilt(n), n, c00, kmat1)
      mat6(2:4:2, 1:3:2) = mat6(2:4:2, 1:3:2) + factor * kmat1(2:4:2, 1:3:2)
    endif
  enddo

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine mat4_multipole (knl, tilt, n, c0, kick_mat)
!
! Subroutine to find the kick matrix (Jacobian) due to a multipole.
! This routine is not meant for general use.
!
! Input:
!   c0   -- Coord_struct: coordinates of particle
!   knl  -- Real(rp): Strength of multipole
!   tilt -- Real(rp): Tilt of multipole
!
! Output:
!   kick_mat(4,4) -- Real(rp): Kick matrix (Jacobian) at c0.
!-


subroutine mat4_multipole (knl, tilt, n, c0, kick_mat)

  implicit none

  real(rp) c0(6)
  real(rp) x_pos, y_pos, x, y, knl, tilt
  real(rp) sin_ang, cos_ang, mat(2,2), rot(2,2)
  real(rp) kick_mat(4,4)

  integer m, n

! init

  kick_mat = 0
  forall (m = 1:4) kick_mat(m,m) = 1

  x_pos = c0(1)
  y_pos = c0(3)
           
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

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine bbi_slice_calc (n_slice, sig_z, z_slice)

  implicit none

  integer :: i, n_slice, n_slice_old = 0 
  real(rp) sig_z, z_slice(:), y
  real(rp), save :: z_norm(100)

!

  if (n_slice == 1) then
    z_slice(1) = 0
  elseif (n_slice > 1) then
    do i = 1, n_slice
      if (n_slice /= n_slice_old) then
        y = (i - 0.5) / n_slice - 0.5
        z_norm(i) = inverse(probability_funct, y, -5.0_rp, 5.0_rp, 1.0e-5_rp)
      endif
      z_slice(i) = sig_z * z_norm(i)
    enddo
    n_slice_old = n_slice
  else
    print *, 'ERROR IN BBI_SLICE_CALC: N_SLICE IS NEGATIVE:', n_slice
    call err_exit
  endif

  z_slice(n_slice+1) = 0

end subroutine

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

  mm(1,1:6) = c * mat6(1,1:6) - s * mat6(3,1:6)
  mm(2,1:6) = c * mat6(2,1:6) - s * mat6(4,1:6)
  mm(3,1:6) = c * mat6(3,1:6) + s * mat6(1,1:6)
  mm(4,1:6) = c * mat6(4,1:6) + s * mat6(2,1:6)
  mm(5,1:6) =     mat6(5,1:6)
  mm(6,1:6) =     mat6(6,1:6)

  mat6(1:6,1) = mm(1:6,1) * c - mm(1:6,3) * s
  mat6(1:6,2) = mm(1:6,2) * c - mm(1:6,4) * s
  mat6(1:6,3) = mm(1:6,3) * c + mm(1:6,1) * s
  mat6(1:6,4) = mm(1:6,4) * c + mm(1:6,2) * s
  mat6(1:6,5) = mm(1:6,5)
  mat6(1:6,6) = mm(1:6,6)
                     
end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine solenoid_mat_calc (ks, length, mat4)
!
! Subroutine to calculate the 4x4 transverse transfer matrix for a solenoid.
! This routine is not meant for general use.
!-

subroutine solenoid_mat_calc (ks, length, mat4)

  implicit none

  real(rp) ks, length, kss, c, s, c2, s2, cs, ll, kl, kl2
  real(rp) mat4(4,4)

!

  kss = ks / 2

  if (abs(length * kss) < 1e-10) then
    ll = length
    kl = kss * length 
    kl2 = kl**2
    mat4(1,:) = (/  1.0_rp,   ll,      kl,      kl*ll    /)
    mat4(2,:) = (/ -kl*kss,   1.0_rp, -kl2*kss, kl       /)
    mat4(3,:) = (/ -kl,      -kl*ll,   1.0_rp,  ll       /)
    mat4(4,:) = (/  kl2*kss, -kl,     -kl*kss,  1.0_rp /)
    return
  endif

  c = cos(kss*length)
  s = sin(kss*length)
  c2 = c*c
  s2 = s*s
  cs = c*s

  mat4(1,1) = c2
  mat4(1,2) = cs / kss
  mat4(1,3) = cs
  mat4(1,4) = s2 / kss
  mat4(2,1) = -kss * cs
  mat4(2,2) = c2
  mat4(2,3) = -kss * s2 
  mat4(2,4) = cs
  mat4(3,1) = -cs
  mat4(3,2) = -s2 / kss
  mat4(3,3) = c2
  mat4(3,4) = cs / kss
  mat4(4,1) = kss * s2
  mat4(4,2) = -cs
  mat4(4,3) = -kss * cs
  mat4(4,4) = c2

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine drift_mat6_calc (mat6, length, start, end)
!
! Subroutine to calculate a drift transfer matrix with a possible kick.
!-

subroutine drift_mat6_calc (mat6, length, start, end)

  implicit none

  real(rp) start(:)
  real(rp), optional :: end(:)
  real(rp) ave(6)
  real(rp) mat6(:,:), length, rel_p, len_p2, len_p3

!

  call mat_make_unit(mat6)

  if (length == 0) return

  if (present(end)) then
    ave = (start + end) / 2
  else
    ave = start
  endif

  rel_p = 1 + ave(6)
  len_p2 = length / rel_p**2
  len_p3 = len_p2 / rel_p

  mat6(1,2) = length / rel_p
  mat6(3,4) = length / rel_p
  mat6(1,6) = -ave(2) * len_p2
  mat6(3,6) = -ave(4) * len_p2
  mat6(5,2) = -ave(2) * len_p2
  mat6(5,4) = -ave(4) * len_p2
  if (present(end)) then
    mat6(5,6) = (start(2)**2 + start(2)*end(2) + end(2)**2 + &
                 start(4)**2 + start(4)*end(4) + end(4)**2) * len_p3 / 3
  else
    mat6(5,6) = (start(2)**2 + start(4)**2) * len_p3
  endif

end subroutine

end module
