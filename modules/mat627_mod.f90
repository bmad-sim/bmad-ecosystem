#include "CESR_platform.inc"

!+
!--------------------------------------
! ******** THIS IS ANTIQUATED. ********
! ********     DO NOT USE.     ********
!--------------------------------------
!-

module mat627_mod

  use bmad_struct
  use bmad_interface
  use make_mat6_mod

! For 6x27 matrices

  integer, parameter :: x11$ = 7, x12$ = 8, x13$ = 9, x14$ = 10, x15$ = 11
  integer, parameter :: x16$ = 12, x22$ = 13, x23$ = 14, x24$ = 15
  integer, parameter :: x25$ = 16, x26$ = 17, x33$ = 18, x34$ = 19, x35$ = 20
  integer, parameter :: x36$ = 21, x44$ = 22, x45$ = 23, x46$ = 24
  integer, parameter :: x55$ = 25, x56$ = 26, x66$ = 27

  type mat627_struct
    real(rp) m(6,27)
  end type

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine track_long (ring, orbit_, ix_start, direction, mats627)
!
! Subroutine to track for 1-turn. This subroutine is ment for long term
! tracking and uses 2nd order transport matrices for some of the tracking.
!
! Moudules Needed:
!   use mat627_mod
!
! Input:
!   ring             -- Ring_struct: Ring.
!     %param%aperture_limit_on -- Logical: Sets whether TRACK_LONG looks to
!                                   see whether a particle is lost or not
!   orbit_(ix_start) -- Coord_struct: Coordinates at start of tracking.
!   ix_start         -- Integer: Start index (See Note).
!   direction        -- Integer: Direction to track.
!                            = +1  -> Track forward
!                            = -1  -> Track backward
!   mats627(0:n_ele_maxx) -- Mat627_struct: array of 6x27 2nd order matrices.
!                       Use RING_MAKE_MAT627 to calculate these matrices 
!                       before calling TRACK_LONG.
!
! Output:
!   ring
!     %param%lost    -- Logical: Set when a particle is lost with the aperture
!                        limit on.
!     %param%ix_lost -- Integer: set to index of element where particle is lost
!   orbit_(0:n_ele_maxx) -- Coord_struct: Coordinates of tracked particle. 
!   orbit_(ix_start) -- Coord_struct: Coordinates at end of tracking.
!
! Note: Starting and ending points are just after the elements with index
!   IX_START. For example, if DIRECTION = +1 then the first element
!   tracked through is element ix_start+1. If DIRECTION = -1 then the first
!   element tracked through is element ix_start.
!
! Note: Remember that if you are tracking backward and you have calculated
!   the orbit by tracking forward (with, for example, CLOSED_ORBIT_AT_START)
!   then you need to reverse velocity signs:
!           orbit_(ix_start)%x_vel = -orbit_(ix_start)%x_vel 
!           orbit_(ix_start)%y_vel = -orbit_(ix_start)%y_vel 
!
! Note: If x_limit (or y_limit) for an element is zero then TRACK_LONG will
!   take x_limit (or y_limit) as infinite (this is standard BMAD).
!-

subroutine track_long (ring, orbit_, ix_start, direction, mats627)

  implicit none

  type (ring_struct) ring
  type (coord_struct) orbit_(0:*)
  type (ele_struct) ele            
  type (mat627_struct) mats627(*)

  integer ix_start, direction
  integer n, i

  real(rdef) x_lim, y_lim

  logical debug / .false. /
                            
! track through elements.

  print *, 'ERROR IN TRACK_LONG: THIS SUBROUTINE IS OBSOLETE. PLEASE SEE DCS.'
  stop

  ring%param%lost = .false.

  if (direction == +1) then
    call track_fwd (ix_start+1, ring%n_ele_ring)
    if (ring%param%lost) return
    orbit_(0) = orbit_(ring%n_ele_ring) 
    call track_fwd (1, ix_start)
  elseif (direction == -1) then
    call track_back (ix_start, 1)
    if (ring%param%lost) return
    orbit_(ring%n_ele_ring) = orbit_(0)
    call track_back (ring%n_ele_ring, ix_start+1)
  else
    print *, 'ERROR IN TRACK_LONG: BAD DIRECTION:', direction
    call err_exit
  endif

contains

!--------------------------------------------------------------------------

subroutine track_fwd (ix1, ix2)

  integer ix1, ix2

  do n = ix1, ix2

    call track1_627 (orbit_(n-1), ring%ele_(n), ring%param, &
                                                 mats627(n)%m, orbit_(n))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n
      return
    endif

    if (debug) then
      print *, ring%ele_(n)%name
      print *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine

!--------------------------------------------------------------------------

subroutine track_back (ix1, ix2)

  integer ix1, ix2
  logical reverse
  real(rdef) mat_save(6,6)

!

  do n = ix1, ix2, -1

    call track1_627 (orbit_(n), ring%ele_(n), ring%param, mats627(n)%m, &
                                                                orbit_(n-1))

! check for lost particles

    if (ring%param%lost) then
      ring%param%ix_lost = n - 1
      return
    endif

    if (debug) then
      print *, ring%ele_(n)%name
      print *, (orbit_(n)%vec(i), i = 1, 6)
    endif

  enddo

end subroutine

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine RING_MAKE_MAT627 (RING, IX_ELE, DIRECTION, MATS627)
!
! Subroutine to make the 6x27 2nd order matrices for long term tracking. 
! Used by, for example, TRACK_LONG.
!
! Modules Needed:
!   use mat627_mod
!
! Input:
!   RING        -- Ring_struct: Ring containing the element.
!   IX_ELE      -- Integer: Index of the element. if < 0 then entire
!                    ring will be made. In this case group elements will
!                    be made up last.
!   DIRECTION   -- Integer: Transport Direction. 
!                   = +1  -> For Forward (normal) particle tracking.
!                   = -1  -> For Backward particle tracking.
!
! Output:
!   MATS627(n_ele_maxx) -- Mat627_struct: Array of 6x27 matrices.
!     %M(6,27)            -- 6x27 matrix.
!-

recursive subroutine ring_make_mat627 (ring, ix_ele, direction, mats627)

  implicit none

  type (ring_struct)  ring
  type (ele_struct)  ele
  type (mat627_struct) mats627(:)

  integer direction
  integer i, ix_ele, i1, i2, i3, ix1, ix2, ix3

! make entire ring if ix_ele < 0

  if (ix_ele < 0) then

    do i = ring%n_ele_ring+1, ring%n_ele_max
      if (ring%ele_(i)%control_type /= group_lord$)  &
                                 call control_bookkeeper (ring, i)
    enddo

    do i = ring%n_ele_ring+1, ring%n_ele_max
      if (ring%ele_(i)%control_type == group_lord$)  &
                                 call control_bookkeeper (ring, i)
    enddo

    do i = 1, ring%n_ele_use
      if (ring%ele_(i)%key /= hybrid$)  &
         call make_mat627(ring%ele_(i), ring%param, direction, mats627(i)%m)
    enddo

    return

  endif

!-----------------------------------------------------------
! otherwise make a single element

  ele = ring%ele_(ix_ele)
  call control_bookkeeper (ring, ix_ele)

! for a regular element

  if (ele%key == hybrid$) return

  if (ix_ele <= ring%n_ele_ring) then
    call make_mat627(ring%ele_(ix_ele), ring%param, direction, &
                                                       mats627(ix_ele)%m)
  endif                        

! for a control element

  do i1 = ring%ele_(ix_ele)%ix1_slave, ring%ele_(ix_ele)%ix2_slave
    i = ring%control_(i1)%ix_slave
    call ring_make_mat627 (ring, i, direction, mats627)
  enddo

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine MAKE_MAT627 (ELE, PARAM, DIRECTION, MAT627)
!
! Subroutine to make the 6x27 2nd order transfer matrix for an element. 
! Note: tilts are not inclueded in the calculation
!
! Modules Needed:
!   use mat627_mod
!
! Input:
!   ELE       -- Ele_struct: Element
!   PARAM     -- Param_struct: Ring parameters
!   DIRECTION -- Integer: Transport Direction. 
!                   = +1  -> For Forward (normal) particle tracking.
!                   = -1  -> For Backward particle tracking.
!
! Output:
!   MAT627(6,27) -- Real(rp): 6x27 transfer matrix.
!
! Elements for whom 2nd order terms are NOT computed:
!     SBEND, BEAMBEAM, etc...
!-

subroutine make_mat627 (ele, param, direction, mat627)

  implicit none

  type (ele_struct), target :: ele
  type (param_struct)  param

  real(rp) mat627(6,27) 
  real(rp) mat6_end(6,6), mat2(2,2), mat4(4,4), kmat1(4,4), kmat2(4,4)
  real(rp) e1, e2, angle, rho, cos_angle, sin_angle, k1, ks, length, kc
  real(rp) phi, k2l, k3l, c2, s2, cs, ks2, del_l
  real(rp) s_pos, s_pos_old, z_slice(100)
  real(rp) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)
  real(rp) r, c_e, c_m, gamma_old, gamma_new, vec_st(4)
  real(rp) sqrt_k, arg, kick2
  real(rp) cx, sx, cy, sy, k2l_2, k2l_3, k2l_4
  real(rp) x_off, y_off, x_pit, y_pit
  real(rp) lcs, lc2s2

  integer i, n, n_slice, direction

!--------------------------------------------------------
! init

  length = ele%value(l$)  
  mat627 = 0
  forall (i = 1:6) mat627(i,i) = 1
                 
!--------------------------------------------------------
! marker

  if (ele%key == marker$) return
  if (ele%key == sbend$) return

!--------------------------------------------------------
! drift or element is off or
! Electric Separator or Kicker.

  if (ele%key == drift$ .or. .not. ele%is_on .or. &
      ele%key == elseparator$ .or. ele%key == kicker$) then
    mat627(1,2) = length
    mat627(3,4) = length
    mat627(1,x26$) = -length
    mat627(3,x46$) = -length
    mat627(5,x22$) = -length / 2
    mat627(5,x44$) = -length / 2
    return
  endif

!--------------------------------------------------------

  select case (ele%key)

!--------------------------------------------------------
! quadrupole

  case (quadrupole$) 

    k1 = ele%value(k1$)

    if (k1 == 0) then
      mat627(1,2)  =  length
      mat627(3,4)  =  length
      mat627(5,x22$)  =  -length / 2.0
      mat627(5,x44$)  =  -length / 2.0
      return
    endif

    call quad_mat_calc (-k1, length, mat627(1:2,1:2))
    call quad_mat_calc ( k1, length, mat627(3:4,3:4))

    cx = mat627(1, 1)
    sx = mat627(1, 2)
    cy = mat627(3, 3)
    sy = mat627(3, 4)

    mat627(1, x16$) = k1 * length * sx / 2
    mat627(1, x26$) = -(sx + length * cx) / 2
    mat627(2, x16$) = -k1 * (sx - length * cx) / 2
    mat627(2, x26$) = mat627(1, x16$)

    mat627(3, x36$) = -k1 * length * sy / 2
    mat627(3, x46$) = -(sy + length * cy) / 2
    mat627(4, x36$) =  k1 * (sy - length * cy) / 2
    mat627(4, x46$) = mat627(3, x36$)

    mat627(5, x11$) =  k1 * (cx*sx - length) / 4 
    mat627(5, x12$) =  k1 * sx**2 / 2            
    mat627(5, x22$) = -(length + sx*cx) / 4      
  
    mat627(5, x33$) =  k1 * (length - sy*cy) / 4 
    mat627(5, x34$) = -k1 * sy**2 / 2            
    mat627(5, x44$) = -(length + sy*cy) / 4      
                     
!--------------------------------------------------------
! Sextupole.

  case (sextupole$) 

    k2l = ele%value(k2$) * length
    k2l_2 = k2l * length
    k2l_3 = k2l * length**2
    k2l_4 = k2l * length**3

    mat627(1, 2) = length
    mat627(3, 4) = length
  
    mat627(1, x11$) = -k2l_2 / 4
    mat627(1, x12$) = -k2l_3 / 6
    mat627(1, x22$) = -k2l_4 / 24
    mat627(1, x33$) =  k2l_2 / 4
    mat627(1, x34$) =  k2l_3 / 6
    mat627(1, x44$) =  k2l_4 / 24

    mat627(2, x11$) = -k2l   / 2
    mat627(2, x12$) = -k2l_2 / 2
    mat627(2, x22$) = -k2l_3 / 6
    mat627(2, x33$) =  k2l   / 2
    mat627(2, x34$) =  k2l_2 / 2
    mat627(2, x44$) =  k2l_3 / 6

    mat627(3, x13$) = k2l_2 / 2
    mat627(3, x14$) = k2l_3 / 6
    mat627(3, x23$) = k2l_3 / 6
    mat627(3, x24$) = k2l_4 / 12

    mat627(4, x13$) = k2l     
    mat627(4, x14$) = k2l_2 / 2
    mat627(4, x23$) = k2l_2 / 2
    mat627(4, x24$) = k2l_3 / 3

    mat627(1,x26$) = -length
    mat627(3,x46$) = -length
    mat627(5,x22$) = -length / 2
    mat627(5,x44$) = -length / 2

!--------------------------------------------------------
! octupole

  case (octupole$) 

    mat627(1,2) = length
    mat627(3,4) = length
    mat627(1,x26$) = -length
    mat627(3,x46$) = -length
    mat627(5,x22$) = -length / 2
    mat627(5,x44$) = -length / 2

!--------------------------------------------------------
! solenoid

  case (solenoid$) 

    ks = ele%value(ks$) * direction

    call solenoid_mat_calc (ks, length, mat627(1:4,1:4))

    c2 = mat627(1,1)
    s2 = mat627(1,4) * ks / 2
    cs = mat627(1,3)

    lcs = length * cs
    lc2s2 = length * (c2 - s2) / 2
           
    arg = length / 2         

    mat627(1, x16$) =  lcs * ks  
    mat627(1, x26$) = -lc2s2 * 2 
    mat627(1, x36$) = -lc2s2 * ks
    mat627(1, x46$) = -lcs * 2   
                      
    mat627(2, x16$) =  lc2s2 * ks**2 / 2 
    mat627(2, x26$) =  lcs * ks          
    mat627(2, x36$) =  lcs * ks**2 / 2   
    mat627(2, x46$) = -lc2s2 * ks        
                                         
    mat627(3, x16$) =  lc2s2 * ks        
    mat627(3, x26$) =  lcs * 2           
    mat627(3, x36$) =  lcs * ks          
    mat627(3, x46$) = -lc2s2 * 2         
                                         
    mat627(4, x16$) = -lcs * ks**2 / 2   
    mat627(4, x26$) =  lc2s2 * ks        
    mat627(4, x36$) =  mat627(2, x16$)             
    mat627(4, x46$) =  lcs * ks          
                      
    mat627(5, x11$) = -arg * (ks/2)**2
    mat627(5, x14$) =  arg * ks       
    mat627(5, x22$) = -arg            
    mat627(5, x23$) = -arg * ks       
    mat627(5, x33$) = -arg * (ks/2)**2
    mat627(5, x44$) = -arg            
                     
!--------------------------------------------------------
! rf cavity, etc. 
! For an RF cavity this is not quite correct.

  case (rfcavity$, beambeam$, multipole$, ab_multipole$, custom$) 

    return

!--------------------------------------------------------
! wiggler

  case (wiggler$) 

    call mat_make_unit (mat627(1:6,1:6))     ! make a unit matrix

    if (param%beam_energy == 0) then
      k1 = 0
    else
      k1 = -0.5 * (c_light * ele%value(b_max$) / param%beam_energy)**2
    endif

    mat627(1, 1) = 1
    mat627(1, 2) = length
    mat627(2, 1) = 0
    mat627(2, 2) = 1

    call quad_mat_calc (k1, length, mat627(3:4,3:4))

    cy = mat627(3, 3) 
    sy = mat627(3, 4) 

    mat627(1, x26$) = -length 

    mat627(3, x36$) = -k1 * length * sy
    mat627(3, x46$) = -length * cy
    mat627(4, x36$) = -k1 * length * cy
    mat627(4, x46$) = mat627(3, x36$)

    mat627(5, x22$) = - length / 2 

    mat627(5, x33$) =  k1 * (length - sy*cy) / 4
    mat627(5, x34$) = -k1 * sy**2 / 2
    mat627(5, x44$) = -(length + sy*cy) / 4
                     
!--------------------------------------------------------
! solenoid/quad

  case (sol_quad$) 

    k1 = ele%value(k1$)
    ks = ele%value(ks$) * direction

    call sol_quad_mat627_calc (ks, k1, length, mat627)

!--------------------------------------------------------
! rbends are not allowed internally

  case (rbend$) 

    print *, 'ERROR IN MAKE_MAT627: RBEND ELEMENTS NOT ALLOWED INTERNALLY!'
    call err_exit

!--------------------------------------------------------
! unrecognized element

  case default

    print *, 'ERROR IN MAKE_MAT627: UNKNOWN ELEMENT KEY',  &
                                    ele%key, '  ', key_name(ele%key)

  end select

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine SOL_QUAD_MAT627_CALC (KS, K1, LENGTH, MAT627)
!
! Subroutine to calculate the transfer matrix for a combination 
! solenoid/quadrupole element.
!
! Modules Needed:
!   use mat627_mod
!
! Input:
!   ks      [Real]       Solenoid strength
!   k1      [Real]       Quadrupole strength
!   length  [Real]       Sol_quad length
!
! Output
!   mat627(6,27) [Real]  Transfer matrix across the sol_quad
!-

subroutine sol_quad_mat627_calc (ks, k1, s_len, m)

  implicit none

  real(rp) ks, k1, s_len
  real(rp) m(6,27)

  integer i, j

  real(rp) ks2, s, c, snh, csh
  real(rp) darg1, alpha, alpha2, beta, beta2, f, q, r, a, b
  real(rp) df, dalpha2, dalpha, dbeta2, dbeta, darg
  real(rp) dC, dCsh, dS, dSnh, dq, dr, da, db
  real(rp) ks3, fp, fm, dfm, dfp, df_f, ug
  real(rp) s1, s2, snh1, snh2, dsnh1, dsnh2, ds1, ds2
  real(rp) coef1, coef2, dcoef1, dcoef2, ks4

  real(rp) t5(4,4), t6(4,4)

! Calc
          
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
  Snh2 = Snh / beta

  coef1 = ks2*r + 4*k1*a
  coef2 = ks2*q + 4*k1*b

  m = 0
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

  df      = -2 * (ks4 + 2*k1**2) / f
  dalpha2 = df/2 - ks2
  dalpha  = (df/2 - ks2)/(2*alpha)
  dbeta2  = ks2 + df/2
  dbeta   = (ks2 + df/2)/(2*beta)
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
  dSnh2 = dSnh / beta - Snh * dbeta / beta2

  dcoef1 = -2*ks2*r + ks2*dr - 4*k1*a + 4*k1*da
  dcoef2 = -2*ks2*q + ks2*dq - 4*k1*b + 4*k1*db                     

  t6(1,1) = m(1,1)*df_f + 2*ug*(fp*dC + C*dfp + fm*dCsh + Csh*dfm)
  t6(1,2) = m(1,2)*df_f + (2*ug/k1) * (dq*S1 + q*dS1 - dr*Snh1 - r*dSnh1)
  t6(1,3) = m(1,3)*df_f + (ks*ug/k1)*(-db*S1 - b*dS1 + da*Snh1 + a*dSnh1)
  t6(1,4) = m(1,4)*(df_f - 2) + 4*ks*ug*(-dC + dCsh) 

  t6(2,1) = m(2,1)*(df_f + 1) - &
              (ug/2)*(dcoef1*S2 + coef1*dS2 + dcoef2*Snh2 + coef2*dSnh2)
  t6(2,2) = t6(1,1)
  t6(2,3) = m(2,3)*(df_f - 2) + ks3*ug*(dC - dCsh) 
  t6(2,4) = m(2,4)*(df_f - 1) + ug*ks*(da*S2 + a*dS2 + db*Snh2 + b*dSnh2)

  t6(3,1) = -t6(2,4)
  t6(3,2) = -t6(1,4)
  t6(3,3) = m(3,3)*df_f + 2*ug*(fm*dC + C*dfm + fp*dCsh + Csh*dfp)
  t6(3,4) = m(3,4)*(df_f - 1) + 2*ug*(dr*S2 + r*dS2 + dq*Snh2 + q*dSnh2)

  t6(4,1) = -t6(2,3)        
  t6(4,2) = -t6(1,3)
  t6(4,3) = m(4,3)*(df_f + 2) + &
               (ug/(2*k1))*(-dcoef2*S1 - coef2*dS1 + dcoef1*Snh1 + coef1*dSnh1)
  t6(4,4) = t6(3,3)

  m(1:4,x16$) = t6(1:4,1)
  m(1:4,x26$) = t6(1:4,2)
  m(1:4,x36$) = t6(1:4,3)
  m(1:4,x46$) = t6(1:4,4)

  do i = 1, 4
    do j = 1, 4
      t5(i,j) = -t6(2,j) * m(1,i) + t6(1,j) * m(2,i) - &
                     t6(4,j) * m(3,i) + t6(3,j) * m(4,i)
      if (i == j) t5(i,j) = t5(i,j) / 2
    enddo
  enddo

  m(5, (/ x11$, x12$, x13$, x14$ /) ) = t5(1,1:4)
  m(5, (/ x12$, x22$, x23$, x24$ /) ) = t5(2,1:4)
  m(5, (/ x13$, x23$, x33$, x34$ /) ) = t5(3,1:4)
  m(5, (/ x14$, x24$, x34$, x44$ /) ) = t5(4,1:4)

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine TRACK1_627 (START, ELE, PARAM, mat627, END)
!
! Particle tracking through a single element.
! Uses the 2nd order 6x27 transport matrices whenever possible.
! This is for long term tracking.
!
! Modules Needed:
!   use mat627_mod
!
! Input:
!   START        -- Coord_struct: Starting position
!   ELE          -- Ele_struct: Element
!   PARAM        -- Param_struct:
!   mat627(6,27) -- Real(rp): 6x27 2nd order transport matrix.
!
! Output:
!   END   -- Coord_struct: End position
!
! Notes:
!
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!-

subroutine track1_627 (start, ele, param, mat627, end)

  implicit none

  type (coord_struct)  start, end
  type (ele_struct)  ele
  type (param_struct)  param

  real(rp) x_kick, y_kick
  real(rp) mat627(6,27)
  real(rp) x_lim, y_lim

!-------------------------------------------------------------------
! some simple cases

  select case (ele%key)
  case (marker$, drift$, elseparator$, kicker$, beambeam$, octupole$, &
        sbend$, rfcavity$, hybrid$, multipole$, ab_multipole$, custom$, &
        wiggler$)
    call track1 (start, ele, param, end)
    return
  end select


!-------------------------------------------------------------------
! 2nd order tracking
! initially set end = start

  end = start     ! transfer start to end

  select case (ele%key)
  case (quadrupole$, sextupole$, solenoid$, sol_quad$) 
    
    call offset_particle (ele, param, end, set$, set_canonical = .false.)
    call track1_order2 (end%vec, ele, mat627, end%vec)
    call offset_particle (ele, param, end, unset$, set_canonical = .false.)

! error

  case default
    print *, 'ERROR IN TRACK1_627: UNKNOWN ELEMENT: ', &
                                       key_name(ele%key), ele%type
    call err_exit
  end select

!-------------------------------------------------------------------
! check for particles outside aperture

  if (param%aperture_limit_on) then

    x_lim = ele%value(x_limit$)
    if (x_lim <= 0) x_lim = 1e10
    if (abs(end%vec(1)) > x_lim) param%lost = .true.

    y_lim = ele%value(y_limit$)
    if (y_lim <= 0) y_lim = 1e10
    if (abs(end%vec(3)) > y_lim) param%lost = .true.

  endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine track1_order2 (vec_begin, ele, mat627, vec_end)

  implicit none

  type (ele_struct) ele

  real(rp) vec_begin(6), vec_end(6), vv(27), mat627(6,27)

!

  vv(1:6)   = vec_begin
  vv(7:12)  = vv(1) * vv(1:6)
  vv(13:17) = vv(2) * vv(2:6)
  vv(18:21) = vv(3) * vv(3:6)
  vv(22:24) = vv(4) * vv(4:6)
  vv(25:26) = vv(5) * vv(5:6)
  vv(27:27) = vv(6) * vv(6:6)

  vec_end = matmul(mat627, vv)      

end subroutine

end module
