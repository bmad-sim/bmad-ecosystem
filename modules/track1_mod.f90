!+
! Module track1_mod
!
! Module of helper routines for track1 routines
!-

#include "CESR_platform.inc"

module track1_mod

  use bmad_struct
  use bmad_interface
  use make_mat6_mod  

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine check_aperture_limit (orb, ele, param)
!
! Subroutine to check if an orbit is outside the aperture
!
! Modules needed:
!   use bmad
!
! Input:
!   orb -- Coord_struct: coordinates of a particle.
!   ele -- Ele_struct: Element holding the aperture
!     %value(x_limit$) -- Horizontal aperture.
!     %value(y_limit$) -- Vertical aparture.
!   param -- Param_struct: Parameter structure
!     %aperture_limit_on -- The aperture limit is only checked if this is true.
!
! Output:
!   param -- Param_struct: Parameter structure:
!     %lost -- Set True if the orbit is outside the aperture. 
!              Note: %lost is NOT set False if the orbit is inside 
!                the aperture.
!-

subroutine check_aperture_limit (orb, ele, param)

  implicit none

  type (coord_struct), intent(in) :: orb
  type (ele_struct),   intent(in)  :: ele
  type (param_struct), intent(inout) :: param

  real(rdef) x_lim, y_lim

!

  if (.not. param%aperture_limit_on) return

  x_lim = ele%value(x_limit$)
  if (x_lim <= 0) x_lim = 1e10
  if (abs(orb%x%pos) > x_lim) param%lost = .true.

  y_lim = ele%value(y_limit$)
  if (y_lim <= 0) y_lim = 1e10
  if (abs(orb%y%pos) > y_lim) param%lost = .true.

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_drift (orb, length)
!
! Subroutine to track through a drift.
!
! Modules needed:
!   use precision_def
!
! Input:
!   orb(6) -- Real(rdef): Orbit at start of the drift
!   length -- Real(rdef): Length of drift.
!
! Output:
!   orb(6) -- Real(rdef): Orbit at end of the drift
!-

subroutine track_a_drift (orb, length)

  implicit none

  real(rdef) orb(6), length, rel_E

  rel_E = 1 + orb(6)

  orb(1) = orb(1) + length * orb(2) / rel_E
  orb(3) = orb(3) + length * orb(4) / rel_E
  orb(5) = orb(5) - length * (orb(2)**2 + orb(4)**2) / (2 * rel_E**2)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine track_a_bend (start, ele, param, end)
!
! Particle tracking through a bend element.
! This subroutine assumes no k1 quadrupole component.
! For e1 or e2 non-zero this subroutine treats the dipole edges as thin quads.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position with x', y' (not cononical)
!   ele    -- Ele_struct: Element
!
! Output:
!   end     -- Coord_struct: End position with x', y' (not cononical)
!   is_lost -- Logical: Set T or F depending upon whether the particle 
!              reaches the exit face.
!-

!$Id$
!$Log$
!Revision 1.3  2002/11/05 18:06:29  dcs
!Minor Mods
!
!Revision 1.2  2002/09/05 14:01:07  dcs
!*** empty log message ***
!
!Revision 1.1  2002/08/20 20:37:06  dcs
!Added symp_lie_bmad / symp_lie_ptc switches.
!
!Revision 1.8  2002/08/05 22:01:33  dcs
!Bug fix to k1 /= 0 stuff
!
!Revision 1.7  2002/08/05 20:04:59  dcs
!Bug fix for sbends calling make_mat6.
!
!Revision 1.5  2002/06/13 14:54:29  dcs
!Interfaced with FPP/PTC
!
!Revision 1.4  2002/02/23 20:32:26  dcs
!Double/Single Real toggle added
!
!Revision 1.3  2002/01/08 21:44:44  dcs
!Aligned with VMS version  -- DCS
!
!Revision 1.2  2001/09/27 18:31:59  rwh24
!UNIX compatibility updates
!

subroutine track_a_bend (start, ele, param, end)

  implicit none

  type (coord_struct)  start, end, start2
  type (ele_struct)  ele
  type (param_struct) param
                    
  real*8 g0, g, r, r0, theta0, del, x1, xp1, zp, x_center, y_center
  real*8 cos0, sin0, xc, ys, b, c, x2, x_exit, y_exit, theta, s_travel
  real*8 cos1, sin1, radix

  real(rdef) k1, kc, mat2(2,2), phi, mat_i6(6), dE, fact, length

! some init

  g0 = ele%value(g_design$) 
  g =  ele%value(g$) / (1 + start%z%vel)

  end = start
  call offset_particle (ele, param, end, set$)

  if (g == 0) then
    end%vec(1) = end%vec(1) + ele%value(l$) * end%vec(2)
    end%vec(3) = end%vec(3) + ele%value(l$) * end%vec(4)
    call offset_particle (ele, param, end, unset$)
    return
  endif

  r0 = 1/ g0
  r = 1/ g
 
  theta0 = ele%value(l$) * g0

! track through the entrence face. Treat as thin lens.

  del = tan(ele%value(e1$)) * g
  end%x%vel = end%x%vel + del * end%x%pos
  end%y%vel = end%y%vel - (del+end%x%vel*g) * end%y%pos

!-----------------------------------------------------------------------
! Track through main body...
! Two cases:
!  If k1 == 0 then do an exact calculation using simple geometry.
!  If k1 /= 0 then essentually just use the 1st order transfer matrix.

  k1 = ele%value(k1$)
  select case (k1 == 0)

  case (.true.)  ! k1 == 0

! We use a local coordinate system (x, y) aligned with
! the entrence face so local x is the same as the particle x and 
! local y is the same as the particle z.
! the local coordinate system origin is the nominal center of rotation
! (that is the center for an on-energy particle with zero offsets).

! For reverse bends with g and theta negative then the calculation is done 
! effectively under the transformation: 
!           g -> -g,  theta -> -theta,  x -> -x,  Px -> -Px

! x,y_center is the center of the actual rotation

    x1  = end%x%pos
    xp1 = end%x%vel
    zp = end%z%vel

    x_center =  (r0 + x1) - r / sqrt(1 + xp1**2)
    y_center = r * xp1 / sqrt(1 + xp1**2)   

    cos0 = cos(theta0)
    sin0 = sin(theta0)

    xc = x_center * cos0
    ys = y_center * sin0

    b = 2 * (r0 - xc - ys) 
    c = x_center**2 + y_center**2 - 2 * r0 * (xc + ys) + (r0**2 - r**2)

    radix = b**2 - 4*c
    if (radix < 0) then
      type *, 'ERROR IN TRACK_A_BEND: TRAJECTORY DOES NOT INTERSECT FACE.'
      type *, '      [THAT IS, THE PARTICLE AMPLITUDE IS TOO LARGE.]'
      param%lost = .true.
      return
    else
      param%lost = .false.
    endif

! x_exit, y_exit is the point where the particle intersects the exit face.
! x2 is the distance from the nominal exit point to the actual exit point.

    x2 = (-b + sign(sqrt(radix), r0)) / 2 

    x_exit = (r0 + x2) * cos0 - x_center
    y_exit = (r0 + x2) * sin0 - y_center

    if (g0 > 0) then
      theta = atan2 (y_exit, x_exit) + atan(xp1)
    else                      ! for reverse bends
      theta = atan2 (y_exit, x_exit) - pi + atan(xp1)
    endif

    s_travel = r * theta

    end%x%pos = x2
    end%x%vel = tan(atan(xp1) + theta0 - theta)
    end%y%pos = end%y%pos + end%y%vel * s_travel
    end%z%pos = end%z%pos + ele%value(l$) - s_travel * sqrt(1 + end%y%vel**2) 


! k1 /= 0

  case default

    dE = start%z%vel
    k1 = k1 / (1 + dE)
    kc = g**2 + k1
    length = ele%value(l$)
    start2 = end  ! Save coords after entrence face

    call quad_mat_calc (-kc, length, mat2)
    end%vec(1:2) = matmul (mat2, end%vec(1:2))

    call quad_mat_calc (k1, length, mat2)
    end%vec(3:4) = matmul (mat2, end%vec(3:4))

    phi = sqrt(abs(kc)) * length
    if (kc < 0) then
      mat_i6(1) = (1 - cosh(phi)) * g / kc
      mat_i6(2) = sinh(phi) * g / sqrt(-kc)
      mat_i6(5) = (phi - sinh(phi)) * g**2 / abs(kc)**1.5
    else
      mat_i6(1) = (1 - cos(phi)) * g / kc
      mat_i6(2) = sin(phi) * g / sqrt(kc)
      mat_i6(5) = (sin(phi) - phi) * g**2 / kc**1.5
    endif

    end%vec(1) = end%vec(1) + mat_i6(1) * dE
    end%vec(2) = end%vec(2) + mat_i6(2) * dE
    end%vec(5) = end%vec(5) + mat_i6(5) * dE

    fact = start2%vec(4) * (mat_i6(5) + length)

    end%vec(5) = end%vec(5) - &
                    mat_i6(2) * start2%vec(1) - mat_i6(1) * start2%vec(2) + &
                    fact * start2%vec(4) / 2
    end%vec(3) = end%vec(3) + fact * dE

  end select
  
!------------------------------------------------------------------------
! Track through the exit face. Treat as thin lens.

  del = tan(ele%value(e2$)) * g
  end%x%vel = end%x%vel + del * end%x%pos
  end%y%vel = end%y%vel - (del-end%x%vel*g) * end%y%pos

  call offset_particle (ele, param, end, unset$)

end subroutine

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

!+
! Subroutine track_a_accel_sol (start, ele, param, end)
!
! Subroutine to track through an accel_sol element
!
! Modules Needed:
!   use bmad
!
! Input:
!   START  -- Coord_struct: Starting position
!   ELE    -- Ele_struct: Element
!   PARAM  -- Param_struct:
!     %APERTURE_LIMIT_ON -- If .true. then %LOST will be set if the
!                 particle is outsile the aperture.
!
! Output:
!   END   -- Coord_struct: End position
!   PARAM
!     %LOST -- Set .true. If the particle is outside the aperture and
!                %APERTURE_LIMIT_ON is set. Also: %LOST is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %APERTURE_LIMIT_ON.
!-
                
subroutine track_a_accel_sol (start, ele, param, end) 

  implicit none

  type (coord_struct)  start, end, c0
  type (ele_struct)  ele
  type (param_struct)  param

  real(rdef) gamma_b, gamma_new, gamma_old, l_over_gamma, ll, ls(5), s_cumul
  real(rdef) s_grand_cum, vec_st(4), x_beg_lim, y_beg_lim, x_lim, y_lim
  real(rdef) x_lim_chng_rate, y_lim_chng_rate, phase, mat4(4,4), length
  real(rdef) b_x(5), b_y(5), beta_b, beta_s, c_e, c_m, en_gain, gam_inv2_b

  real(rdef), parameter :: beta_crit$ = 0.999

  integer i, j

  logical check_x_ap, check_y_ap, inside_segment

! linac stuff
! beta_b is the total speed in units of c_light (before entering the element)
! beta_s is the longitudinal speed in units of c_light

  gamma_b = param%energy * (end%z%vel + 1) / e_mass
  beta_b = sqrt(1 - 1 / gamma_b**2)
  gam_inv2_b = 1.0 / gamma_b**2
  if (gam_inv2_b <= 0.001) then
    beta_s = 1/sqrt(1 + end%x%vel**2 + end%y%vel**2) *(1-gam_inv2_b/2)
  else
    beta_s = sqrt((1 - gam_inv2_b) / (1 + end%x%vel**2 + end%y%vel**2))
  endif

! Calculation of end.z.vel

  end = start
  length = ele%value(l$)

  call offset_particle (ele, param, end, set$)

  if (ele%value(volt$) /= 0) then
    phase = ele%value(lag$) + end%z%pos  &
                      / (ele%value(rf_wavelength$) * beta_s)
    en_gain = ele%value(volt$) * sin(twopi * phase) / 1.e9
    if ((en_gain + gamma_b * e_mass) <= e_mass) then
      param%lost = .true.
      return
    else
      end%z%vel = end%z%vel + en_gain / param%energy
      c_e = en_gain / (e_mass * length)
    endif
  else
    c_e = 0.0
  endif

! Beginning fringe

  c_m = param%particle * c_light * ele%value(b_z$) / (e_mass * 1.e9)
  call mat_unit(mat4, 4, 4)
  mat4(2,3) = c_m / 2 *  &
                  sqrt((1 + end%x%vel**2 + end%y%vel**2) / (gamma_b**2 - 1))
  mat4(4,1) = -mat4(2,3)
  end%vec(1:4) = matmul(mat4, end%vec(1:4))

! Segment before first steerings:

  ls(1) = ele%value(s_st1$)
  b_x(1) = 0
  b_y(1) = 0

! Segment with first steerings:

  ls(2) = ele%value(l_st1$)
  b_x(2) = ele%value(b_x1$)
  b_y(2) = ele%value(b_y1$)

! Segment between steerings:

  ls(3) = ele%value(s_st2$) - (ele%value(s_st1$) + ele%value(l_st1$))
  b_x(3) = 0
  b_y(3) = 0

! Segment with second steerings:

  ls(4)= ele%value(l_st2$)
  b_x(4) = ele%value(b_x2$)
  b_y(4) = ele%value(b_y2$)

! Segment after second steerings:

  ls(5) = length - (ele%value(s_st2$) + ele%value(l_st2$))
  b_x(5) = 0
  b_y(5) = 0

  gamma_old = gamma_b
  gamma_new = gamma_old

  s_grand_cum = 0
  x_beg_lim = ele%value(x_beg_limit$)
  y_beg_lim = ele%value(y_beg_limit$)
  x_lim_chng_rate = (ele%value(x_limit$) - x_beg_lim) / length
  y_lim_chng_rate = (ele%value(y_limit$) - y_beg_lim) / length
  check_x_ap = .false.
  check_y_ap = .false.
  if (abs(c_m) > 0.001) then
    l_over_gamma = pi / (10 * abs(c_m))
    if (ele%value(x_limit$) * x_beg_lim /= 0) check_x_ap = .true.
    if (ele%value(y_limit$) * y_beg_lim /= 0) check_y_ap = .true.
  else
    l_over_gamma = 100 * pi
  endif

  do j = 1, 5

    s_cumul = 0
    inside_segment = .true.
    do while (inside_segment)
      gamma_old = gamma_new
      ll = l_over_gamma * gamma_old
      if (s_cumul + ll >= ls(j)) then
        ll = ls(j) - s_cumul
        inside_segment = .false.
      endif
      s_cumul = s_cumul + ll
      s_grand_cum = s_grand_cum + ll

      if (ll > 1.e-5) then
        gamma_new = gamma_old + c_e * ll
        call accel_sol_mat_calc (ll, c_m, c_e, gamma_old, gamma_new,  &
                                    b_x(j), b_y(j), end%vec, mat4, vec_st)
        end%vec(1:4) = matmul(mat4, end%vec(1:4))
        do i = 1, 4
          end%vec(i) = end%vec(i) + param%particle * vec_st(i)
        enddo

! Calculation of end.z.pos

        if (abs(c_e) > 0.001) then
          end%z%pos = end%z%pos - ll + (sqrt(1 + (c_e * ll  &
                          + sqrt(gamma_old**2 - 1))**2) - gamma_old) / c_e
        else
          end%z%pos = end%z%pos + ll * (sqrt(1 - 1/gamma_old**2) - 1  &
              - c_e * ll * (gamma_old**2 - 0.5)/gamma_old**3)
        endif

        if (param%aperture_limit_on) then
          if (check_x_ap) then
            x_lim = x_lim_chng_rate * s_grand_cum + x_beg_lim
            if (abs(end%x%pos) > x_lim) param%lost = .true.
          endif
          if (check_y_ap) then
            y_lim = y_lim_chng_rate * s_grand_cum + y_beg_lim
            if (abs(end%y%pos) > y_lim) param%lost = .true.
          endif
          if (param%lost) return
        endif

      endif
    enddo

  enddo

! Ending fringe

  call mat_unit(mat4, 4, 4)
  mat4(4,1) = c_m / 2 *  &
                sqrt((1 + end%x%vel**2 + end%y%vel**2) / (gamma_new**2 - 1))
  mat4(2,3) = -mat4(4,1)
  end%vec(1:4) = matmul(mat4, end%vec(1:4))

  call offset_particle (ele, param, end, unset$)

end subroutine

end module
