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
! Subroutine to check if an orbit is outside the aperture.
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
!               The exception is when the orbit is larger than 
!               bmad_com%max_aperture_limit. In this case param%lost will
!               be set to True.
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

  real(rp) x_lim, y_lim, x_beam, y_beam, l2

!

  x_lim = ele%value(x_limit$)
  if (x_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    x_lim = bmad_com%max_aperture_limit

  y_lim = ele%value(y_limit$)
  if (y_lim <= 0 .or. .not. param%aperture_limit_on) &
                                    y_lim = bmad_com%max_aperture_limit

  if (x_lim == 0 .and. y_lim == 0) return

  l2 = ele%value(l$) / 2
  x_beam = orb%vec(1) - ele%value(x_offset_tot$)  ! - ele%value(x_pitch$) * l2
  y_beam = orb%vec(3) - ele%value(y_offset_tot$)  ! - ele%value(y_pitch$) * l2

  if (ele%key == ecollimator$) then
    if (x_lim == 0 .or. y_lim == 0) then
      print *, 'ERROR IN CHECK_APERTURE_LIMIT: ECOLLIMATOR HAS ONE LIMIT ZERO'
      print *, '      AND THE OTHER NOT: ', ele%name
      call err_exit
    endif
    if ((x_beam / x_lim)**2 + (y_beam / y_lim)**2 > 1) param%lost = .true.
  else
    if (abs(x_beam) > x_lim) param%lost = .true.
    if (abs(y_beam) > y_lim) param%lost = .true.
  endif

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
!   orb(6) -- Real(rp): Orbit at start of the drift
!   length -- Real(rp): Length of drift.
!
! Output:
!   orb(6) -- Real(rp): Orbit at end of the drift
!-

subroutine track_a_drift (orb, length)

  implicit none

  real(rp) orb(6), length, rel_E

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
! If the k1 quadrupole component is zero then the tracking through the body
! is an exact (gometrical) calculation.
! For e1 or e2 non-zero this subroutine treats the edges as thin quads.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position with x', y' (not cononical).
!   ele    -- Ele_struct: Element
!
! Output:
!   end     -- Coord_struct: End position with x', y' (not cononical).
!   is_lost -- Logical: Set T or F depending upon whether the particle 
!              reaches the exit face.
!-

subroutine track_a_bend (start, ele, param, end)

  implicit none

  type (coord_struct)  start, end, start2
  type (ele_struct)  ele
  type (param_struct) param

  real(8) g0, g, r, r0, theta0, del, x1, xp1, zp, x_center, y_center
  real(8) cos0, sin0, xc, ys, b, c, x2, x_exit, y_exit, theta, s_travel
  real(8) cos1, sin1, radix, f, length

  real(rp) k1, kc, mat2(2,2), phi, mat_i6(6), dE, fact

! some init

  g0 = ele%value(g$) 
  g =  (ele%value(g$) + ele%value(delta_g$)) / (1 + start%vec(6))
  length = ele%value(l$)

  end = start
  call offset_particle (ele, param, end, set$)

  if (g == 0) then
    end%vec(1) = end%vec(1) + length * end%vec(2)
    end%vec(3) = end%vec(3) + length * end%vec(4)
    end%vec(5) = end%vec(5) - length * (end%vec(2)**2 + end%vec(4)**2) / 2 
    call offset_particle (ele, param, end, unset$)
    return
  endif

  r0 = 1 / g0
  r  = 1 / g
 
  theta0 = length * g0

! Track through the entrence face. Treat as thin lens.
! The second order terms come from the Hamiltonian term:
!       H = (g * sec^2(e1$) / 2) * p_x * y^2
! (See the MAD8 Physics writeup)

  del = tan(ele%value(e1$)) * g
  end%vec(2) = end%vec(2) + del * end%vec(1)
  end%vec(4) = end%vec(4) - del * end%vec(3)

  f = g / (2 * cos(ele%value(e1$))**2)

  end%vec(1) = end%vec(1) + f * end%vec(3)**2
  end%vec(4) = end%vec(4) - 2 * f * end%vec(2) * end%vec(3)

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

    x1  = end%vec(1)
    xp1 = end%vec(2)
    zp = end%vec(6)

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
      print *, 'ERROR IN TRACK_A_BEND: TRAJECTORY DOES NOT INTERSECT FACE.'
      print *, '      [THAT IS, THE PARTICLE AMPLITUDE IS TOO LARGE.]'
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

    end%vec(1) = x2
    end%vec(2) = tan(atan(xp1) + theta0 - theta)
    end%vec(3) = end%vec(3) + end%vec(4) * s_travel
    end%vec(5) = end%vec(5) + length - s_travel * sqrt(1 + end%vec(4)**2) 


! k1 /= 0

  case default

    dE = start%vec(6)
    k1 = k1 / (1 + dE)
    kc = g**2 + k1
    start2 = end  ! Save coords after entrence face

    call quad_mat2_calc (-kc, length, mat2)
    end%vec(1:2) = matmul (mat2, end%vec(1:2))

    call quad_mat2_calc (k1, length, mat2)
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

  f = g / (2 * cos(ele%value(e2$))**2)

  end%vec(1) = end%vec(1) - f * end%vec(3)**2
  end%vec(4) = end%vec(4) + 2 * f * end%vec(2) * end%vec(3)

  del = tan(ele%value(e2$)) * g
  end%vec(2) = end%vec(2) + del * end%vec(1)
  end%vec(4) = end%vec(4) - del * end%vec(3)

  call offset_particle (ele, param, end, unset$)

end subroutine

end module
