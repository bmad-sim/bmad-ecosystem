!+
! Subroutine track_bend (start, ele, param, end)
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
!Revision 1.8  2002/08/05 22:01:33  dcs
!Bug fix to k1 /= 0 stuff
!
!Revision 1.7  2002/08/05 20:04:59  dcs
!Bug fix for sbends calling make_mat6.
!
!Revision 1.6  2002/07/16 20:44:02  dcs
!*** empty log message ***
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

#include "CESR_platform.inc"

subroutine track_bend (start, ele, param, end)

  use bmad_struct
  use bmad_interface
  use make_mat6_mod

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
      type *, 'ERROR IN TRACK_BEND: TRAJECTORY DOES NOT INTERSECT FACE.'
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
    end%z%pos = end%z%pos + ele%value(l$) - s_travel 


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
