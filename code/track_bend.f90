!+
! Subroutine track_bend (start, ele, end, is_lost)
!
! Particle tracking through a bend element.
! Subroutine assumes no k1 quadrupole component
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
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

subroutine track_bend (start, ele, end, is_lost)

  use bmad_struct
  implicit none

  type (coord_struct)  start, end
  type (ele_struct)  ele
                    
  real*8 r0, r, theta0, del, x1, xp1, zp, x_center, y_center
  real*8 cos0, sin0, xc, ys, b, c, x2, x, y, theta, s_travel
  real*8 cos1, sin1, radix

  logical is_lost

! some init

  r0 = ele%value(rho$) 
  r =  r0 * (1 + start%z%vel) 
  theta0 = ele%value(angle$)

! track through the entrence face. Treat as thin lens.

  end = start

  del = tan(ele%value(e1$)) / r
  end%x%vel = end%x%vel + del * end%x%pos
  end%y%vel = end%y%vel - (del+ end%x%vel/r) * end%y%pos

! Track through main body...
! We use a local coordinate system (x, y) aligned with
! the entrence face so local x is the same as the particle x and 
! local y is the same as the particle z.
! the local coordinate system origin is the nominal center of rotation

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
    is_lost = .true.
    return
  else
    is_lost = .false.
  endif

  x2 = (-b + sqrt(radix)) / 2

  x = (r0 + x2) * cos0 - x_center
  y = (r0 + x2) * sin0 - y_center

  theta = atan2 (y, x) + atan(xp1)
  s_travel = r * theta

  end%x%pos = x2
  end%x%vel = tan(atan(xp1) + theta0 - theta)
  end%y%pos = end%y%pos + end%y%vel * s_travel
  end%z%pos = end%z%pos + ele%value(l$) - s_travel 

! Track through the exit face. Treat as thin lens.

  del = tan(ele%value(e2$)) / r
  end%x%vel = end%x%vel + del * end%x%pos
  end%y%vel = end%y%vel - (del-end%x%vel/r) * end%y%pos

end subroutine
