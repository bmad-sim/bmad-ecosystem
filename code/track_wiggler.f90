!+
! Subroutine track_wiggler (start, ele, param, end, is_lost)
!
! Subroutine to track through the body of a wiggler. This routine is
! used by track1.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!
! Output:
!   end   -- Coord_struct: End position
!   is_lost -- Logical: Set T or F depending upon whether the particle
!              reaches the exit face.!-

! For wiggler tracking use the "hard edge" model of Wiedemann, 
! "Part. Acc. Ph. II", pg. 65. Note that in Wiedemann fig 2.8, l_h is 
! mistakenly drawn as being the full pole width (it is the half pole width).

subroutine track_wiggler (start, ele, param, end, is_lost)

  use bmad_struct
  use bmad_interface

  implicit none

  type (coord_struct)  start, end, before
  type (ele_struct)  ele, bend
  type (param_struct)  param

  real length, k_z, factor, dx
  real const1, const3, tan_theta
  real x_lim, y_lim, l_period, l_bend, l_drift, rho_bend, angle
  real knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

  integer i, j, n, n_slice, n_pole

  logical, parameter :: set$ = .true., unset$ = .false.    
  logical is_lost, pole_multipoles_on

!

  is_lost = .false.      ! assume success
  end = start     ! transfer start to end

! check to see if the wiggler is off. 
! if so then it just looks like a drift.

  if (ele%value(rho$) == 0) then
    end%vec(1) = end%vec(1) + length * end%vec(2) * (1 - end%vec(6))
    end%vec(3) = end%vec(3) + length * end%vec(4) * (1 - end%vec(6))
    return
  endif

! normal case of wiggler on.

  call offset_coords_m (ele, param, end, set$, .true., .true.)
  call tilt_coords (ele%value(tilt$), end%vec, set$)
             
  length = ele%value(l$)
  n_pole = nint(ele%value(n_pole$))

  if (n_pole == 1) then
    print *, 'ERROR IN TRACK_WIGGLER: WIGGLER HAS ONLY 1 POLE FOR: ', ele%name
    call err_exit
  endif

  l_period = length / n_pole 
  l_bend = 8 * l_period / pi**2
  l_drift = l_period - l_bend
  rho_bend = 4 * ele%value(rho$) / pi
  angle = l_bend / rho_bend

  k_z = pi / l_period

  const1 = 1 / (rho_bend * (1 + end%z%vel))
  const3 =  (pi/l_period)**2 / 6 

  bend = ele
  bend%value(ix1_m$:ix2_m$) = 0
  bend%value(a2$) = ele%value(ap2$);  bend%value(b2$) = ele%value(bp2$)
  bend%value(a4$) = ele%value(ap4$);  bend%value(b4$) = ele%value(bp4$) 
  bend%value(a6$) = ele%value(ap6$);  bend%value(b6$) = ele%value(bp6$) 
  bend%value(a8$) = ele%value(ap8$);  bend%value(b8$) = ele%value(bp8$) 
  call multipole_to_vecs(bend, param%particle, knl, tilt)
  knl = knl / 2

  if (any(knl /= 0)) then
    pole_multipoles_on = .true.
  else
    pole_multipoles_on = .false.
  endif

! If the wiggler has an odd number of poles then the end poles get their
! strength changed so that a particle entering at the origin leaves at
! the origin

  if (mod(n_pole, 2) == 1) then
    factor = sqrt(rho_bend**2 - (l_bend/2)**2)
    dx = 2 * (rho_bend - factor) + l_drift * l_bend / (2 * factor)
    factor = 1 / (1 - 2 * rho_bend * dx / (l_bend * (length-l_bend/2)))
  endif

! wiggler track

  do i = 1, n_pole 
    call track_period (i, l_bend, rho_bend, l_drift, factor) 
    if (is_lost) return
  enddo 

! untilt coords

  call tilt_coords (ele%value(tilt$), end%vec, unset$)
  call offset_coords_m (ele, param, end, unset$, .true., .true.)

  return

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains

subroutine track_period (i_pole, l_bend, rho_bend, l_drift, factor)

  real l_bend, rho_bend, l_drift, y_ave, factor

  real*8 denom, s_center, x_center, a, b, c, descrim
  real*8 xp_out, del_s, rho, x_out, radix
        
  integer i_pole, flip

! every other pole bends in the opposite direction. 
! flip =  1  -->  First pole bends particles "upward" (+x direction)
! flip = -1  -->  First pole bends "downward" (-x direction)
                     
  flip = 1 - 2*mod(i_pole, 2) 

! 1/2 of the multipole

  if (pole_multipoles_on) then
    do n = 2, 8, 2          
      call multipole_kick (-flip*knl(n), tilt(n), n, end)
    enddo
  endif

! track 1/2 pole
! (s,x) is the center of rotation wrt the end boundary

  y_ave = end%y%pos + end%y%vel * l_bend/4
  if (abs(k_z*y_ave) > 70) then
    rho = 0
  else
    rho = flip * rho_bend * (1 + end%z%vel) / cosh(k_z*y_ave)
  endif

  if (mod(n_pole, 2) == 1 .and. i_pole == 1) rho = rho * factor

  denom = sqrt(1 + end%x%vel**2)
  s_center = -rho * end%x%vel / denom  - l_bend/2
  x_center = rho / denom + end%x%pos

  radix = rho**2 - s_center**2
  if (radix < 0) then
    type *, 'ERROR IN TRACK_WIGGLER: TRAJECTORY DOES NOT INTERSECT FACE.'
    type *, '      [THAT IS, THE PARTICLE AMPLITUDE IS TOO LARGE.]'
    is_lost = .true.
    return
  endif

  x_out = x_center - sign (sqrt(radix), rho)
  xp_out = -s_center / (x_center - x_out)

  del_s = abs ((atan(xp_out) - atan(end%x%vel)) * rho)
  end%y%pos = end%y%pos + end%y%vel * del_s

! vertical edge focusing from leaving pole

  end%y%vel = end%y%vel - flip * const1 * xp_out * sinh(k_z * end%y%pos) / k_z

! track the drift

  x_out = (x_out + xp_out * l_drift) 
  end%y%pos = end%y%pos + end%y%vel * l_drift

! vertical edge focusing entering pole

  end%y%vel = end%y%vel - flip * const1 * xp_out * sinh(k_z * end%y%pos) / k_z

! track 1/2 pole

  y_ave = end%y%pos + end%y%vel * l_bend/4
  if (abs(k_z*y_ave) > 70) then
    rho = 0
  else
    rho = -flip * rho_bend * (1 + end%z%vel) / cosh(k_z*y_ave)
  endif
  if (mod(n_pole, 2) == 1 .and. i_pole == n_pole) rho = rho * factor

  denom = sqrt(1 + xp_out**2)
  s_center = -rho * xp_out / denom  - l_bend/2
  x_center = rho / denom + x_out

  radix = rho**2 - s_center**2
  if (radix < 0) then
    type *, 'ERROR IN TRACK_WIGGLER: TRAJECTORY DOES NOT INTERSECT FACE.'
    type *, '      [THAT IS, THE PARTICLE AMPLITUDE IS TOO LARGE.]'
    is_lost = .true.
    return
  endif

  end%x%pos = x_center - sign(sqrt(radix), rho)
  end%x%vel = -s_center / (x_center - end%x%pos)

  del_s = abs ((atan(xp_out) - atan(end%x%vel)) * rho)
  end%y%pos = end%y%pos + end%y%vel * del_s

! 1/2 of the multipole

  if (pole_multipoles_on) then
    do n = 2, 8, 2
      call multipole_kick (flip*knl(n), tilt(n), n, end)
    enddo
  endif
                         
end subroutine

end subroutine
