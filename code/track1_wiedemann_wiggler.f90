!+
! Subroutine track1_wiedemann_wiggler (start, ele, param, end)
!
! Subroutine to track through the body of a wiggler. This routine is
! used by track1.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!
! Output:
!   end    -- Coord_struct: End position
!   param  -- Param_struct:
!     %lost  -- Logical: Set T or F depending upon whether the particle
!               reaches the exit face.
!
! For wiggler tracking use the "hard edge" model of Wiedemann, 
! "Part. Acc. Ph. II", pg. 65. Note that in Wiedemann fig 2.8, l_h is 
! mistakenly drawn as being the full pole width (it is the half pole width).
!-


#include "CESR_platform.inc"

subroutine track1_wiedemann_wiggler (start, ele, param, end)
                             
  use bmad

  implicit none

  type (coord_struct)  start, end, before
  type (ele_struct)  ele
  type (param_struct)  param

  real(rdef) length, k_z, factor, dx
  real(rdef) l_original, l_start, l_end
  real(rdef) const1, const3, tan_theta
  real(rdef) x_lim, y_lim, l_period, l_bend, l_drift, rho_bend, angle
  real(rdef) knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

  integer i, j, n, n_slice, n_pole

  logical pole_multipoles_on

!

  param%lost = .false.      ! assume success
  end = start     ! transfer start to end

! check to see if the wiggler is off. 
! if so then it just looks like a drift.

  if (ele%value(rho$) == 0) then
    end%vec(1) = end%vec(1) + ele%value(l$) * end%vec(2) * (1 - end%vec(6))
    end%vec(3) = end%vec(3) + ele%value(l$) * end%vec(4) * (1 - end%vec(6))
    return
  endif

! normal case of wiggler on.
! l_original is the original length of the wiggler if the wiggler
! has been split into pieces by split_ring.
! l_start & l_end are the starting and stopping points for the tracking
! (these have been set by split_ring).

  call offset_particle (ele, param, end, set$, set_multipoles=.false.)
             
  if (ele%value(l_original$) == 0) then ! no split has been done
    length = ele%value(l$)
    n_pole = nint(ele%value(n_pole$))
    l_start = 0
    l_end = length
  else
    length = ele%value(l_original$)
    n_pole = nint(ele%value(n_pole$) * ele%value(l_original$) / ele%value(l$))
    l_start = ele%value(l_start$)
    l_end = ele%value(l_end$)
  endif

  if (n_pole == 1) then
    print *, 'ERROR IN TRACK1_WIEDEMANN_WIGGLER: WIGGLER HAS ONLY 1 POLE FOR: ', ele%name
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

  call multipole_ele_to_kt(ele, param%particle, knl, tilt, .true.)
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
  else
    factor = 1
  endif

! wiggler track

  do i = 1, n_pole 
    if (i*l_period < l_start) cycle
    call track_period (i, l_bend, rho_bend, l_drift, factor) 
    if (param%lost) return
    if (i*l_period > l_end) exit
  enddo 

! untilt coords

  call offset_particle (ele, param, end, unset$, set_multipoles=.false.)

  return

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
contains

subroutine track_period (i_pole, l_bend, rho_bend, l_drift, factor)

  type pos_vel_struct8
    real*8 pos, vel                            ! position and velocity
  end type

  type coord_struct8
    union
      map
        type (pos_vel_struct8)  x, y, z
      endmap
      map
        real*8 vec(6)
      endmap
    endunion
  end type

  type (coord_struct8) end8

  real(rdef) l_bend, rho_bend, l_drift, y_ave, factor, l_1, l_2, l_track

  real*8 denom, s_center, x_center, a, b, c, descrim
  real*8 del_s, rho, radix, x_vel_old
        
  integer i_pole, flip

! every other pole bends in the opposite direction. 
! flip =  1  -->  First pole bends particles "upward" (+x direction)
! flip = -1  -->  First pole bends "downward" (-x direction)
                     
  flip = 1 - 2*mod(i_pole, 2) 

! l_edge1 is the distance from the wiggler start to the start of the pole
! l_1 and l_2 define the first 1/2 pole

  l_1 = (i_pole - 1) * l_period  
  l_2 = l_1 + l_bend / 2

! 1/2 of the multipole

  if (l_start <= l_1 .and. l_1 < l_end) then
    if (pole_multipoles_on) then
      do n = 2, 8, 2          
         call multipole_kick (-flip*knl(n), tilt(n), n, end)
      enddo
    endif
  endif

  end8%vec = end%vec

! track 1/2 pole
! (s,x) is the center of rotation wrt the end boundary
! l_track is the actual s-distance to track through

  l_track = l_bend / 2  
  if (l_start > l_1) l_track = l_track - (l_start - l_1)
  if (l_2 > l_end)   l_track = l_track - (l_2 - l_end)

  if (l_track > 0) then
    x_vel_old = end8%x%vel
    y_ave = end8%y%pos + end8%y%vel * l_track / 2

    if (abs(k_z*y_ave) > 70) then
      rho = 0
    else
      rho = flip * rho_bend * (1 + end8%z%vel) / cosh(k_z*y_ave)
    endif

    if (i_pole == 1) rho = rho * factor

    denom = sqrt(1 + end8%x%vel**2)
    s_center = -rho * end8%x%vel / denom  - l_track
    x_center = rho / denom + end8%x%pos

    radix = rho**2 - s_center**2
    if (radix < 0) then
      type *, 'ERROR IN TRACK1_WIEDEMANN_WIGGLER: TRAJECTORY DOES NOT INTERSECT FACE.'
      type *, '      [THAT IS, THE PARTICLE AMPLITUDE IS TOO LARGE.]'
      param%lost = .true.
      return
    endif

    end8%x%pos = x_center - sign (sqrt(radix), rho)
    end8%x%vel = -s_center / (x_center - end8%x%pos)

    del_s = abs ((atan(end8%x%vel) - atan(end%x%vel)) * rho)
    end8%y%pos = end8%y%pos + end8%y%vel * del_s
  endif

! vertical edge focusing from leaving pole

  if (l_start <= l_2 .and. l_2 < l_end) then
    end8%y%vel = end8%y%vel - &
                  flip * const1 * end8%x%vel * sinh(k_z * end8%y%pos) / k_z
  endif

! track the drift

  l_1 = l_2
  l_2 = l_1 + l_drift
  l_track = l_drift  
  if (l_start > l_1) l_track = l_track - (l_start - l_1)
  if (l_2 > l_end)   l_track = l_track - (l_2 - l_end)

  if (l_track > 0) then
    end8%x%pos = end8%x%pos + end8%x%vel * l_track
    end8%y%pos = end8%y%pos + end8%y%vel * l_track
  endif

! vertical edge focusing entering pole

  l_2 = i_pole * l_period  
  l_1 = l_2 - l_bend / 2 

  if (l_start <= l_1 .and. l_1 < l_end) then
    end8%y%vel = end8%y%vel - &
                  flip * const1 * end8%x%vel * sinh(k_z * end8%y%pos) / k_z
  endif

! track 1/2 pole

  l_track = l_bend / 2  
  if (l_start > l_1) l_track = l_track - (l_start - l_1)
  if (l_2 > l_end)   l_track = l_track - (l_2 - l_end)

  if (l_track > 0) then
    x_vel_old = end8%x%vel
    y_ave = end8%y%pos + end8%y%vel * l_track/2

    if (abs(k_z*y_ave) > 70) then
      rho = 0
    else
      rho = -flip * rho_bend * (1 + end8%z%vel) / cosh(k_z*y_ave)
    endif
    if (i_pole == n_pole) rho = rho * factor

    denom = sqrt(1 + end8%x%vel**2)
    s_center = -rho * end8%x%vel / denom  - l_track
    x_center = rho / denom + end8%x%pos

    radix = rho**2 - s_center**2
    if (radix < 0) then
      type *, 'ERROR IN TRACK1_WIEDEMANN_WIGGLER: TRAJECTORY DOES NOT INTERSECT FACE.'
      type *, '      [THAT IS, THE PARTICLE AMPLITUDE IS TOO LARGE.]'
      param%lost = .true.
      return
    endif

    end8%x%pos = x_center - sign(sqrt(radix), rho)
    end8%x%vel = -s_center / (x_center - end8%x%pos)

    del_s = abs ((atan(x_vel_old) - atan(end8%x%vel)) * rho)
    end8%y%pos = end8%y%pos + end8%y%vel * del_s
  endif

! 1/2 of the multipole

  end%vec = end8%vec

  if (l_start <= l_2 .and. l_2 <= l_end) then
    if (pole_multipoles_on) then
      do n = 2, 8, 2
        call multipole_kick (flip*knl(n), tilt(n), n, end)
      enddo
    endif
  endif
                      
end subroutine

end subroutine
