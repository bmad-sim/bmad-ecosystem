!+
! Subroutine TRACK1_627 (START, ELE, PARAM, mat627, END)
!
! Particle tracking through a single element.
! Uses the 2nd order 6x27 transport matrices whenever possible.
! This is for long term tracking.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   START        -- Coord_struct: Starting position
!   ELE          -- Ele_struct: Element
!   PARAM        -- Param_struct:
!   mat627(6,27) -- Real: 6x27 2nd order transport matrix.
!
! Output:
!   END   -- Coord_struct: End position
!
! Notes:
!
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:58  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine track1_627 (start, ele, param, mat627, end)

  use bmad_struct
  use bmad_interface

  implicit none

  type (coord_struct)  start, end
  type (ele_struct)  ele
  type (param_struct)  param

  real x_kick, y_kick
  real mat627(6,27)
  real x_lim, y_lim

  logical, parameter :: set$ = .true., unset$ = .false.    

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
! initially set end = start

  end = start     ! transfer start to end

! Use half the kick at the start
! Note: Change in %vel is NOT dependent upon energy since we are using
! cononical momentum

  if (.not. ele%is_on) then
    x_kick = 0
    y_kick = 0
  else
    x_kick = ele%value(hkick$)
    y_kick = ele%value(vkick$)
    if (x_kick /= 0) end%x%vel = end%x%vel + x_kick / 2
    if (y_kick /= 0) end%y%vel = end%y%vel + y_kick / 2
  endif

!-------------------------------------------------------------------
! 2nd order tracking

  select case (ele%key)
  case (quadrupole$, sextupole$, solenoid$, sol_quad$) 
    
    call offset_coords_m (ele, param, end, set$, .false., .true.)
    call tilt_coords (ele%value(tilt$), end%vec, set$)

    call track1_order2 (end%vec, ele, mat627, end%vec)

    call tilt_coords (ele%value(tilt$), end%vec, unset$)
    call offset_coords_m (ele, param, end, unset$, .false., .true.)

!-------------------------------------------------------------------
! error

  case default
    type *, 'ERROR IN TRACK1_627: UNKNOWN ELEMENT: ', &
                                       key_name(ele%key), ele%type
    call err_exit
  end select

!-------------------------------------------------------------------
! finally end with half the kick

  if (x_kick /= 0) end%x%vel = end%x%vel + x_kick / 2
  if (y_kick /= 0) end%y%vel = end%y%vel + y_kick / 2

! check for particles outside aperture

  if (param%aperture_limit_on) then

    x_lim = ele%value(x_limit$)
    if (x_lim <= 0) x_lim = 1e10
    if (abs(end%x%pos) > x_lim) param%lost = .true.

    y_lim = ele%value(y_limit$)
    if (y_lim <= 0) y_lim = 1e10
    if (abs(end%y%pos) > y_lim) param%lost = .true.

  endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------

subroutine track1_order2 (vec_begin, ele, mat627, vec_end)

  use bmad_struct

  implicit none

  type (ele_struct) ele

  real vec_begin(6), vec_end(6), vv(27), mat627(6,27)

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


