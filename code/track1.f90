!+
! Subroutine track1 (start, ele, param, end)
!
! Particle tracking through a single element. 
! This routine is simply chooses what method to track through.
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!     %aperture_limit_on -- If .true. then param%lost will be set if the
!                 particle is outside the aperture.
!
! Output:
!   end   -- Coord_struct: End position
!   param
!     %lost -- Set .true. If the particle is outside the aperture and
!                %aperture_limit_on is set. Also: %lost is set .true. if
!                the particle does not make it through a bend irregardless
!                of the the setting of %aperture_limit_on.
!
! Notes:
!
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!-

#include "CESR_platform.inc"

subroutine track1 (start, ele, param, end)

  use bmad
  use mad_mod

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (coord_struct) orb
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

  integer tracking_method

! Radiation damping and/or fluctuations for the 1st half of the element.

  if (sr_com%damping_on .or. sr_com%fluctuations_on) then
    call track1_radiation (start, ele, param, orb, start_edge$) 
  else
    orb = start
  endif

! bmad_standard handles the case when the element is turned off.

  tracking_method = ele%tracking_method
  if (.not. ele%is_on) tracking_method = bmad_standard$

  select case (tracking_method)

  case (bmad_standard$) 
    call track1_bmad (orb, ele, param, end)

  case (runge_kutta$) 
    call track1_runge_kutta (orb, ele, param, end)

  case (linear$) 
    call track1_linear (orb, ele, param, end)

  case (custom$) 
    call track1_custom (orb, ele, param, end)

  case (taylor$) 
    call track1_taylor (orb, ele, param, end)

  case (symp_map$) 
    call track1_symp_map (orb, ele, param, end)

  case (symp_lie_bmad$) 
    call symp_lie_bmad (ele, param, orb, end, .false.)

  case (symp_lie_ptc$) 
    call track1_symp_lie_ptc (orb, ele, param, end)

  case (wiedemann$) 
    call track1_wiedemann_wiggler (orb, ele, param, end)

  case (adaptive_boris$) 
    call track1_adaptive_boris (orb, ele, param, end)

  case (boris$) 
    call track1_boris (orb, ele, param, end)

  case (mad$)
    call track1_mad (orb, ele, param, end)

  case default
    print *, 'ERROR IN TRACK1: UNKNOWN TRACKING_METHOD: ', ele%tracking_method
    call err_exit

  end select

! Radiation damping and/or fluctuations for the last half of the element

  if (sr_com%damping_on .or. sr_com%fluctuations_on) then
    call track1_radiation (end, ele, param, end, end_edge$) 
  endif

! check for particles outside aperture

  call check_aperture_limit (end, ele, param)

end subroutine
