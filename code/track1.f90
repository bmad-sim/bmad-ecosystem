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
!     %aperture_limit_on -- If True then a particle may be lost by going outside
!                of the element aperture. 
!
! Output:
!   end   -- Coord_struct: End position
!   param
!     %lost -- Set True If the particle cannot make it through an element.
!                Set False otherwise.
!
! Notes:
! It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
!-

#include "CESR_platform.inc"

subroutine track1 (start, ele, param, end)

  use bmad
  use mad_mod
  use boris_mod

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (coord_struct) orb
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

  integer tracking_method

! Init

  param%lost = .false.  ! assume everything will be OK
  call attribute_bookkeeper (ele, param)

! check for particles outside aperture

  if (ele%aperture_at == entrance_end$ .or. ele%aperture_at == both_ends$) &
                                  call check_aperture_limit (start, ele, param)
  if (param%lost) then
    param%end_lost_at = entrance_end$
    end%vec = 0       ! it never got to the end so zero this.
    return
  endif

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
    call track1_runge_kutta (orb, ele, param, end, track_com)

  case (linear$) 
    call track1_linear (orb, ele, param, end)

  case (custom$) 
    call track1_custom (orb, ele, param, end)

  case (taylor$) 
    call track1_taylor (orb, ele, param, end)

  case (symp_map$) 
    call track1_symp_map (orb, ele, param, end)

  case (symp_lie_bmad$) 
    call symp_lie_bmad (ele, param, orb, end, .false., track_com)

  case (symp_lie_ptc$) 
    call track1_symp_lie_ptc (orb, ele, param, end)

  case (wiedemann$) 
    call track1_wiedemann_wiggler (orb, ele, param, end)

  case (adaptive_boris$) 
    call track1_adaptive_boris (orb, ele, param, end, track_com)

  case (boris$) 
    call track1_boris (orb, ele, param, end, track_com)

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

  if (ele%aperture_at == exit_end$ .or. ele%aperture_at == both_ends$) &
                                call check_aperture_limit (end, ele, param)
  if (param%lost) then
    param%end_lost_at = exit_end$
    return
  endif

end subroutine
