!+
! Subroutine track1_symp_map (start, ele, param, end)
!
! Particle tracking through a single element using a partially inverted taylor
! map (In PTC/FPP this is called a genfield). 
!
! Note: It is assumed that HKICK and VKICK are the kicks in the horizontal
! and vertical kicks irregardless of the value for TILT.
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
!   end   -- Coord_struct: End position
!-

#include "CESR_platform.inc"

subroutine track1_symp_map (start, ele, param, end)

  use accelerator

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

  real(dp) re(lnv)

  integer i

! Make the genfield map if needed.

  if (.not. (associated(ele%gen_field) .and. &
                              associated(ele%taylor(1)%term))) then
    if (.not. associated(ele%taylor(1)%term)) &
                              call ele_to_taylor(ele, start, param)
    call kill_gen_field (ele%gen_field)  ! clean up if necessary
    allocate (ele%gen_field)
    call taylor_to_genfield (ele%taylor, ele%gen_field, ele%gen0)
  endif

! track and add the constant term back in

  call vec_bmad_to_ptc (start%vec, re(1:6))
  re = ele%gen_field * re
  call vec_ptc_to_bmad (re(1:6), end%vec)
  end%vec = end%vec + ele%gen0

end subroutine
