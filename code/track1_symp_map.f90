!+
! Subroutine track1_symp_map (start, ele, param, end)
!
! Particle tracking through a single element using a partially inverted taylor
! map (In PTC/FPP this is called a genfield). 
!
! Modules Needed:
!   use bmad
!
! Input:
!   start  -- Coord_struct: Starting position
!   ele    -- Ele_struct: Element
!   param  -- Param_struct:
!     %aperture_limit_on -- If .true. then %LOST will be set if the
!                 particle is outsile the aperture.
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

subroutine track1_symp_map (start, ele, param, end)

  use accelerator

  implicit none

  type (coord_struct), intent(in)  :: start
  type (coord_struct), intent(out) :: end
  type (ele_struct),   intent(inout)  :: ele
  type (param_struct), intent(inout) :: param

  real(dp) re(lnv)
  real(rdef) r0(6)

  integer i

! Make the genfield map if needed.

  if (.not. associated(ele%gen_field)) then
    allocate (ele%gen_field)
    call taylor_to_genfield (ele%taylor, ele%gen_field, r0)
  endif

! track and add the constant term back in

  call vec_bmad_to_ptc (start%vec, re(1:6))
  re = ele%gen_field * re
  call vec_ptc_to_bmad (re(1:6), end%vec)
  end%vec = end%vec + r0

end subroutine
