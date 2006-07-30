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

use ptc_interface_mod, except => track1_symp_map
use tpsalie_analysis, only: assignment(=), operator(*), lnv

implicit none

type (coord_struct) :: start
type (coord_struct) :: end
type (coord_struct) start2
type (ele_struct) :: ele
type (param_struct) :: param

real(dp) re(lnv)

! Put in offsets if needed.

start2 = start

if (ele%map_with_offsets) then  ! simple case
  call track1_this_body

else
  call offset_particle (ele, param, start2, set$, &
                          set_canonical = .false., set_multipoles = .false.)
  call track1_this_body
  call offset_particle (ele, param, end, unset$, &
                          set_canonical = .false., set_multipoles = .false.)
endif


!---------------------------------------------------------------------
contains

subroutine track1_this_body

! Make the genfield map if needed.

if (.not. (associated(ele%gen_field) .and. &
                              associated(ele%taylor(1)%term))) then
  if (.not. associated(ele%taylor(1)%term)) &
                              call ele_to_taylor(ele, param, start2)
  call kill_gen_field (ele%gen_field)  ! clean up if necessary
  allocate (ele%gen_field)
  call taylor_to_genfield (ele%taylor, ele%gen_field, ele%gen0)
endif

! track and add the constant term back in

call vec_bmad_to_ptc (start2%vec, re(1:6))
re = ele%gen_field * re
call vec_ptc_to_bmad (re(1:6), end%vec)
end%vec = end%vec + ele%gen0

end subroutine

end subroutine
