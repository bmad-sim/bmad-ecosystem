!+
! Function pointer_to_fibre (ele) result (assoc_fibre)
!
! Routine to return the reference PTC fibre associated with a given Bmad element.
!
! The reference fibre is the fibre whose upstream edge corresponds to the downstream end 
! of the bmad element. This is generally ele%ptc_fibre%next.
!
! The reference fibre is so chosen since the reference edge of a Bmad element, where such things 
! as Twiss parameters are computed, is the downstream edge while a PTC fibre uses the upstream 
! edge for the reference.
!
! Input:
!   ele   -- ele_struct: Bmad element
!
! Output:
!   assoc_fibre -- fibre, pointer: Pointer to the associated fibre.
!-

function pointer_to_fibre (ele) result (assoc_fibre)

use bmad_routine_interface, dummy => pointer_to_fibre

implicit none

type (ele_struct), target :: ele
type (fibre), pointer :: assoc_fibre
type (branch_struct), pointer :: branch

!

branch => pointer_to_branch(ele)

if (ele%ix_ele == 0) then
  assoc_fibre => branch%ele(1)%ptc_fibre
else
  assoc_fibre => ele%ptc_fibre%next
endif

end function pointer_to_fibre

