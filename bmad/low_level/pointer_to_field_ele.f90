!+
! function pointer_to_field_ele(ele, ix_field_ele, dz_offset) result (field_ele)
!
! Routine to return a pointer to one of the "field elements" associated with a given element.
!
! A "field element" associated with a given element is an element that field info for the given element. 
! For example: a slice_slave of a super_slave of a multipass_slave will have its field info
! stored in the multipass_lord.
!
! The number of associated field elements can be determined by the routine num_field_eles.
!
! Note: groups, overlays, and girders will never have field info.
! Note: An element like a quadrupole with no lords will have one associated field element which is itself.
!
! Input:
!   ele           -- ele_struct: Element with sum number of associated field elements.
!   ix_field_ele  -- integer: Index of the field element to point to. This index runs from
!                     1 to num_field_eles(ele).
!
! Output:
!   field_ele     -- ele_struct: Pointer to the field element with index ix_field_ele.
!                     Will point to null if ix_field_ele is out of range.
!   dz_offset     -- real(rp), optional: Longitudinal offset of ele upstream edge from the field ele pointed to.
!-

function pointer_to_field_ele(ele, ix_field_ele, dz_offset) result (field_ele)

use bmad_routine_interface, dummy => pointer_to_field_ele

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: field_ele, fele
integer ix_field_ele, ix
real(rp), optional :: dz_offset
real(rp) offset

!

nullify(field_ele)
if (ix_field_ele < 1) return

ix = 0
offset = 0
call iterate_over_field_eles(ele, ix, ix_field_ele, field_ele, offset)
if (present(dz_offset)) dz_offset = offset

!---------------------------------------
contains

recursive subroutine iterate_over_field_eles(ele, ixf, ix_field_ele, field_ele, offset)

type (ele_struct), target :: ele
type (ele_struct), pointer :: this_ele, field_ele
real(rp) offset
integer ixf, ix_field_ele, i

!

select case (ele%key)
case (overlay$, group$, girder$, null_ele$); return
end select

!

if (ele%slave_status == slice_slave$) then
    call iterate_over_field_eles (ele%lord, ixf, ix_field_ele, field_ele, offset)
    offset = ele%s_start - field_ele%s_start

elseif (ele%field_calc == refer_to_lords$) then
  do i = 1, ele%n_lord
    this_ele => pointer_to_lord(ele, i)

    select case (this_ele%key)
    case (overlay$, group$, girder$); cycle
    end select

    if (this_ele%lord_status /= multipass_lord$) then
      offset = offset + ele%s_start - this_ele%s_start
    endif

    call iterate_over_field_eles (this_ele, ixf, ix_field_ele, field_ele, offset)
    if (associated(field_ele)) return
  enddo

else
  ixf = ixf + 1
  if (ixf /= ix_field_ele) return
  field_ele => ele
endif

end subroutine iterate_over_field_eles

end function pointer_to_field_ele
