!+
! Subroutine get_field_ele_list (ele, field_eles, dz_offset, n_field_ele)
!
! Subroutine to get the list of elements that specify the field for the given element.
!
! This is a list of elements that possibly have field info for a given element. 
! For example: a slice_slave under a super_slave under a multipass_slave. In this case
! the multipass_lord will store the field info.
! Note: groups, overlays, and girders will never have field info.
!
! If the given element does not have any lords, 
! the element list will just consist of the element itself.
!
! This routine will increase the size of field_ele_list if needed but will
! not decrease it.
!
! Input:
!   lat   -- lat_struct: Lattice
!   ele   -- Ele_struct: Element whose fields are to be evaluated.
!
! Output:
!   field_eles(:) -- Ele_pointer_struct, allocatable :: Array of field_eles.
!   dz_offset(:)  -- Offsets of ele from field elements.
!   n_field_ele   -- Integer: Number of field_eles.
!-

subroutine get_field_ele_list (ele, field_eles, dz_offset, n_field_ele)

use bmad_routine_interface, dummy => get_field_ele_list

implicit none

type (ele_struct), target :: ele
type (ele_pointer_struct), allocatable :: field_eles(:)

real(rp), allocatable :: dz_offset(:)
real(rp) offset
integer n_field_ele

!

n_field_ele = 0
if (.not. allocated(field_eles)) call re_allocate_eles (field_eles, 3)
allocate (dz_offset(3))

offset = 0
call get_field_eles (ele, offset)

!--------------------------------------------------------------------------
contains

recursive subroutine get_field_eles (this_ele, thiz_offset)

type (ele_struct), target :: this_ele
type (ele_struct), pointer :: field_ele
real(rp) thiz_offset, thiz_offset2
integer i, ix

!

if (this_ele%field_calc == refer_to_lords$) then
  do i = 1, this_ele%n_lord
    field_ele => pointer_to_lord(this_ele, i)

    select case (field_ele%key)
    case (overlay$, group$, girder$); cycle
    end select

    if (field_ele%lord_status == multipass_lord$) cycle

    thiz_offset2 = thiz_offset + this_ele%s_start - field_ele%s_start

    call get_field_eles (field_ele, thiz_offset2)

  enddo

else
  n_field_ele = n_field_ele + 1
  if (n_field_ele > size(field_eles) .or. n_field_ele > size(dz_offset)) then
    call re_allocate_eles(field_eles, n_field_ele + 10, .true.)
    call re_allocate (dz_offset, n_field_ele + 10, .false.)
  endif
  field_eles(n_field_ele)%ele => this_ele
  dz_offset(n_field_ele) = thiz_offset
endif

end subroutine get_field_eles

end subroutine get_field_ele_list
