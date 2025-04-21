!+
! Subroutine allocate_element_array (ele, upper_bound)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
!
! Note: Use allocate_lat_ele_array instead for all ele(:) arrays that
!       are part of a lattice.
!   
! Input:
!   ele(:)      -- Ele_struct, pointer: Element array.
!   upper_bound -- Integer, Optional: Optional desired upper bound.
!                    Default: 1.3*ubound(ele(:)) or 10 if ele is not allocated.
!
! Output:
!   ele(:)      -- Ele_struct, pointer: Allocated element array.
!-

subroutine allocate_element_array (ele, upper_bound)

use bmad_routine_interface, dummy => allocate_element_array

implicit none

type (ele_struct), pointer :: ele(:)
type (ele_struct), pointer :: temp_ele(:)

integer, optional :: upper_bound
integer curr_ub, ub, i

! get new size

ub = 10
if (associated (ele)) ub = max (int(1.3*size(ele)), ub)
if (present(upper_bound))  ub = upper_bound

!  save ele if present

if (associated (ele)) then
  if (ub == ubound(ele, 1)) return
  curr_ub = min(ub, ubound(ele, 1))
  do i = curr_ub+1, ubound(ele, 1)
    call deallocate_ele_pointers(ele(i))
  enddo
  temp_ele => ele
  allocate(ele(0:ub))
  call transfer_eles (temp_ele(0:curr_ub), ele(0:curr_ub))
  deallocate (temp_ele)
else
  curr_ub = -1
  allocate(ele(0:ub))
endif

! 

do i = curr_ub+1, ub
  call init_ele (ele(i))
  ele(i)%ix_ele = i
end do

end subroutine allocate_element_array

