!+
! Function pointer_to_wake_ele (ele, delta_s) result (wake_ele)
!
! Routine to return a pointer to the element having the wake info associated with a given element.
! For most elements, wake_ele and ele are idential. The exceptions are super_slave and slice_slave
! elements. For super_slave and slice_slave elements, there will be an associated wake if a lord has
! a wake and the wake's location (which is the center of the ord element) is within the slave element
!
! Input:
!   ele         -- ele_struct: Lattice element.
!
! Output:
!   wake_ele    -- ele_struct: Element having the associated wake.
!                   wake_ele will be nullified if there is no associated wake.
!   delta_s     -- real(rp), optional: distance of wake locaiton from beginning of ele.
!-

function pointer_to_wake_ele (ele, delta_s) result (wake_ele)

use bmad_routine_interface, except_dummy => pointer_to_wake_ele

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: wake_ele

real(rp), optional :: delta_s
real(rp) ds
integer i

! It is assumed that a super_slave has at most one lord with wakes.
! To prevent double counting due to roundoff errors when the super_lord has been split evenly in two so 
! that the split comes right at the wake location, we shift the wake location by an insignificant length.

select case (ele%slave_status)
case (super_slave$, slice_slave$)
  do i = 1, ele%n_lord
    wake_ele => pointer_to_lord(ele, i)
    if (.not. associated(wake_ele%wake)) cycle
    ds = (1.0_rp - 1e-10_rp) * wake_ele%value(l$) / 2 + wake_ele%s_start - ele%s_start
    if (ds < 0  .or. ele%value(l$) <= ds) cycle
    if (present(delta_s)) delta_s = ds
    return
  enddo
  
  nullify(wake_ele)

!

case default
  if (associated(ele%wake)) then
    wake_ele => ele
    if (present(delta_s)) delta_s = ele%value(l$) / 2
  else
    nullify(wake_ele)
  endif
end select

end function pointer_to_wake_ele 
