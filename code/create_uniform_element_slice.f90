!+
! Subroutine create_uniform_element_slice (ele, param, i_slice, n_slice_tot, sliced_ele)
!
! Routine to create an element that represents a slice of another element.
! This routine can be used for detailed tracking through an element.
! Slices are of uniform length. For a general element slice maker see:
!   create_element_slice
!
! i_slice should range from 1 to n_slice_tot.
! The longitudinal position from the beginning of the element of 
! the entrance face of the i^th slice is at:
!     s_start = length_ele * (i_slice-1) / n_slice_tot 
! The exit face of the i^th slice is at: 
!     s_end = length_ele * i_slice / n_slice_tot 
!
! It is assumed that this routine will be called in a loop so the set:
!     sliced_ele = ele
! will only be done when i_slice = 1.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele         -- Ele_struct: Element to slice and dice.
!   param       -- lat_param_struct: Lattice parameters
!   i_slice     -- Integer: Slice index
!   n_slice_tot -- Integer: Total number of slices.
!
! Output:
!   slice_ele -- Ele_struct: Sliced portion of ele.
!-

subroutine create_uniform_element_slice (ele, param, i_slice, n_slice_tot, sliced_ele)

use bookkeeper_mod, except_dummy => create_uniform_element_slice

implicit none

type (ele_struct) ele, sliced_ele
type (lat_param_struct) param

integer i_slice, n_slice_tot
real(rp) l_slice
logical at_entrance_end, at_exit_end, err_exit

!

if (i_slice == 1) then
  call transfer_ele (ele, sliced_ele, .true.)
endif

at_entrance_end = (i_slice == 1)
at_exit_end = (i_slice == n_slice_tot)

l_slice = ele%value(l$) / n_slice_tot
call create_element_slice (sliced_ele, ele, l_slice, l_slice * (i_slice - 1), &
                                      param, at_entrance_end, at_exit_end, err_exit)

end subroutine
