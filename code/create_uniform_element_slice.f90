!+
! Subroutine create_uniform_element_slice (ele, param, i_slice, n_slice_tot, &
!                                                         sliced_ele, s_start, s_end)
!
! Routine to create an element that represents a slice of another element.
! This routine can be used for detailed tracking through an element.
! The default, if s_start and s_end are not specified, is to create slices of uniform length. 
! For a general element slice maker see:
!   create_element_slice
!
! When tracking through an element, this routine should be repeatedly called with
! i_slice starting at 1 and ending at n_slice_tot.
! The longitudinal position from the beginning of the element of 
! the entrance face of the i^th slice is at:
!     s_enter = s_start + ds * (i_slice-1)
! The exit face of the i^th slice is at: 
!     s_exit = s_start + ds * i_slice
! where
!     ds = (s_end - s_start) / n_slice_tot 
!
! It is assumed that this routine will be called in a loop so the set:
!     sliced_ele = ele
! will only be done when i_slice = 1.
!
! Input:
!   ele         -- Ele_struct: Element to slice and dice.
!   param       -- lat_param_struct: Lattice parameters
!   i_slice     -- Integer: Slice index
!   n_slice_tot -- Integer: Total number of slices.
!   s_start     -- real(rp), optional: Starting edge of 1st slice relative to
!                    beginning of element. Default is 0.
!   s_end       -- real(rp), optional: Ending edge of last slice relative to
!                    beginning of element. Default is element length.
!
! Output:
!   slice_ele -- Ele_struct: Sliced portion of ele.
!-

subroutine create_uniform_element_slice (ele, param, i_slice, n_slice_tot, &
                                                         sliced_ele, s_start, s_end)

use equal_mod, except_dummy => create_uniform_element_slice

implicit none

type (ele_struct) ele, sliced_ele
type (lat_param_struct) param

real(rp), optional :: s_start, s_end
real(rp) l_slice, s0, s1
integer i_slice, n_slice_tot
logical at_upstream_end, at_downstream_end, err_flag

!

s0 = real_option(0.0_rp, s_start)
s1 = real_option(ele%value(l$), s_end)

if (i_slice == 1 .and. s0 == 0 .and. s1 == ele%value(l$)) then
  call transfer_ele (ele, sliced_ele, .true.)
endif

at_upstream_end = (i_slice == 1 .and. s0 == 0)
at_downstream_end = (i_slice == n_slice_tot .and. s1 == ele%value(l$))

l_slice = (s1 - s0) / n_slice_tot
call create_element_slice (sliced_ele, ele, l_slice, s0 + l_slice * (i_slice - 1), &
                                      param, at_upstream_end, at_downstream_end, err_flag, sliced_ele)

end subroutine
