!+
! Subroutine element_slice_iterator (ele, param, i_slice, n_slice_tot, sliced_ele, s_start, s_end)
!
! Routine to create a sequence of element slices starting with a slice representing the beginning region of 
! the element and ending with a slice representing the end element region.
!
! This routine can be used for detailed tracking through an element.
! The default, if s_start and s_end are not specified, is to create slices of uniform length. 
! For a general element slice maker see:
!   create_element_slice
!
! It is assumed that:
!   * If s_start is present and i_slice = 1, s_start = 0
!   * If s_end is present and i_slice = n_slice_tot, s_end = element_length
!
! When tracking through an element, if s_start/s_end are not present, this routine should be 
! repeatedly called with i_slice starting at 1 and ending at n_slice_tot.
! If s_start/s_end are not present, the longitudinal position from the beginning of 
! the element of the slice is:
!     [ds * (i_slice-1), ds * i_slice]
! where
!     ds = element_length / n_slice_tot 
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
!   s_start     -- real(rp), optional: Starting edge of slice relative to beginning of element.
!   s_end       -- real(rp), optional: Ending edge of slice relative to beginning of element.
!
! Output:
!   slice_ele -- Ele_struct: Sliced portion of ele.
!-

subroutine element_slice_iterator (ele, param, i_slice, n_slice_tot, sliced_ele, s_start, s_end)

use bmad_routine_interface, except_dummy => element_slice_iterator

implicit none

type (ele_struct) ele, sliced_ele, ref_ele
type (lat_param_struct) param

real(rp), optional :: s_start, s_end
real(rp) s0, s1, ds
integer i_slice, n_slice_tot
logical at_upstream_end, at_downstream_end, err_flag

!

ds = ele%value(l$) / n_slice_tot
s0 = real_option((i_slice-1)*ds, s_start)
s1 = real_option(i_slice*ds, s_end)

if (i_slice == 1) then
  call transfer_ele (ele, sliced_ele, .true.)
endif

at_upstream_end = (i_slice == 1)
at_downstream_end = (i_slice == n_slice_tot)

if (i_slice == 1) then
  call create_element_slice (sliced_ele, ele, s1-s0, s0, &
                                param, at_upstream_end, at_downstream_end, err_flag, pointer_to_next_ele(ele, -1))
else
  call create_element_slice (sliced_ele, ele, s1-s0, s0, &
                                param, at_upstream_end, at_downstream_end, err_flag, sliced_ele)
endif

end subroutine
