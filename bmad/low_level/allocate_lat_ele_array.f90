!+
! Subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch, do_ramper_slave_setup)
!
! Subroutine to allocate or re-allocate an element array.
! The old information is saved.
! The lower bound is always 0.
!
! Input:
!   lat                     -- Lat_struct: Lattice with element array.
!     %branch(ix_branch)%ele(:)  -- Element array to reallocate.
!   upper_bound             -- Integer, Optional: Optional desired upper bound.
!                             Default: 1.3*ubound(ele(:)) or 10 if ele is not allocated.
!   ix_branch               -- Integer, optional: Branch index. Default is 0.
!   do_ramper_slave_setup   -- logical, optional: Default False. If true, setup ramper slaves.
!                             Generally this needs to be done if reallocating with a fully formed lattice.
!
! Output:
!   lat         -- Lat_struct: Lattice with element array.
!     %branch(ix_branch)%ele(:) -- Ele_struct, pointer: Resized element array.
!-

subroutine allocate_lat_ele_array (lat, upper_bound, ix_branch, do_ramper_slave_setup)

use bmad_routine_interface, dummy => allocate_lat_ele_array

implicit none

type (lat_struct), target :: lat
integer, optional :: upper_bound
integer, optional :: ix_branch
integer ix_br, i
logical, optional :: do_ramper_slave_setup

!

ix_br = integer_option (0, ix_branch)

call allocate_element_array (lat%branch(ix_br)%ele, upper_bound)
if (ix_br == 0) lat%ele => lat%branch(0)%ele

do i = 0, ubound(lat%branch(ix_br)%ele, 1)
  lat%branch(ix_br)%ele(i)%branch => lat%branch(ix_br)
enddo

lat%branch(ix_br)%ele%ix_branch = ix_br

if (logic_option(.false., do_ramper_slave_setup)) call ramper_slave_setup(lat)

end subroutine allocate_lat_ele_array

