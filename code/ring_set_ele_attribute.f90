!+
! Subroutine ring_set_ele_attribute (ring, i_ele, attrib_name,
!                           attrib_value, err_flag, make_mat6_flag, orbit_)
!
! Subroutine to set the attribute of an element, propagate the change to
! any slave elements, and optionally to remake the 6x6 transfer matrix. 
! This routine is essentually equivalent to changing the attribute directly:
!     ring%ele_(i_ele)%value(attribute$) = attrib_value
! and then calling ring_make_mat6. The difference is that 
! this routine also does error checking to make sure everything is OK. 
! For example, the routine checks that you are not trying to change 
! an attribute of a super_slave.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring            -- Ring_struct: Ring structure containing the element 
!                                                            to be modified.
!     %ele_(i_ele)      -- Element to be modified
!   i_ele           -- Integer: Index of element in the ring structure
!   attrib_name     -- Character*16: Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   attrib_value    -- Real: Attribute value.
!   make_mat6_flag  -- Logical: If True then make the 6x6 transfer matrix.
!   orbit_(0:n_ele_maxx) -- Coord_struct: [Optional] closed orbit about
!                            which the 6x6 matrices are made.
!
! Output:
!   ring      -- Ring_struct: Modified ring.
!     %ele_(i_ele)  -- Modified element.
!   err_flag  -- Logical: Set True if there is an error. Otherwise set False.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine ring_set_ele_attribute (ring, i_ele, attrib_name, &
                            attrib_value, err_flag, make_mat6_flag, orbit_)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) :: ring
  type (coord_struct), optional :: orbit_(0:n_ele_maxx)

  real attrib_value

  integer i_ele
  integer i, ix, ir, ix_attrib

  character*(*) attrib_name

  logical make_mat6_flag, err_flag

! error checks

  call check_ele_attribute_set (ring, i_ele, attrib_name, ix_attrib, &
                                                              err_flag, .true.)
  if (err_flag) return

! setting the attribute value is trivial

  ring%ele_(i_ele)%value(ix_attrib) = attrib_value

! bookkeeping 

  if (make_mat6_flag) then
    call ring_make_mat6 (ring, i_ele, orbit_)
  else
    call control_bookkeeper (ring, i_ele)
  endif

  err_flag = .false.

end subroutine
