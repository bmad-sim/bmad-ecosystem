!+
! Subroutine set_ele_attribute (ring, i_ele, attrib_name,
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
!   use bmad
!
! Input:
!   ring            -- Ring_struct: Ring structure containing the element 
!                                                            to be modified.
!     %ele_(i_ele)      -- Element to be modified
!   i_ele           -- Integer: Index of element in the ring structure
!   attrib_name     -- Character*16: Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   attrib_value    -- Real(rdef): Attribute value.
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
!Revision 1.6  2003/03/04 16:03:29  dcs
!VMS port
!
!Revision 1.5  2003/02/12 18:15:45  dcs
!Added pointer to pointer_to_attrubute
!
!Revision 1.4  2003/01/27 14:40:42  dcs
!bmad_version = 56
!
!Revision 1.3  2002/11/06 06:48:32  dcs
!Changed arg array
!
!Revision 1.1  2002/06/13 15:07:21  dcs
!Merged with FPP/PTC
!
!Revision 1.3  2002/02/23 20:32:24  dcs
!Double/Single Real toggle added
!
!Revision 1.2  2001/09/27 18:31:57  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"

subroutine set_ele_attribute (ring, i_ele, attrib_name, &
                                attrib_value, err_flag, make_mat6_flag, orbit_)

  use bmad

  implicit none

  type (ring_struct) :: ring
  type (coord_struct), optional :: orbit_(0:)

  real(rdef) attrib_value
  real(rdef), pointer :: ptr_attrib

  integer i_ele
  integer i, ix, ir, ix_attrib

  character*(*) attrib_name

  logical make_mat6_flag, err_flag

! error checks

  if (i_ele < 1 .or. i_ele > ring%n_ele_max) then
    print *, 'ERROR IN SET_ELE_ATTRIBUTE: ELEMENT INDEX OUT OF BOUNDS:', i_ele
    err_flag = .true.
    return
  endif

  call pointer_to_attribute (ring%ele_(i_ele), attrib_name, .true., &
                                             ptr_attrib, ix_attrib, err_flag)
  if (err_flag) return

  call check_attrib_free (ring%ele_(i_ele), ix_attrib, ring, err_flag)
  if (err_flag) return

! setting the attribute value is trivial

  ptr_attrib = attrib_value

! bookkeeping

  if (associated(ring%ele_(i_ele)%taylor(1)%term)) &
                              call kill_taylor(ring%ele_(i_ele)%taylor)
  if (associated(ring%ele_(i_ele)%gen_field)) &
                              call kill_gen_field (ring%ele_(i_ele)%gen_field)

  if (make_mat6_flag) then
    call ring_make_mat6 (ring, i_ele, orbit_)
  else
    call control_bookkeeper (ring, i_ele)
  endif

  err_flag = .false.

end subroutine
