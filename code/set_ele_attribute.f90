!+
! Subroutine set_ele_attribute (ring, i_ele, attrib_name,
!               attrib_value, err_flag, make_mat6_flag, orbit_)
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
!   attrib_name     -- Character(16): Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   attrib_value    -- Real(rp): Attribute value.
!   make_mat6_flag  -- Logical: If True then make the 6x6 transfer matrix.
!   orbit_(0:)      -- Coord_struct, optional: closed orbit about
!                            which the 6x6 matrices are made.
!
! Output:
!   ring            -- Ring_struct: Modified ring.
!     %ele_(i_ele)    -- Modified element.
!   err_flag        -- Logical: Set True if there is an error. 
!                        Otherwise set False.
!-

#include "CESR_platform.inc"

subroutine set_ele_attribute (ring, i_ele, attrib_name, &
                                attrib_value, err_flag, make_mat6_flag, orbit_)

  use bmad_struct
  use bmad_interface, except => set_ele_attribute
  use bookkeeper_mod, only: control_bookkeeper

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele
  type (coord_struct), optional :: orbit_(0:)

  real(rp) attrib_value
  real(rp), pointer :: ptr_attrib

  integer i_ele, ix_attrib, key

  character(*) attrib_name

  logical make_mat6_flag, err_flag, do_bookkeeping

! error checks

  if (i_ele < 1 .or. i_ele > ring%n_ele_max) then
    print *, 'ERROR IN SET_ELE_ATTRIBUTE: ELEMENT INDEX OUT OF BOUNDS:', i_ele
    err_flag = .true.
    return
  endif

  ele => ring%ele_(i_ele)

  call pointer_to_attribute (ele, attrib_name, .true., &
                                             ptr_attrib, ix_attrib, err_flag)
  if (err_flag) return
  if (.not. attribute_free (ele, ix_attrib, ring)) return

! Setting the attribute value is trivial

  ptr_attrib = attrib_value

! bookkeeping

  if (ele%key == lcavity$ .or. ele%key == custom$) &
                                      call compute_element_energy (ring)
  
  call control_bookkeeper (ring, i_ele)

  if (make_mat6_flag) then
    do_bookkeeping = bmad_com%auto_bookkeeper
    bmad_com%auto_bookkeeper = .false.
    call ring_make_mat6 (ring, i_ele, orbit_)
    bmad_com%auto_bookkeeper = do_bookkeeping
  endif
  
  err_flag = .false.

end subroutine
