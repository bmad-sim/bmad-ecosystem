!+
! Subroutine check_ele_attribute_set (ring, i_ele, attrib_name,
!                                    ix_attrib, err_flag, err_print_flag)
!
! Subroutine to check whether a particular attribute of an element
! can be changed directly. Attributes that cannot be changed directly are
! super_slave attributes (since these attributes are controlled by their
! super_lords) and attributes that are controlled by an overlay_lord.
! the routine also checks that i_ele is not out of bounds.
!
! Modules needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ring            -- Ring_struct: Ring structure.
!     %ele_(i_ele)      -- Element whose attribute is to be checked.
!   i_ele           -- Integer: Index of element in the ring structure.
!   attrib_name     -- Character*16: Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   err_print_flag  -- Logical: If True then print error message on error.
!
! Output:
!   ix_attrib   -- Ineger: Index to the attribute if everything is OK.
!   err_flag    -- Logical: Set False if not a valid set otherwise set true
!-

Subroutine check_ele_attribute_set (ring, i_ele, attrib_name, &
                                      ix_attrib, err_flag, err_print_flag)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  integer i_ele, ix_attrib, i, ir, ix

  character*(*) attrib_name

  logical err_print_flag, err_flag

!

  err_flag = .true.

  if (i_ele < 1 .or. i_ele > ring%n_ele_max) then
    print *, 'ERROR IN RING_SET_ELE_ATTRIBUTE: ELEMENT INDEX OUT OF BOUNDS:', &
                                                                i_ele
    return
  endif

  ele => ring%ele_(i_ele)
  ix_attrib = attribute_index (ele, attrib_name)

  if (ix_attrib < 1 .or. ix_attrib > n_attrib_maxx) then
    print *, 'ERROR IN RING_SET_ELE_ATTRIBUTE: INVALID ATTRIBUTE: ', attrib_name
    print *, '      FOR THIS ELEMENT: ', ele%name
    return
  endif

  if (ele%control_type == super_slave$) then
    print *, 'ERROR IN RING_SET_ELE_ATTRIBUTE:'
    print *, '      TRYING TO VARY AN ATTRIBUTE OF: ', ele%name
    print *, '      WHICH IS A SUPER_SLAVE WILL NOT WORK.'
    print *, '      VARY THE ATTRIBUTE OF ONE OF ITS SUPER_LORDS INSTEAD.'
    return
  endif

  do i = ele%ic1_lord, ele%ic2_lord
    ix = ring%ic_(i)
    ir = ring%control_(ix)%ix_lord
    if (ring%ele_(ir)%control_type == overlay_lord$) then
      if (ring%control_(ix)%ix_attrib == ix_attrib) then
        print '(a)', &
          'ERROR IN RING_SET_ELE_ATTRIBUTE: THE ATTRIBUTE: ', attrib_name, &
          '      OF ELEMENT: ' // ele%name, &
          '      IS CONTROLLED BY THE OVERLAY_LORD: ' // ring%ele_(ir)%name, &
          '      YOU CANNOT VARY THIS ATTRIBUTE DIRECTLY.'
        return
      endif
    endif
  enddo

  err_flag = .false.

end subroutine
