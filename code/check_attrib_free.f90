!+
! Subroutine check_attrib_free (ele, ix_attrib, ring, err_flag, err_print_flag)
!
! Subroutine to check if an attribute is free to vary.
! Attributes that cannot be changed directly are super_slave attributes (since
! these attributes are controlled by their super_lords) and attributes that
! are controlled by an overlay_lord.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: Element
!   ix_attrib       -- Integer: Index to the attribute in ele%value() array.
!   ring            -- Ring_struct: Ring structure.
!   err_print_flag  -- Logical, optional: If present and False then supress
!                       printing of an error message on error.
!
! Output:
!   err_flag   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

#include "CESR_platform.inc"

subroutine check_attrib_free (ele, ix_attrib, ring, err_flag, err_print_flag)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct) :: ring
  type (ele_struct) :: ele

  integer ix_attrib, i, ir, ix

  logical err_flag, do_print
  logical, optional :: err_print_flag

! check that attribute can be adjusted.

  err_flag = .true.

  do_print = .true.
  if (present(err_print_flag)) do_print = err_print_flag

  if (ele%control_type == super_slave$) then
    if (do_print) then
      print *, 'ERROR IN CHECK_ATTRIB_FREE:'
      print *, '      TRYING TO VARY AN ATTRIBUTE OF: ', ele%name
      print *, '      WHICH IS A SUPER_SLAVE WILL NOT WORK.'
      print *, '      VARY THE ATTRIBUTE OF ONE OF ITS SUPER_LORDS INSTEAD.'
    endif
    return
  endif

  do i = ele%ic1_lord, ele%ic2_lord
    ix = ring%ic_(i)
    ir = ring%control_(ix)%ix_lord
    if (ring%ele_(ir)%control_type == overlay_lord$) then
      if (ring%control_(ix)%ix_attrib == ix_attrib) then
        if (do_print) then
          print '((1x, a))', &
            'ERROR IN CHECK_ATTRIB_FREE. THE ATTRIBUTE: ' // &
                                             attribute_name(ele, ix_attrib), &
            '      OF ELEMENT: ' // ele%name, &
            '      IS CONTROLLED BY OVERLAY_LORD: ' // ring%ele_(ir)%name, &
            '      YOU CANNOT VARY THIS ATTRIBUTE DIRECTLY.'
        endif
        return
      endif
    endif
  enddo

  err_flag = .false.

end subroutine
