!+
! Subroutine pointer_to_attribute (ring, i_ele, attrib_name, do_allocation,
!                            ptr_attrib, ix_attrib, err_flag, err_print_flag)
!
! Returns a pointer to an attribute of an element with name attrib_name.
! Also checks whether the attribute can be changed directly.
!
! Attributes that cannot be changed directly are super_slave attributes (since
! these attributes are controlled by their super_lords) and attributes that
! are controlled by an overlay_lord.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring            -- Ring_struct: Ring structure.
!     %ele_(i_ele)      -- Element whose attribute is to be checked.
!   i_ele           -- Integer: Index of element in the ring structure.
!   attrib_name     -- Character*16: Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then supress
!                       printing of an error message on error.
!
! Output:
!   ptr_attrib -- Real(rdef), pointer: Pointer to the attribute.
!                     Pointer will be deassociated if there is a problem.
!   ix_attrib  -- Ineger: Index to the attribute.
!   err_flag   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

#include "CESR_platform.inc"

Subroutine pointer_to_attribute (ring, i_ele, attrib_name, do_allocation, &
                  ptr_attrib, ix_attrib, err_flag, err_print_flag)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ring_struct), target :: ring
  type (ele_struct), pointer :: ele

  real(rdef), pointer :: ptr_attrib

  integer i_ele, ix_attrib, i, ir, ix

  character*(*) attrib_name

  logical err_flag, do_allocation
  logical, optional :: err_print_flag

! init check

  err_flag = .true.
  nullify (ptr_attrib)

  if (i_ele < 1 .or. i_ele > ring%n_ele_max) then
    if (err_print_flag) print *, 'ERROR IN POINTER_TO_ATTRIBUTE: ', &
                                       'ELEMENT INDEX OUT OF BOUNDS:', i_ele
    return
  endif

! Get attribute index

  ele => ring%ele_(i_ele)
  ix_attrib = attribute_index (ele, attrib_name)

! multipole?

  if (ix_attrib >= a0$) then   ! multipole attribute

    if (.not. associated(ele%a)) then
      if (do_allocation) then
        allocate(ele%a(0:n_pole_maxx), ele%b(0:n_pole_maxx))
        ele%a = 0; ele%b = 0
      else
        if (err_print_flag) print *, 'ERROR IN POINTER_TO_ATTRIBUTE: ', &
                            'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ', ele%name
        return
      endif
    endif

    if (ix_attrib >= b0$) then
      ptr_attrib => ele%b(ix_attrib-b0$)
    else
      ptr_attrib => ele%a(ix_attrib-a0$)
    endif

  elseif (ix_attrib < 1 .or. ix_attrib > n_attrib_maxx) then
    if (err_print_flag) then
      print *, 'ERROR IN POINTER_TO_ATTRIBUTE: INVALID ATTRIBUTE: ', attrib_name
      print *, '      FOR THIS ELEMENT: ', ele%name
    endif
    return

! otherwise must be in ele%value(:) array

  else
    ptr_attrib => ele%value(ix_attrib)
  endif

! check that attribute can be adjusted.

  if (ele%control_type == super_slave$) then
    if (err_print_flag) then
      print *, 'ERROR IN POINTER_TO_ATTRIBUTE:'
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
        if (err_print_flag) then
          print '((1x, a))', &
            'ERROR IN POINTER_TO_ATTRIBUTE. THE ATTRIBUTE: ' // attrib_name, &
            '      OF ELEMENT: ' // ele%name, &
            '      IS CONTROLLED BY THE OVERLAY_LORD: ' // ring%ele_(ir)%name, &
            '      YOU CANNOT VARY THIS ATTRIBUTE DIRECTLY.'
        endif
        return
      endif
    endif
  enddo

  err_flag = .false.

end subroutine
