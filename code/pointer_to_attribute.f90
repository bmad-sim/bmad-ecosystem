!+
! Subroutine pointer_to_attribute (ele, attrib_name, do_allocation,
!                            ptr_attrib, ix_attrib, err_flag, err_print_flag)
!
! Returns a pointer to an attribute of an element with name attrib_name.
! Note: use check_attrib_free to see if the attribute may be
!   varied independently.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   attrib_name     -- Character*16: Name of attribute. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then supress
!                       printing of an error message on error.
!
! Output:
!   ptr_attrib -- Real(rp), pointer: Pointer to the attribute.
!                     Pointer will be deassociated if there is a problem.
!   ix_attrib  -- Ineger: Index to the attribute.
!   err_flag   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

#include "CESR_platform.inc"

Subroutine pointer_to_attribute (ele, attrib_name, do_allocation, &
                  ptr_attrib, ix_attrib, err_flag, err_print_flag)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct), target :: ele

  real(rp), pointer :: ptr_attrib

  integer i_ele, ix_attrib, i, ir, ix

  character*(*) attrib_name

  logical err_flag, do_allocation, do_print
  logical, optional :: err_print_flag

! init check

  err_flag = .true.
  nullify (ptr_attrib)

  do_print = .true.
  if (present(err_print_flag)) do_print = err_print_flag

! Init_ele is special in that its attributes go in special places.

  if (ele%key == init_ele$) then
    select case (attrib_name)
    case ('X_POSITION')
      ptr_attrib => ele%x_position
    case ('Y_POSITION')
      ptr_attrib => ele%y_position
    case ('Z_POSITION')
      ptr_attrib => ele%z_position
    case ('THETA_POSITION')
      ptr_attrib => ele%theta_position
    case ('PHI_POSITION')
      ptr_attrib => ele%phi_position
    case ('BETA_X')
      ptr_attrib => ele%x%beta
    case ('ALPHA_X')
      ptr_attrib => ele%x%alpha
    case ('PHI_X')
      ptr_attrib => ele%x%phi
    case ('ETA_X')
      ptr_attrib => ele%x%eta
    case ('ETAP_X')
      ptr_attrib => ele%x%etap
    case ('BETA_Y')
      ptr_attrib => ele%y%beta
    case ('ALPHA_Y')
      ptr_attrib => ele%y%alpha
    case ('PHI_Y')
      ptr_attrib => ele%y%phi
    case ('ETA_Y')
      ptr_attrib => ele%y%eta
    case ('ETAP_Y')
      ptr_attrib => ele%y%etap
    case ('C11')
      ptr_attrib => ele%c_mat(1,1)
    case ('C12')
      ptr_attrib => ele%c_mat(1,2)
    case ('C21')
      ptr_attrib => ele%c_mat(2,1)
    case ('C22')
      ptr_attrib => ele%c_mat(2,2)
    case ('S')
      ptr_attrib => ele%s
    case default
      if (do_print) then
        print *, 'ERROR IN POINTER_TO_ATTRIBUTE: BAD ATTRIBTE NAME: ', &
                                                               attrib_name
        print *, '      FOR: ', trim(ele%name)
        return
      endif
    end select
    err_flag = .false.
    return
  endif

! Get attribute index

  ix_attrib = attribute_index (ele, attrib_name)

! multipole?

  if (ix_attrib >= a0$) then   ! multipole attribute

    if (.not. associated(ele%a)) then
      if (do_allocation) then
        allocate(ele%a(0:n_pole_maxx), ele%b(0:n_pole_maxx))
        ele%a = 0; ele%b = 0
      else
        if (do_print) print *, 'ERROR IN POINTER_TO_ATTRIBUTE: ', &
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
    if (do_print) then
      print *, 'ERROR IN POINTER_TO_ATTRIBUTE: INVALID ATTRIBUTE: ', attrib_name
      print *, '      FOR THIS ELEMENT: ', ele%name
    endif
    return

! otherwise must be in ele%value(:) array

  else
    ptr_attrib => ele%value(ix_attrib)
  endif


  err_flag = .false.

end subroutine
