!+
! Subroutine pointer_to_attribute (ele, attrib_name, do_allocation,
!                            ptr_attrib, ix_attrib, err_flag, err_print_flag)
!
! Returns a pointer to an attribute of an element with name attrib_name.
! Note: Use attribute_free to see if the attribute may be
!   varied independently.
!
! Modules needed:
!   use bmad
!
! Input:
!   ele             -- Ele_struct: After this routine finishes Ptr_attrib 
!                        will point to a variable within this element.
!   attrib_name     -- Character(40): Name of attribute. Must be uppercase.
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
!   ix_attrib  -- Ineger: If applicable then this is the index to the 
!                     attribute in the ele%value(:) array.
!   err_flag   -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!-

#include "CESR_platform.inc"

Subroutine pointer_to_attribute (ele, attrib_name, do_allocation, &
                  ptr_attrib, ix_attrib, err_flag, err_print_flag)

  use bmad_struct
  use bmad_interface, except => pointer_to_attribute
  use multipole_mod

  implicit none

  type (ele_struct), target :: ele
  type (lr_wake_struct), allocatable :: lr(:)

  real(rp), pointer :: ptr_attrib

  integer ix_attrib, ix_d, n, ios, n_lr

  character(*) attrib_name
  character(40) a_name
  character(24) :: r_name = 'pointer_to_attribute'

  logical err_flag, do_allocation, do_print
  logical, optional :: err_print_flag

! init check

  err_flag = .true.
  nullify (ptr_attrib)

  do_print = logic_option (.true., err_print_flag)
  call str_upcase (a_name, attrib_name)
  ix_attrib = 0

! Check to see if the attribute is a long-range wake

  if (a_name(1:3) == 'LR(') then
    ix_d = index(a_name, ').')
    if (ix_d == 0) goto 9000 ! Error message and return
    read (a_name(4:ix_d-1), *, iostat = ios) n
    if (ios /= 0) goto 9000
    if (n < 0) goto 9000

    if (.not. associated (ele%wake)) then
      if (.not. do_allocation) goto 9100
      call init_wake (ele%wake, 0, 0, 0, n)
    endif

    n_lr = size(ele%wake%lr)
    if (n_lr < n) then
      if (.not. do_allocation) goto 9100
      allocate (lr(n_lr))
      lr = ele%wake%lr
      deallocate (ele%wake%lr)
      allocate (ele%wake%lr(n))
      ele%wake%lr = lr_wake_struct (0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                    0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0, .false.)
      ele%wake%lr(1:n_lr) = lr
      deallocate (lr)
    endif

    select case (a_name(ix_d+2:))
    case ('FREQ');      ptr_attrib => ele%wake%lr(n)%freq
    case ('R_OVER_Q');  ptr_attrib => ele%wake%lr(n)%r_over_q
    case ('Q');         ptr_attrib => ele%wake%lr(n)%q
!    case ('M');         ptr_attrib => ele%wake%lr(n)%m
    case ('ANGLE');     ptr_attrib => ele%wake%lr(n)%angle
    case default;       goto 9000
    end select    

    err_flag = .false.
    return

  endif

! Init_ele is special in that its attributes go in special places.

  if (ele%key == init_ele$) then
    select case (a_name)
    case ('X_POSITION')
      ptr_attrib => ele%floor%x
    case ('Y_POSITION')
      ptr_attrib => ele%floor%y
    case ('Z_POSITION')
      ptr_attrib => ele%floor%z
    case ('THETA_POSITION')
      ptr_attrib => ele%floor%theta
    case ('PHI_POSITION')
      ptr_attrib => ele%floor%phi
    case ('PSI_POSITION')
      ptr_attrib => ele%floor%psi
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
    case ('BEAM_ENERGY')
      ptr_attrib => ele%value(beam_energy$)
    case ('S')
      ptr_attrib => ele%s
    case default
      goto 9000 ! Error message and return
    end select
    err_flag = .false.
    return
  endif

! Get attribute index

  ix_attrib = attribute_index (ele, a_name)

! multipole?

  if (ix_attrib >= a0$ .and. ix_attrib <= b20$) then   ! multipole attribute

    if (.not. associated(ele%a)) then
      if (do_allocation) then
        call multipole_init (ele)
      else
        if (do_print) call out_io (s_error$, r_name, &
                        'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
        return
      endif
    endif

    if (ix_attrib >= b0$) then
      ptr_attrib => ele%b(ix_attrib-b0$)
    else
      ptr_attrib => ele%a(ix_attrib-a0$)
    endif

  elseif (ix_attrib < 1 .or. ix_attrib > n_attrib_maxx) then
    goto 9000 ! Error message and return

! otherwise must be in ele%value(:) array

  else
    ptr_attrib => ele%value(ix_attrib)
  endif


  err_flag = .false.
  return

! Error message and return

9000 continue
  if (do_print) call out_io (s_error$, r_name, &
          'INVALID ATTRIBUTE: ' // a_name, 'FOR THIS ELEMENT: ' // ele%name)
  return

9100 continue
  if (do_print) call out_io (s_error$, r_name, &
                 'WAKE ATTRIBUTE NOT ALLOCATED: ' // a_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
  return

end subroutine
