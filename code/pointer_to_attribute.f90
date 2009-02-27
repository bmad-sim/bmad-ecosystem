!+
! Subroutine pointer_to_attribute (ele, attrib_name, do_allocation,
!                            ptr_attrib, err_flag, err_print_flag, ix_attrib)
!
! Returns a pointer to an attribute of an element ele with attribute name attrib_name.
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: Alternatively consider the routine pointers_to_attribute.
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
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!   ix_attrib  -- Integer, optional: If applicable then this is the index to the 
!                     attribute in the ele%value(:), ele%a_pole(:) or ele%b_pole arrays.
!-

#include "CESR_platform.inc"

Subroutine pointer_to_attribute (ele, attrib_name, do_allocation, &
                  ptr_attrib, err_flag, err_print_flag, ix_attrib)

use bmad_struct
use bmad_interface, except_dummy => pointer_to_attribute
use multipole_mod

implicit none

type (ele_struct), target :: ele
type (lr_wake_struct), allocatable :: lr(:)

real(rp), pointer :: ptr_attrib

integer, optional :: ix_attrib
integer ix_d, n, ios, n_lr, ix_a, ix1, ix2

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
if (present(ix_attrib)) ix_attrib = 0

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
                                  0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0, .false.)
    ele%wake%lr(1:n_lr) = lr
    deallocate (lr)
  endif

  select case (a_name(ix_d+2:))
  case ('FREQ');      ptr_attrib => ele%wake%lr(n)%freq
  case ('R_OVER_Q');  ptr_attrib => ele%wake%lr(n)%r_over_q
  case ('Q');         ptr_attrib => ele%wake%lr(n)%q
  case ('ANGLE');     ptr_attrib => ele%wake%lr(n)%angle
  case default;       goto 9000
  end select    

  err_flag = .false.
  return

endif

! Special cases

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
case ('BETA_A')
  ptr_attrib => ele%a%beta
case ('ALPHA_A')
  ptr_attrib => ele%a%alpha
case ('PHI_A')
  ptr_attrib => ele%a%phi
case ('ETA_A')
  ptr_attrib => ele%a%eta
case ('ETAP_A')
  ptr_attrib => ele%a%etap
case ('BETA_B')
  ptr_attrib => ele%b%beta
case ('ALPHA_B')
  ptr_attrib => ele%b%alpha
case ('PHI_B')
  ptr_attrib => ele%b%phi
case ('ETA_B')
  ptr_attrib => ele%b%eta
case ('ETAP_B')
  ptr_attrib => ele%b%etap
case ('CMAT_11')
  ptr_attrib => ele%c_mat(1,1)
case ('CMAT_12')
  ptr_attrib => ele%c_mat(1,2)
case ('CMAT_21')
  ptr_attrib => ele%c_mat(2,1)
case ('CMAT_22')
  ptr_attrib => ele%c_mat(2,2)
case ('S')
  ptr_attrib => ele%s
case ('REF_TIME')
  ptr_attrib => ele%ref_time
end select

if (len(a_name) >= 6) then
  ix1 = index('123456', a_name(6:6))
  ix2 = index('123456', a_name(7:7))
  if (a_name(1:5) == "XMAT_" .and. ix1 /= 0 .and. ix2 /= 0) then
    ptr_attrib =>ele%mat6(ix1,ix2)
  endif
endif

if (associated(ptr_attrib)) then
  err_flag = .false.
  return
endif

! Must be an indexed attribute

ix_a = attribute_index (ele, a_name)
if (present(ix_attrib)) ix_attrib = ix_a
call pointer_to_indexed_attribute (ele, ix_a, do_allocation, &
                                      ptr_attrib, err_flag, err_print_flag)
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
