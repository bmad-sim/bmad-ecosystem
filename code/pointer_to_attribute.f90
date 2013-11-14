!+
! Subroutine pointer_to_attribute (ele, attrib_name, do_allocation,
!                            ptr_attrib, err_flag, err_print_flag, ix_attrib)
!
! Returns a pointer to an attribute of an element ele with attribute name attrib_name.
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: Alternatively consider the routines:
!     pointers_to_attribute
!     ele_attribute_value
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
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   ptr_attrib -- Real(rp), pointer: Pointer to the attribute.
!                     Pointer will be deassociated if there is a problem.
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!   ix_attrib  -- Integer, optional: If applicable then this is the index to the 
!                     attribute in the ele%value(:), ele%a_pole(:) or ele%b_pole arrays.
!-

subroutine pointer_to_attribute (ele, attrib_name, do_allocation, &
                  ptr_attrib, err_flag, err_print_flag, ix_attrib)

use bmad_interface, except_dummy => pointer_to_attribute

implicit none

type (ele_struct), target :: ele
type (rf_wake_lr_struct), allocatable :: lr(:)

real(rp), pointer :: ptr_attrib

integer, optional :: ix_attrib
integer ix_d, n, ios, n_lr, ix_a, ix1, ix2, n_cc, n_coef, n_v, ix, iy

character(*) attrib_name
character(40) a_name
character(24) :: r_name = 'pointer_to_attribute'

logical err_flag, do_allocation, do_print, err, out_of_bounds
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

  if (.not. associated (ele%rf_wake)) then
    if (.not. do_allocation) goto 9100
    call init_wake (ele%rf_wake, 0, 0, n)
  endif

  n_lr = size(ele%rf_wake%lr)
  if (n_lr < n) then
    if (.not. do_allocation) goto 9100
    allocate (lr(n_lr))
    lr = ele%rf_wake%lr
    deallocate (ele%rf_wake%lr)
    allocate (ele%rf_wake%lr(n))
    ele%rf_wake%lr = rf_wake_lr_struct (0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, &
                                        0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0, .false.)
    ele%rf_wake%lr(1:n_lr) = lr
    deallocate (lr)
  endif

  select case (a_name(ix_d+2:))
  case ('FREQ');      ptr_attrib => ele%rf_wake%lr(n)%freq
  case ('R_OVER_Q');  ptr_attrib => ele%rf_wake%lr(n)%r_over_q
  case ('Q');         ptr_attrib => ele%rf_wake%lr(n)%q
  case ('ANGLE');     ptr_attrib => ele%rf_wake%lr(n)%angle
  case default;       goto 9000
  end select    

  err_flag = .false.
  return

endif

! Wall3d section: s_spline

out_of_bounds = .false.

! Wall3d section: dr_ds1

! Wall3d section

if (a_name(1:10) == 'WALL.SECTION') then
  if (.not. associated(ele%wall3d)) goto 9210
  n_cc = get_cross_index(a_name, 11, err, 1, size(ele%wall3d%section))
  if (err) goto 9200

  if (a_name == 'S') then
    if (n_cc == 1) goto 9210  ! must have s = 0
    ptr_attrib => ele%wall3d%section(n_cc)%s
    err_flag = .false.
    return
  endif

  if (a_name(1:11) == 'WALL.DR_DS') then
    ptr_attrib => ele%wall3d%section(n_cc)%dr_ds
    err_flag = .false.
    return
  endif

  if (a_name(1:1) == 'V') then
    n_v = get_cross_index(a_name, 2, err, 1, size(ele%wall3d%section(n_cc)%v))
    if (err) goto 9200
    select case (a_name)
    case ('.X');        ptr_attrib => ele%wall3d%section(n_cc)%v(n_v)%x
    case ('.Y');        ptr_attrib => ele%wall3d%section(n_cc)%v(n_v)%y
    case ('.RADIUS_X'); ptr_attrib => ele%wall3d%section(n_cc)%v(n_v)%radius_x
    case ('.RADIUS_Y'); ptr_attrib => ele%wall3d%section(n_cc)%v(n_v)%radius_y
    case ('.TILT');     ptr_attrib => ele%wall3d%section(n_cc)%v(n_v)%tilt
    case default;       goto 9200
    err_flag = .false.
    end select
    return
  endif

  goto 9200
endif

! Special cases

select case (a_name)
case ('X_POSITION')
  ptr_attrib => ele%floor%r(1)
case ('Y_POSITION')
  ptr_attrib => ele%floor%r(2)
case ('Z_POSITION')
  ptr_attrib => ele%floor%r(3)
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
case ('ETA_X')
  ptr_attrib => ele%x%eta
case ('ETAP_X')
  ptr_attrib => ele%x%etap
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
case ('ETA_Y')
  ptr_attrib => ele%y%eta
case ('ETAP_Y')
  ptr_attrib => ele%y%etap
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

if (a_name(1:11) == 'CURVATURE_X' .and. a_name(13:14) == '_Y' .and. a_name(16:) == '') then
  ix = index('0123456789', a_name(12:12)) - 1
  iy = index('0123456789', a_name(15:15)) - 1
  if (ix == -1 .or. iy == -1) return
  if (ix > ubound(ele%photon%surface%curvature_xy, 1)) return
  if (iy > ubound(ele%photon%surface%curvature_xy, 2)) return
  ptr_attrib => ele%photon%surface%curvature_xy(ix,iy)
endif

if (a_name(1:5) == "XMAT_") then
  if (len(a_name) >= 7) then
    ix1 = index('123456', a_name(6:6))
    ix2 = index('123456', a_name(7:7))
    if (ix1 > 0 .and. ix2 > 0) ptr_attrib => ele%mat6(ix1,ix2)
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

!----------------------------------------
! Error message and return

9000 continue
if (do_print) call out_io (s_error$, r_name, &
          'INVALID ATTRIBUTE: ' // a_name, 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9100 continue
if (do_print) call out_io (s_error$, r_name, &
                 'WAKE ATTRIBUTE NOT ALLOCATED: ' // a_name, &
                 'FOR THIS ELEMENT: ' // ele%name)
return

!----------------------------------------
9200 continue
if (do_print) then
  if (out_of_bounds) then
    call out_io (s_error$, r_name, &
        'INDEX OUT OF BOUNDS IN ATTRIBUTE: ' // attrib_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  else
    call out_io (s_error$, r_name, &
        'MALFORMED ATTRIBUTE: ' // attrib_name, &
        'FOR THIS ELEMENT: ' // ele%name)
  endif
endif
return

!----------------------------------------
9210 continue
if (do_print) call out_io (s_error$, r_name, &
        'CROSS-SECTION NOT DEFINED SO CANNOT SET ATTRIBUTE: ' // attrib_name, &
        'FOR THIS ELEMENT: ' // ele%name)
return

!---------------------------------------------------------------
contains

!+
! Function reads number of the form "...(num)" and checks to
! see if num is between n_min and n_max.
! Function also chops "...(num)" from name.
!-

function get_cross_index(name, ix_name, err, n_min, n_max) result (ixc)

character(*) name

integer ix_name, n_min, n_max, ixc, ios

logical err

!

err = .true.

if (name(ix_name:ix_name) /= '(') return
name = name(ix_name+1:)

ix = index(name, ')')
if (ix < 2) return

read (name(1:ix-1), *, iostat = ios) ixc
if (ios /= 0 .or. name(1:ix-1) == '') return
name = name(ix+1:)

if (ixc < n_min .or. ixc > n_max) then
  out_of_bounds = .true.
  return
endif

err = .false.

end function

end subroutine
