!+
! Subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation, a_ptr, err_flag, err_print_flag)
!
! DEPRECATED ROUTINE! DO NOT USE!
! Consider instead: pointer_to_attribute.
!
! Returns a pointer to an attribute of an element ele with attribute index ix_attrib.
!
! Use of this routine is restricted to attributes that have an index like k1$, tracking_method$, etc.
!
! Alternatively, consider the routine pointers_to_attribute.
! Note: Use attribute_free to see if the attribute may be varied independently.
!
! Input:
!   ele             -- Ele_struct: After this routine finishes A_ptr 
!                        will point to a variable within this element.
!   ix_attrib       -- Integer, Attribute index.
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   a_ptr      -- all_pointer_struct: Pointer to the attribute. 
!     %r           -- pointer to real attribute. Nullified if error or attribute is not real.               
!     %i           -- pointer to integer attribute. Nullified if error or attribute is not integer.
!     %l           -- pointer to logical attribute. Nullified if error or attribute is not logical.               
!   err_flag   -- Logical: Set True if attribtute not found. False otherwise.
!-

subroutine pointer_to_indexed_attribute (ele, ix_attrib, do_allocation, a_ptr, err_flag, err_print_flag)

use attribute_mod, dummy => pointer_to_indexed_attribute
use taylor_mod, only: taylor_term_index

implicit none

type (ele_struct), target :: ele
type (all_pointer_struct) :: a_ptr

integer :: ix_attrib
integer i, j, n, ix, iy, expn(6)

character(40) a_name
character(10) str
character(*), parameter :: r_name = 'pointer_to_indexed_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag

! Init

err_flag = .true.
nullify (a_ptr%r, a_ptr%i, a_ptr%l)
do_print = logic_option (.true., err_print_flag)

! Taylor

if (ele%key == taylor$ .and. ix_attrib > taylor_offset$) then
  write (str, '(i0)') ix_attrib - taylor_offset$
  n = index('123456', str(1:1))
  if (.not. associated(ele%taylor(1)%term)) then
    if (.not. do_allocation) return
    do i = 1, 6
      call init_taylor_series(ele%taylor(i), 0)
    enddo
  endif

  expn = 0
  do i = 2, len_trim(str)
    j = index('123456', str(i:i))
    expn(j) = expn(j) + 1
  enddo

  i = taylor_term_index(ele%taylor(n), expn, do_allocation)
  if (i /= 0) then
    a_ptr%r => ele%taylor(n)%term(i)%coef
    err_flag = .false.
  endif
  return
endif

! Overlay or Group

if (ele%key == overlay$ .or. ele%key == group$) then
  if (is_attribute(ix_attrib, control_var$)) then
    ix = ix_attrib - var_offset$ 
    if (ix > size(ele%control%var)) return
    a_ptr%r => ele%control%var(ix)%value
    err_flag = .false.
    return
  endif

  if (is_attribute(ix_attrib, old_control_var$)) then
    ix = ix_attrib - old_control_var_offset$ 
    if (ix > size(ele%control%var)) return
    a_ptr%r => ele%control%var(ix)%old_value
    err_flag = .false.
    return
  endif
endif

! Multipole or curvature attribute
! Note that elements that have surface curvature automatically have ele%photon allocated.

if (ix_attrib >= a0$ .and. ix_attrib <= b21$) then   
  a_name = attribute_name(ele, ix_attrib)

  ! Curvature
  if (a_name(1:4) == 'CURV') then
    read (a_name(12:12), *) ix
    read (a_name(15:15), *) iy
    if (ix > ubound(ele%photon%curvature%xy, 1) .or. iy > ubound(ele%photon%curvature%xy, 2)) return
    a_ptr%r => ele%photon%curvature%xy(ix,iy)

  ! Multipole
  else
    if (.not. associated(ele%a_pole)) then
      if (do_allocation) then
        call multipole_init (ele, magnetic$)
      else
        if (do_print) call out_io (s_error$, r_name, 'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
        return
      endif
    endif

    if (ix_attrib >= b0$) then
      a_ptr%r => ele%b_pole(ix_attrib-b0$)
    else
      a_ptr%r => ele%a_pole(ix_attrib-a0$)
    endif
  endif

! Electric Multipole 

if (ix_attrib >= a0_elec$ .and. ix_attrib <= b21_elec$) then   
  a_name = attribute_name(ele, ix_attrib)

  if (.not. associated(ele%a_pole_elec)) then
    if (do_allocation) then
      call multipole_init (ele, electric$)
    else
      if (do_print) call out_io (s_error$, r_name, 'MULTIPOLE NOT ALLOCATED FOR ELEMENT: ' // ele%name)
      return
    endif
  endif

  if (ix_attrib >= b0_elec$) then
    a_ptr%r => ele%b_pole_elec(ix_attrib-b0_elec$)
  else
    a_ptr%r => ele%a_pole_elec(ix_attrib-a0_elec$)
  endif
endif

! If none of the above

else
  select case (ix_attrib)
  case (is_on$);                          a_ptr%l => ele%is_on
  case (symplectify$);                    a_ptr%l => ele%symplectify
  case (mat6_calc_method$);               a_ptr%i => ele%mat6_calc_method
  case (tracking_method$);                a_ptr%i => ele%tracking_method

  case default
    ! Out of bounds
    if (ix_attrib < 1 .or. ix_attrib > num_ele_attrib$) then
      if (do_print) call out_io (s_error$, r_name, 'INVALID ATTRIBUTE INDEX: \i0\ ', 'FOR THIS ELEMENT: ' // ele%name, &
          i_array = [ix_attrib])
      return
    endif
    a_ptr%r => ele%value(ix_attrib)
  end select
endif

err_flag = .false.

end subroutine pointer_to_indexed_attribute 
