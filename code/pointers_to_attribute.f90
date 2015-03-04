!+
! Subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation,
!                     ptr_array, err_flag, err_print_flag, eles, ix_attrib)
!
! Returns an array of pointers to an attribute with name attrib_name within 
! elements with name ele_name.
! Note: ele_name = 'BEAM_START' corresponds to the lat%beam_start substructure. 
! Note: ele_name can be a list of element indices. For example:
!           ele_name = "3:5"
!  This sets elements 3, 4, and 5 in the lat%ele(:) array.
! Note: ele_name can be in key:name format and include the wild card characters "*" and "%".
!       For example: "quad:q*"
! Note: Use attribute_free to see if the attribute may be varied independently.
! Note: When using wild cards, it is *not* an error if some of the matched elements do 
!       not have the the required attribute as long as at least one does.
! Note: Alternatively consider the routines:
!     ele_attribute_value
!     set_ele_attribute
!
! Modules needed:
!   use bmad
!
! Input:
!   lat             -- lat_struct: Lattice.
!   ele_name        -- Character(*): Element name. Must be uppercase
!   attrib_name     -- Character(*): Attribute name. Must be uppercase.
!                       For example: "HKICK".
!   do_allocation   -- Logical: If True then do an allocation if needed.
!                       EG: The multipole An and Bn arrays need to be allocated
!                       before their use.
!   err_print_flag  -- Logical, optional: If present and False then suppress
!                       printing of an error message on error.
!
! Output:
!   ptr_array(:) -- Real_pointer_struct, allocatable: Pointer to the attribute.
!                     Pointer will be deassociated if there is a problem.
!   err_flag     -- Logical: Set True if attribtute not found or attriubte
!                     cannot be changed directly.
!   eles(:)      -- Ele_pointer_struct, optional, allocatable: Array of element pointers.
!   ix_attrib    -- Integer, optional: If applicable then this is the index to the 
!                     attribute in the ele%value(:) array.
!-

Subroutine pointers_to_attribute (lat, ele_name, attrib_name, do_allocation, &
                        ptr_array, err_flag, err_print_flag, eles, ix_attrib)

use bmad_interface, except_dummy => pointers_to_attribute

implicit none

type (lat_struct), target :: lat
type (ele_struct), target :: beam_start
type (real_pointer_struct), allocatable :: ptr_array(:)
type (real_pointer_struct), allocatable :: ptrs(:)
type (ele_pointer_struct), optional, allocatable :: eles(:)
type (ele_pointer_struct), allocatable :: eles2(:)
type (all_pointer_struct) a_ptr

integer, optional :: ix_attrib
integer n, i, ix, key, ix_a, n_loc

character(*) ele_name
character(100) ele_name_temp
character(*) attrib_name
character(24) :: r_name = 'pointers_to_attribute'

logical err_flag, do_allocation, do_print
logical, optional :: err_print_flag

! init

err_flag = .false.
do_print = logic_option (.true., err_print_flag)

select case (ele_name)

! Selected parameters in bmad_com

case ('BMAD_COM')

  if (present(eles)) then
    call re_allocate_eles (eles, 0)
  endif

  call re_allocate (ptr_array, 1)

  select case(attrib_name)
  case ('MAX_APERTURE_LIMIT');          ptr_array(1)%r => bmad_com%max_aperture_limit
  case ('DEFAULT_DS_STEP');             ptr_array(1)%r => bmad_com%default_ds_step
  case ('SIGNIFICANT_LENGTH');          ptr_array(1)%r => bmad_com%significant_length
  case ('REL_TOL_TRACKING');            ptr_array(1)%r => bmad_com%rel_tol_tracking
  case ('ABS_TOL_TRACKING');            ptr_array(1)%r => bmad_com%abs_tol_tracking
  case ('REL_TOL_ADAPTIVE_TRACKING');   ptr_array(1)%r => bmad_com%rel_tol_adaptive_tracking
  case ('ABS_TOL_ADAPTIVE_TRACKING');   ptr_array(1)%r => bmad_com%abs_tol_adaptive_tracking
  case ('INIT_DS_ADAPTIVE_TRACKING');   ptr_array(1)%r => bmad_com%init_ds_adaptive_tracking
  case ('MIN_DS_ADAPTIVE_TRACKING');    ptr_array(1)%r => bmad_com%min_ds_adaptive_tracking
  case ('FATAL_DS_ADAPTIVE_TRACKING');  ptr_array(1)%r => bmad_com%fatal_ds_adaptive_tracking
  case default
    if (do_print) call out_io (s_error$, r_name, &
             'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
    deallocate (ptr_array)
    err_flag = .true.
  end select

  return

! like beam_start

case ('BEAM_START')

  if (present(eles)) then
    call re_allocate_eles (eles, 0)
  endif

  call re_allocate (ptr_array, 1)

  select case(attrib_name)
  case ('EMITTANCE_A'); ptr_array(1)%r => lat%a%emit 
  case ('EMITTANCE_B'); ptr_array(1)%r => lat%b%emit
  case ('EMITTANCE_Z'); ptr_array(1)%r => lat%z%emit

  case ('SPIN_X', 'SPIN_Y', 'SPIN_Z', 'SPINOR_POLARIZATION', 'SPINOR_THETA', 'SPINOR_PHI', 'SPINOR_XI')
    i = attribute_index(lat%beam_start_ele, attrib_name)
    ptr_array(1)%r => lat%beam_start_ele%value(i)

  case default
    beam_start%key = def_beam_start$
    ix = attribute_index (beam_start, attrib_name)
    if (ix < 1) then
      if (do_print) call out_io (s_error$, r_name, &
             'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
      deallocate (ptr_array)
      err_flag = .true.
      return
    endif
    if (present(ix_attrib)) ix_attrib = ix
    select case (ix)
    case (x$)
      ptr_array(1)%r => lat%beam_start%vec(1)
    case (px$)
      ptr_array(1)%r => lat%beam_start%vec(2)
    case (y$)
      ptr_array(1)%r => lat%beam_start%vec(3)
    case (py$)
      ptr_array(1)%r => lat%beam_start%vec(4)
    case (z$)
      ptr_array(1)%r => lat%beam_start%vec(5)
    case (pz$)
      ptr_array(1)%r => lat%beam_start%vec(6)
    case (field_x$)
      ptr_array(1)%r => lat%beam_start%field(1)
    case (field_y$)
      ptr_array(1)%r => lat%beam_start%field(2)
    case (phase_x$)
      ptr_array(1)%r => lat%beam_start%phase(1)
    case (phase_y$)
      ptr_array(1)%r => lat%beam_start%phase(2)
    case (t$)
      ptr_array(1)%r => lat%beam_start%t
    case (e_photon$)
      ptr_array(1)%r => lat%beam_start%p0c
    case default
      if (do_print) call out_io (s_error$, r_name, &
               'INVALID ATTRIBUTE: ' // attrib_name, 'FOR ELEMENT: ' // ele_name)
      deallocate (ptr_array)
      err_flag = .true.
    end select
  end select

  return

! This section is for things like "parameter[n_part]" whose value is stored
! in a non-standard location. Things like "parameter[e_tot]" are handled in 
! the usual way in the code section after this one.

case ('PARAMETER')

  select case(attrib_name)
  case ('N_PART')
    if (present(eles)) call re_allocate_eles (eles, 0)
    call re_allocate (ptr_array, 1)
    ptr_array(1)%r => lat%param%n_part
    return
  end select

end select

! Locate elements

call lat_ele_locator (ele_name, lat, eles2, n_loc, err_flag)
if (n_loc == 0) then
  if (do_print) call out_io (s_error$, r_name, 'ELEMENT NOT FOUND: ' // ele_name)
  if (allocated(ptr_array)) deallocate (ptr_array)
  err_flag = .true.
  return  
endif

! Locate attributes

call re_allocate (ptrs, n_loc)
n = 0
do i = 1, n_loc
  call pointer_to_attribute (eles2(i)%ele, &
          attrib_name, do_allocation, a_ptr, err_flag, .false., ix_a)
  if (err_flag .or. .not. associated(a_ptr%r)) cycle
  n = n + 1
  ptrs(n)%r => a_ptr%r
  eles2(n)%ele => eles2(i)%ele
  if (present(ix_attrib)) ix_attrib = ix_a
enddo

if (n == 0) then
  if (do_print) call out_io (s_error$, r_name, 'ATTRIBUTE: ' // attrib_name, &
                                               'NOT FOUND FOR: ' // ele_name)
  if (allocated(ptr_array)) deallocate (ptr_array)
  err_flag = .true.
  return  
endif

! Transfer pointers to ptr_array

if (present(eles)) then
  call re_allocate_eles (eles, n)
  eles = eles2(1:n)
endif

call re_allocate (ptr_array, n)
do i = 1, n
  ptr_array(i)%r => ptrs(i)%r
enddo

err_flag = .false.

end subroutine

