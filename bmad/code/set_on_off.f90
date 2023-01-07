!+
! Subroutine set_on_off (key, lat, switch, orb, use_ref_orb, ix_branch, saved_values, attribute)
!
! Routine to turn on or off a set of elements (quadrupoles, rfcavities, etc.) in a lattice. 
! An element that is turned off acts like a drift.
! lat_make_mat6 will be called to remake lat%ele()%mat6.
!
! This routine can also be used to turn on/off a specific attribute in a set of elements using the
! attribute and saved_values arguments. For example, to turn on/off the K2 component of all bends. 
! For real and integer attributes, switch = off$ means to set to zero.
! For real and integer attributes, switch = on$ and switch = restore_state$ have the same meaning 
! and the value that is set is optained from the saved_values argument (which must be present.)
!
! Input:
!   key             -- integer: Class name of elements to be turned on or off. [quadrupole$, etc.]
!   lat             -- lat_struct: lattice structure holding the elements.
!   switch          -- integer: 
!                       on$            => Turn elements on. If saved_values argument is present, use this.
!                                          If not present (only for logical attributes), set to True.
!                       off$           => Turn elements off (but will not store the present state).
!                       off_and_save$  => Save on/off state and then turn elements off.
!                       save_state$    => Save present on/off state. No turning on or off is done.
!                       restore_state$ => Restore saved on/off state from saved_values argument.
!   orb(0:)         -- coord_struct, optional: Needed for lat_make_mat6
!   use_ref_orb     -- logical, optional: If present and true then use 
!                        ele%map_ref_orb for the reference orbit for
!                        calculating %mat6. Default is false.
!   ix_branch       -- integer, optional: If present then only set for 
!                        this lattice branch.
!   saved_values(:) -- real(rp), allocatable, optional: Element-by element saved values of the component.
!                       Must be present if needed (EG if switch = restore_state$, etc.).
!   attribute       -- character(*), optional: Attribute to turn on/off. Eg: 'K2', 'MULTIPOLE_ON', etc.
!                       Default is 'IS_ON'. Must be upper case.
!                   
!
! Output:
!   lat             -- lat_struct: Modified lattice.
!   saved_values(:) -- real(rp), allocatable, optional: Saved values of the component.
!                       Must be present if needed (EG if switch = off_and_save$, etc.).
!-

subroutine set_on_off (key, lat, switch, orb, use_ref_orb, ix_branch, saved_values, attribute)

use changed_attribute_bookkeeper, dummy => set_on_off

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct), optional :: orb(0:)
type (all_pointer_struct) a_ptr

real(rp), optional, allocatable :: saved_values(:)
real(rp) old_state

integer :: key, switch
integer, optional :: ix_branch
integer i, ib, n_set, n_save_old

logical, optional :: use_ref_orb
logical err_flag

character(*), optional :: attribute
character(20) :: r_name = 'set_on_off'

!

n_save_old = 0
if (present(saved_values)) n_save_old = size(saved_values)

select case (switch)
case (restore_state$, off_and_save$, save_state$)
  if (present(saved_values)) then
    if (.not. allocated(saved_values)) allocate(saved_values(20))
  else
    call out_io (s_fatal$, r_name, 'SAVED_VALUES ARGUMENT NOT PRESENT.')
    if (global_com%exit_on_error) call err_exit
    return
  endif
end select

! Do not set super_slave elements since the appropriate setting can only be done in the lords.
! For example, setting a multipole parameter in a super_slave is a no-no.

do ib = 0, ubound(lat%branch, 1)

  if (present(ix_branch)) then
    if (ix_branch /= ib) cycle
  endif

  branch => lat%branch(ib)

  n_set = 0
  do i = 1, branch%n_ele_max

    ele => branch%ele(i)
    if (ele%key /= key) cycle
    if (ele%slave_status == super_slave$) cycle

    n_set = n_set + 1

    if (present(attribute)) then
      call pointer_to_attribute (ele, attribute, .true., a_ptr, err_flag)
      if (err_flag) return
    else
      a_ptr%l => ele%is_on
    endif

    if (present(saved_values)) then
      if (n_set > size(saved_values)) call re_allocate(saved_values, 2*n_set, .false.)
    endif

    old_state = value_of_all_ptr(a_ptr)

    select case (switch)
    case (on$)
      if (present(saved_values)) then
        call set_value_of(a_ptr, saved_values(n_set))
      else
        call set_value_of(a_ptr, true$)
      endif

    case (off$)
      call set_value_of(a_ptr, false$)  ! false$ = 0

    case (save_state$)
      saved_values(n_set) = value_of_all_ptr(a_ptr)
      cycle

    case (restore_state$)
      call set_value_of(a_ptr, saved_values(n_set))

    case (off_and_save$)
      saved_values(n_set) = value_of_all_ptr(a_ptr)
      call set_value_of(a_ptr, false$)

    case default
      call out_io (s_abort$, r_name, 'BAD SWITCH: \i\ ', switch)
      if (global_com%exit_on_error) call err_exit
    end select

    if (old_state == value_of_all_ptr(a_ptr)) cycle

    call set_flags_for_changed_attribute (ele, a_ptr)

    if (logic_option (.false., use_ref_orb)) then
      call make_mat6(ele, branch%param, ele%map_ref_orb_in)
    else
      call lat_make_mat6(lat, i, orb, ib)
    endif

  enddo

  if (key == rfcavity$) then   ! Reset 1-turn maps
    branch%param%t1_with_rf = 0  
    branch%param%t1_no_rf = 0    
  endif

enddo

! Error check

if (n_set == 0) return ! No setting done so no error checking needed

if ((switch == on$ .and. (associated(a_ptr%r) .or. associated(a_ptr%i))) .or. &
    (switch == on$ .and. associated(a_ptr%l) .and. present(saved_values)) .or. &
     switch == restore_state$) then
  if (n_save_old < n_set) then
    call out_io (s_fatal$, r_name, 'SAVED_VALUES IS BEING USED BEFORE BEING SET!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

call lattice_bookkeeper(lat)

!-------------------------------------------------------------------------
contains

subroutine set_value_of (a_ptr, val)

type (all_pointer_struct) a_ptr
real(rp) val

!

if (associated(a_ptr%r)) then
  a_ptr%r = val
elseif (associated(a_ptr%i)) then
  a_ptr%i = nint(val)
else
  a_ptr%l = is_true(val)
endif

end subroutine set_value_of

end subroutine set_on_off

