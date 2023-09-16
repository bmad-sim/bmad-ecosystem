!+
! Subroutine control_bookkeeper (lat, ele, err_flag)
!
! Subroutine to transfer attibute information from lord to slave elements.
!
! If ele argument is present, bookkeeping will include all the slaves of ele
! but none of the lords.
!
! Note: This subroutine will call attribute_bookkeeper.
! Note: To do a complete bookkeeping job on a lattice use:
!   lattice_bookkeeper
!
! Input:
!   lat      -- lat_struct: lattice to be used
!   ele      -- ele_struct, optional: Element whose attribute values 
!                 have been changed. If not present bookkeeping will be done 
!                 for all elements.
!   err_flag -- logical, optional: Set True if there is an error. False otherwise.
!-

subroutine control_bookkeeper (lat, ele, err_flag)

use bookkeeper_mod, dummy => control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (ele_struct), optional :: ele
type (ele_struct), pointer :: slave, lord, branch_ele, ele2
type (branch_struct), pointer :: branch

integer i, j, ie, ib, n1, n2

logical, optional :: err_flag
logical err

character(*), parameter :: r_name = 'control_bookkeeper'

!----------------------------------------------------------------
! If ele is present we only do bookkeeping for this one element and its slaves

if (present(ele)) then
  call control_bookkeeper1 (lat, ele, .true., err)
  if (present(err_flag)) err_flag = err
  return
endif

!----------------------------------------------------------------
! Else we need to make up all the lords...
! First mark all the elements needing bookkeeping

if (present(err_flag)) err_flag = .false.

if (bmad_com%auto_bookkeeper) then
  lat%ele(:)%bookkeeping_state%control = stale$  ! Bookkeeping done on this element yet?
endif

! Bookkkeeping is done from the top level down.
! The top level elements are those lord elements that have no lords on top of them.

ie_loop: do ie = lat%n_ele_track+1, lat%n_ele_max
  ele2 => lat%ele(ie)
  if (ele2%key == null_ele$) then
    ele2%bookkeeping_state%control = ok$
    ele2%bookkeeping_state%attributes = ok$
    cycle
  endif
  if (ele2%n_lord > 0) cycle
  call control_bookkeeper1 (lat, ele2, .false., err)
  if (err .and. present(err_flag)) err_flag = .true.
enddo ie_loop

! And now bookkeeping for the elements in the tracking lattice

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (.not. bmad_com%auto_bookkeeper .and. branch%param%bookkeeping_state%control /= stale$ .and. &
                                           branch%param%bookkeeping_state%attributes /= stale$) cycle

  do ie = 0, branch%n_ele_track
    ele2 => branch%ele(ie)
    if (ele2%bookkeeping_state%control /= stale$ .and. ele2%bookkeeping_state%attributes /= stale$) cycle
    call attribute_bookkeeper (ele2)
    ele2%bookkeeping_state%control = ok$
  enddo

  branch%param%bookkeeping_state%attributes = ok$
  branch%param%bookkeeping_state%control = ok$
enddo

lat%lord_state%control = ok$
lat%lord_state%attributes = ok$

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine control_bookkeeper1 (lat, ele, force_bookkeeping, err_flag)
!
! This routine is for control bookkeeping for a single element.
! This subroutine is only to be called from control_bookkeeper and is
! not meant for general use.
!-

recursive subroutine control_bookkeeper1 (lat, ele, force_bookkeeping, err_flag)

type (lat_struct), target :: lat
type (ele_struct) ele
type (ele_struct), pointer :: slave

integer i

logical call_a_bookkeeper, force_bookkeeping
logical err_flag

! Only do bookkeeping on this element if it is stale or bookkeeping is forced by the calling routine.

err_flag = .false.

if (ele%bookkeeping_state%control == stale$ .or. ele%bookkeeping_state%attributes == stale$ .or. force_bookkeeping) then

  ! First make sure the attribute bookkeeping for this element is correct since
  ! the makeup_*_slave routines may need it.

  call attribute_bookkeeper (ele, force_bookkeeping)

  ! Slave bookkeeping

  call_a_bookkeeper = .false.

  if (ele%slave_status == super_slave$) then
    ! Attrubute bookkeeping is done in the makeup_super_slave
    call makeup_super_slave (lat, ele, err_flag)

  elseif (ele%slave_status == multipass_slave$) then
    call makeup_multipass_slave (lat, ele, err_flag)
    if (ele%n_lord > 1) call makeup_control_slave (lat, ele, err_flag)
    call_a_bookkeeper = .true.

  elseif (ele%n_lord > 0 .and. ele%slave_status /= slice_slave$) then
    call makeup_control_slave (lat, ele, err_flag)
    call_a_bookkeeper = .true.

  endif

  ! Lord bookkeeping

  if (ele%key == group$) then
    call makeup_group_lord (lat, ele, err_flag)
    call_a_bookkeeper = .true.
  endif

  ! If bookkeeping has been done by a makeup_*_slave routine then
  ! attribute_bookkeeper must be called again.
  ! This is true even if the lattice is static since a slave element
  ! can have its lord's dependent attribute values.
  ! Example: super_slave will, at this point, have its lord's num_steps value but 
  ! num_steps in the slave is different from the lord due to differences in length.

  if (call_a_bookkeeper) call attribute_bookkeeper (ele, force_bookkeeping)

  ele%bookkeeping_state%control = ok$

endif

! Recursively call this routine on the slaves

do i = 1, ele%n_slave
  if (err_flag) return
  slave => pointer_to_slave (ele, i)
  call control_bookkeeper1 (lat, slave, force_bookkeeping, err_flag)
enddo

end subroutine control_bookkeeper1

end subroutine control_bookkeeper

