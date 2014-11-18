!+
! Subroutine create_group (lat, ix_lord, contrl, err, err_print_flag)
!
! Subroutine to create a group control element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat          -- lat_struct: Lattice
!     %ele(ix_lord)%value(command$) -- Real(rp): Initial command value.
!   ix_lord       -- Integer: Index of group lord element (see below).
!   contrl(:)     -- Control_struct: What to control.
!     %ix_slave       -- Integer: Index in lat%branch()%ele() of element controlled.
!     %ix_branch      -- Integer: Branch index.  
!     %ix_attrib      -- Integer: Index in %value() array of
!                                 attribute controlled.
!     %coef           -- Real(rp): Coefficient.
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!   err_print_flag -- Logical, optional: If present and False then suppress.
!                       printing of an error message if attribute is not free.  
!
! Output:
!     lat -- lat_struct: Appropriate values are set in the LAT structure.
!
! Note: Use NEW_CONTROL to get an index for the group element
!
! Example:
!   call new_control (lat, ix_lord)   ! get IX_LORD index
!   lord => lat%ele(ix_lord)   
!   lord%name = 'GROUP1'      ! group name
!   lord%value(command$) = 0  ! start at zero
!   n_control = 2             ! control 2 elements
!
!   contrl(1)%ix_branch = 0   ! Main lattice branch
!   contrl(1)%ix_slave = 10   ! lat%ele(10) is Q01W say.
!   contrl(1)%ix_attrib = k1$ ! The group controls the quadrupole strength.
!   contrl(1)%coef = 0.1      ! A change in the group value of 1 produces
!                             !    a change of 0.1 in k1 of element 10.
!
!   contrl(1)%ix_branch = 1   ! branch #1
!   contrl(2)%ix_slave = 790  ! lat%branch(1)%ele(790) is Q01E say.
!   contrl(2)%ix_attrib = k1$ ! The group controls the quadrupole strength.
!   contrl(2)%coef = -0.1     ! make changes antisymmetric.
!
!   call create_group (lat, ix_lord, 2, contrl)  ! create the group
!
! Notes:
!   A) The value of the group is stored in lat%ele(IX_LORD)%VALUE(COMMAND$)
!   B) Only changes from the previous value are significant. The
!      old value is stored in lat%ele(ix_lord)%value(old_command$).
!   C) Use CONTROL_BOOKKEEPER to update the attributes of the elements
!      controlled by the group element.
!   D) Use lat_make_mat6 to update the attributes AND remake MAT6 for the
!     elements controlled by the group element.
!
! Proceeding with the previous example:
!   lord%value(command$) = 10            ! put in a value for the group
!   call lat_make_mat6 (lat, ix_lord)    ! update the k1's for the 2 quads
!                                        !   AND remake the MAT6's
!   call control_bookkeeper (lat, lord)  ! use this instead
!                                        !   to only update the k1's
!   lord%value(command$) = 13            ! put in a new value for the group
!   call lat_make_mat6 (lat, ix_lord)    ! the change now in k1's is
!                                        !   based upon the delta: 13 - 10
!
! Note: You can control an element's position by setting:
!       contrl(i)%ix_attrib = start_edge$      or
!                           = end_edge$        or
!                           = accordion_edge$  or
!                           = s_position$      or
!                           = lord_pad1$       or
!                           = lord_pad2$
!
! %ix_attrib = start_edge$ and %ix_attrib = end_edge$ controls the
! placement of the edges of an element keeping the lat total length invariant.
! this is done by lengthening and shortening the elements to either side
! keeping the total lat length invariant.
!
! %ix_attrib = accordion_edge$ and %ix_attrib = s_position$ moves both
! the start and end edges simultaneously.
! accordion_edge$ moves the edges antisymmetrically (and thus there is a length
! change of the element).
! s_position$ moves the edges symmetrically creating a z-offset with no
! length change.
!-

subroutine create_group (lat, ix_lord, contrl, err, err_print_flag)

use bmad_interface, except_dummy => create_group

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: slave, lord
type (control_struct)  contrl(:)
type (branch_struct), pointer :: branch

integer i, ix_lord, ix_attrib, n_control, n_con
integer ix1, ix2, ix_min, ix_max, ix_slave, ix_branch

logical err, free
logical, optional :: err_print_flag

character(16) :: r_name = 'create_group'

! Error check

n_control = size(contrl)

do i = 1, n_control
  ix_slave  = contrl(i)%ix_slave
  ix_branch = contrl(i)%ix_branch

  if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) then
    call out_io (s_fatal$, r_name, 'BRANCH INDEX OUT OF BOUNDS. \i0\ ', ix_branch)
    if (global_com%exit_on_error) call err_exit
  endif

  if (ix_slave <= 0 .or. ix_slave > ubound(lat%branch(ix_branch)%ele, 1)) then
    call out_io (s_fatal$, r_name, 'INDEX OUT OF BOUNDS. \i0\ ', ix_slave)
    if (global_com%exit_on_error) call err_exit
  endif
enddo

! init

call check_controller_controls (contrl, lat%ele(ix_lord)%name, err)
if (err) return

lord => lat%ele(ix_lord)
lord%lord_status = group_lord$
lord%key = group$
lord%ix_value = command$
call set_ele_defaults (lord)

if (n_control == 0) return ! If no slaves then nothing to do.

err = .true.
n_con = lat%n_control_max
lord%ix1_slave = n_con + 1

! loop over all controlled elements

do i = 1, n_control

  ! For position control: We need to figure out the elements that
  ! need to be controlled.
  ! Find beginning and ending positions of element
  ! if a super_lord then we must go to the slave elements to find the ends
  ! else not a super lord so finding the ends is simple

  ix_slave = contrl(i)%ix_slave
  ix_attrib = contrl(i)%ix_attrib
  ix_branch = contrl(i)%ix_branch
  branch => lat%branch(ix_branch)
  slave => branch%ele(ix_slave)

  ! If the slave attribute is a multipole component, make sure it exists.
  if (ix_attrib > num_ele_attrib$ .and. .not. associated (slave%a_pole)) then
    call multipole_init(slave)
  endif

  ! Varying the length of a super_slave is permitted so do not check in this case.

  select case (ix_attrib)
  case (start_edge$, end_edge$, accordion_edge$, s_position$)
    free = attribute_free (slave, 'L', err_print_flag)
  case default
    free = attribute_free (slave, attribute_name(slave, ix_attrib), err_print_flag)
  end select

  if (.not. free) then
    if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, &
          'SLAVE ATTRIBUTE NOT FREE TO VARY FOR GROUP LORD: ' // lord%name)
    err = .true.
  endif

  !

  n_con = n_con + 1
  if (n_con > size(lat%control)) call reallocate_control (lat, n_con+100)
  lat%control(n_con) = contrl(i)
  lat%control(n_con)%ix_lord = ix_lord

  ! Update controller info for the slave element

  slave%n_lord = slave%n_lord + 1
  call add_lattice_control_structs (lat, slave)
  lat%ic(slave%ic2_lord) = n_con
  if (slave%slave_status == free$) slave%slave_status = control_slave$

enddo

! End stuff

lord%ix2_slave = n_con
lord%n_slave = n_con - lord%ix1_slave + 1
lat%n_control_max = n_con

err = .false.

end subroutine
