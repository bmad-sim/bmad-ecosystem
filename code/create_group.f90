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
!   err_print_flag -- Logical, optional: If present and False then supress.
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
!       CONTRL(i)%IX_ATTRIB = START_EDGE$      or
!                           = END_EDGE$        or
!                           = ACCORDION_EDGE$  or
!                           = SYMMETRIC_EDGE$
!
! %IX_ATTRIB = START_EDGE$ and %IX_ATTRIB = END_EDGE$ controls the
! placement of the edges of an element keeping the lat total length invariant.
! this is done by lengthening and shortening the elements to either side
! keeping the total lat length invariant.
!
! %IX_ATTRIB = ACCORDION_EDGE$ and %IX_ATTRIB = SYMMETRIC_EDGE$ moves both
! the start and end edges simultaneously.
! ACCORDION_EDGE$ moves the edges antisymmetrically (and thus there is a length
! change of the element).
! SYMMETRIC_EDGE$ moves the edges symmetrically creating a z-offset with no
! length change.
!-

subroutine create_group (lat, ix_lord, contrl, err, err_print_flag)

use bmad_struct
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
    print *, 'ERROR IN CREATE_OVERLAY: BRANCH INDEX OUT OF BOUNDS.', ix_branch
    call err_exit
  endif

  if (ix_slave <= 0 .or. ix_slave > ubound(lat%branch(ix_branch)%ele, 1)) then
    print *, 'ERROR IN CREATE_OVERLAY: INDEX OUT OF BOUNDS.', ix_slave
    call err_exit
  endif
enddo

! init

call check_controller_controls (contrl, lat%ele(ix_lord)%name, err)
if (err) return

lord => lat%ele(ix_lord)
lord%lord_status = group_lord$
lord%key = group$
lord%ix_value = command$

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

  if (ix_attrib == start_edge$ .or. ix_attrib == end_edge$ .or. &
      ix_attrib == symmetric_edge$ .or. ix_attrib == accordion_edge$) then

    if (slave%lord_status == super_lord$) then
      ix_min = lat%control(slave%ix1_slave)%ix_slave
      ix_max = lat%control(slave%ix2_slave)%ix_slave
    elseif (ix_slave < lat%n_ele_track) then
      ix_min = ix_slave
      ix_max = ix_slave
    else
      call out_io (s_error$, r_name, &
                    'A GROUP IS NOT ALLOWED TO CONTROL', &
                    'A ' // control_name(slave%slave_status), &
                    'YOU TRIED TO CONTROL: ' // slave%name)
      return
    endif

    ! now that we have the ends we find the elements to either side whose length
    ! the group can adjust

    if (ix_attrib /= end_edge$) then
      ix1 = ix_min - 1
      do 
        if (branch%ele(ix1)%value(l$) /= 0) exit
        ix1 = ix1 - 1
        if (ix1 < 0) then
          call out_io (s_error$, r_name, &
                        'START_EDGE OF CONTROLED', &
                        'ELEMENT IS AT BEGINNING OF LAT AND CANNOT BE', &
                        'VARIED FOR GROUP: ' // lord%name)
          return
        endif
      enddo
    endif

    if (ix_attrib /= start_edge$) then
      ix2 = ix_max + 1 
      do
        if (branch%ele(ix2)%value(l$) /= 0) exit
        ix2 = ix2 + 1
        if (ix2 > branch%n_ele_track) then
          call out_io (s_error$, r_name, &
                        'END_EDGE OF CONTROLED', &
                        'ELEMENT IS AT END OF LAT AND CANNOT BE', &
                        'VARIED FOR GROUP: ' // lord%name)
          return
        endif
      enddo
    endif

    ! put in coefficients

    select case (ix_attrib)

    case (start_edge$)
      call bookit (ix1, 1)
      call bookit (ix_min, -1)

    case (end_edge$)
      call bookit (ix_max, 1)
      call bookit (ix2, -1)

    case (accordion_edge$)
      call bookit (ix1, -1)
      if (ix_min == ix_max) then
        call bookit (ix_min, 2)
      else
        call bookit (ix_min, 1)
        call bookit (ix_max, 1)
      endif
      call bookit (ix2, -1)

    case (symmetric_edge$)
      call bookit (ix1, 1)
      call bookit (ix2, -1)

    end select

  ! x_limit and y_limit

  elseif (ix_attrib == x_limit$ .or. ix_attrib == y_limit$) then

    if (n_con+2 > size(lat%control)) call reallocate_control (lat, n_con+100)
    lat%control(n_con+1) = contrl(i)
    lat%control(n_con+2) = contrl(i)
    lat%control(n_con+1)%ix_lord = ix_lord
    lat%control(n_con+2)%ix_lord = ix_lord
    if (ix_attrib == x_limit$) then
      lat%control(n_con+1)%ix_attrib = x1_limit$
      lat%control(n_con+2)%ix_attrib = x2_limit$
    else
      lat%control(n_con+1)%ix_attrib = y1_limit$
      lat%control(n_con+2)%ix_attrib = y2_limit$
    endif
    n_con = n_con + 2

    ! Update controller info for the slave element

    slave%n_lord = slave%n_lord + 2
    call add_lattice_control_structs (lat, slave)
    lat%ic(slave%ic2_lord-1) = n_con - 1
    lat%ic(slave%ic2_lord-0) = n_con - 0
    if (slave%slave_status == free$) slave%slave_status = group_slave$

  ! x_limit and y_limit

  elseif (ix_attrib == aperture$) then

    if (n_con+4 > size(lat%control)) call reallocate_control (lat, n_con+100)
    lat%control(n_con+1:n_con+4) = contrl(i)
    lat%control(n_con+1:n_con+4)%ix_lord = ix_lord
    lat%control(n_con+1)%ix_attrib = x1_limit$
    lat%control(n_con+2)%ix_attrib = x2_limit$
    lat%control(n_con+3)%ix_attrib = y1_limit$
    lat%control(n_con+4)%ix_attrib = y2_limit$
    n_con = n_con + 4

    ! Update controller info for the slave element

    slave%n_lord = slave%n_lord + 4
    call add_lattice_control_structs (lat, slave)
    lat%ic(slave%ic2_lord-3) = n_con - 3
    lat%ic(slave%ic2_lord-2) = n_con - 2
    lat%ic(slave%ic2_lord-1) = n_con - 1
    lat%ic(slave%ic2_lord-0) = n_con - 0
    if (slave%slave_status == free$) slave%slave_status = group_slave$

  ! For all else without position control the group setup is simple.

  else

    ! If the slave attribute is a multipole component, make sure it exists.
    if (ix_attrib > n_attrib_maxx .and. .not. associated (slave%a_pole)) then
      call multipole_init(slave)
    endif

    free = attribute_free (slave, attribute_name(slave, ix_attrib), lat, err_print_flag)
    if (.not. free) then
      if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, &
            'SLAVE ATTRIBUTE NOT FREE TO VARY FOR GROUP LORD: ' // lord%name)
      err = .true.
    endif
    n_con = n_con + 1
    if (n_con > size(lat%control)) call reallocate_control (lat, n_con+100)
    lat%control(n_con) = contrl(i)
    lat%control(n_con)%ix_lord = ix_lord

    ! Update controller info for the slave element

    slave%n_lord = slave%n_lord + 1
    call add_lattice_control_structs (lat, slave)
    lat%ic(slave%ic2_lord) = n_con
    if (slave%slave_status == free$) slave%slave_status = group_slave$

  endif

enddo

! final bookkeping

lord%ix2_slave = n_con
lord%n_slave = n_con - lord%ix1_slave + 1
lat%n_control_max = n_con
err = .false.

!---------------------------------------------------------------------------

contains

subroutine bookit (i_ele, scale)

integer scale, i_ele

n_con = n_con + 1
if (n_con > size(lat%control)) call reallocate_control (lat, n_con+100)
lat%control(n_con)%ix_lord = ix_lord
lat%control(n_con)%ix_slave = i_ele
lat%control(n_con)%ix_attrib = l$
lat%control(n_con)%coef = scale * contrl(i)%coef

slave => branch%ele(i_ele)
slave%n_lord = slave%n_lord + 1
call add_lattice_control_structs (lat, slave)
lat%ic(slave%ic2_lord) = n_con
if (slave%slave_status == free$) slave%slave_status = group_slave$ 

end subroutine

end subroutine
