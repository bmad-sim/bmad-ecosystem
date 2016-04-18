!+
! Subroutine remove_eles_from_lat (lat, check_sanity)
!
! Subroutine to compress the ele(:), control(:), and ic(:) arrays to remove
! elements no longer used. Note: to mark an element for removal use:
!     lat%branch(ib)%ele(i)%key = -1
! The lat%control and lat%ic arrays will be appropriately adjusted.
!
! Individual lat%control(i) elements, along with the corresponding lat%ic(j), can 
! be removed by setting:
!     lat%control(i)%ix_attrib = int_garbage$
! In this case, the appropriate lord%n_slave value must have been adjusted for 
! the appropriate lord element.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat            -- lat_struct: Lattice to compress.
!   check_sanity   -- Logical, optional: If True (default) then call lat_sanity_check
!                       after the remove to make sure everything is ok.
!
! Output:
!   lat -- lat_struct: Compressed lattice.
!-

subroutine remove_eles_from_lat (lat, check_sanity)

use bmad_interface, except => remove_eles_from_lat

implicit none
                         
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, ele2, slave
type (branch_struct), pointer :: branch
type (control_struct), pointer :: ctl
type (lat_ele_loc_struct), pointer :: loc

type ele_index_temp
  type (lat_ele_loc_struct), allocatable :: new(:)  ! new(old_ele_index) => new_ele_index
end type
type (ele_index_temp), allocatable :: ibr(:)

integer i, j, ib, ix, i1, i2
integer :: ic(lat%n_ic_max), control(lat%n_control_max), control_to_ic(lat%n_control_max)

logical, optional :: check_sanity
logical err_flag

! Allocate

allocate (ibr(0:ubound(lat%branch, 1)) )
do i = 0, ubound(lat%branch, 1)
  allocate (ibr(i)%new(lat%branch(i)%n_ele_max))
enddo

control = 0
ic = 0

! Mark entries in control and ic arrays for deletion.

do i = 1, lat%n_ic_max
  control_to_ic(lat%ic(i)) = i
enddo

do i = 1, lat%n_control_max
  ctl => lat%control(i)
  slave => pointer_to_ele(lat, ctl%slave)
  if (slave%key == -1) control(i) = -1
  lord => pointer_to_ele(lat, ctl%lord)
  if (lord%key == -1) control(i) = -1
  if (ctl%ix_attrib == int_garbage$) control(i) = -1
  if (control(i) == -1 .and. control_to_ic(i) /= -1) ic(control_to_ic(i)) = -1
enddo

do i = 1, lat%n_ic_max
  if (control(lat%ic(i)) == -1) ic(i) = -1
enddo

! Compress branch%ele(:) array and fill in ibr(:) array.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%ix_from_branch > -1) branch%ix_from_ele = ibr(branch%ix_from_branch)%new(branch%ix_from_ele)%ix_ele
  if (branch%ix_from_ele < 0) branch%ix_from_branch = -1

  i2 = 0
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%key == -1) then
      ibr(ib)%new(i)%ix_ele    = -1
      ibr(ib)%new(i)%ix_branch = -1
    else
      i2 = i2 + 1
      ibr(ib)%new(i)%ix_ele    = i2
      ibr(ib)%new(i)%ix_branch = ib
      if (i2 /= i) branch%ele(i2) = ele
    endif
    if (i == branch%n_ele_track) then
       branch%n_ele_track = i2
    endif
  enddo

  do i = i2+1, branch%n_ele_max
    call init_ele(branch%ele(i), ix_ele = i, branch = branch)
  enddo

  branch%n_ele_max = i2

enddo

! Compress lat%control() array and correct %lord and %slave pointers.

i2 = 0
do i = 1, lat%n_control_max
  if (control(i) == -1) cycle
  i2 = i2 + 1
  control(i) = i2
  ctl => lat%control(i)
  lat%control(i2) = ctl
  lat%control(i2)%lord%ix_ele  = ibr(0)%new(ctl%lord%ix_ele)%ix_ele
  lat%control(i2)%slave%ix_ele = ibr(ctl%slave%ix_branch)%new(ctl%slave%ix_ele)%ix_ele
enddo

lat%n_control_max = i2

! Compress lat%ic() array

i2 = 0
do i = 1, lat%n_ic_max
  if (ic(i) == -1) cycle
  i2 = i2 + 1
  ic(i) = i2
  lat%ic(i2) = control(lat%ic(i))
enddo

lat%n_ic_max = i2

! Correct slave info.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    i1 = ele%ix1_slave; i2 = ele%ix1_slave+ele%n_slave+ele%n_slave_field-1
    if (i1 < 1) cycle
    if (control(i1) == i1 .and. control(i2) == i2) cycle

    ele%ix1_slave = 0  ! Start with no slaves
    ele%n_slave = 0
    do j = i1, i2
      if (control(j) /= -1 .and. ele%ix1_slave == 0) ele%ix1_slave = control(j)
      if (control(j) /= -1) ele%n_slave = control(j) - ele%ix1_slave + 1
    enddo
  enddo

enddo

! Correct lord info.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max
    ele => branch%ele(i)

    ! Update fork element

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      ele2 => pointer_to_next_ele(ele, 1, follow_fork = .true.)
      if (ele2%ix_ele /= 0) ele%value(ix_to_element$) = ibr(ele2%ix_branch)%new(ele2%ix_ele)%ix_ele
    endif

    ! Don't do anything if nothing needs to be modified

    i1 = ele%ic1_lord; i2 = ele%ic1_lord + ele%n_lord + ele%n_lord_field - 1
    if (i1 < 1) cycle
    if (ic(i1) == i1 .and. ic(i2) == i2) cycle

    ! Correct ic and n_lord info.

    ele%ic1_lord = 0  ! Start with no lords
    ele%n_lord = 0
    ele%n_lord_field = 0
    do j = i1, i2
      if (ic(j) /= -1 .and. ele%ic1_lord == 0) ele%ic1_lord = ic(j)
      if (ic(j) /= -1) ele%n_lord = ic(j) - ele%ic1_lord + 1
    enddo

    ! Correct slave_status

    ele%slave_status = not_a_child$
    do j = 1, ele%n_lord
      lord => pointer_to_lord(ele, j)
      select case(lord%lord_status)
      case (group_lord$, overlay_lord$, girder_lord$)
        ! Nothing to do
      case (multipass_lord$)
        ele%slave_status = multipass_slave$
      case (super_lord$)
        ele%slave_status = super_slave$
      end select
    enddo

  enddo

enddo

! 

if (allocated(ele_loc_com%branch)) then
  do i = lbound(ele_loc_com%branch, 1), ubound(ele_loc_com%branch, 1)
    do j = lbound(ele_loc_com%branch(i)%ele, 1), ubound(ele_loc_com%branch(i)%ele, 1)
      loc => ele_loc_com%branch(i)%ele(j)
      loc = ibr(loc%ix_branch)%new(loc%ix_ele)
    enddo
  enddo
endif

! do a check

if (logic_option(.true., check_sanity)) call lat_sanity_check (lat, err_flag)

end subroutine
          
