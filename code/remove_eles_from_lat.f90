!+
! Subroutine remove_eles_from_lat (lat, check_controls)
!
! Subroutine to compress the ele(:), control(:), and ic(:) arrays to remove
! elements no longer used. Note: to mark an element for removal use:
!     lat%branch(ib)%ele(i)%key = -1
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat            -- lat_struct: Lattice to compress.
!   check_controls -- Logical, optional: If True (default) then call check_lat_controls
!                       after the split to make sure everything is ok.
!
! Output:
!   lat -- lat_struct: Compressed lattice.
!-

subroutine remove_eles_from_lat (lat, check_controls)

use bmad_struct
use bmad_interface, except => remove_eles_from_lat

implicit none
                         
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (control_struct), pointer :: ctl

type ele_index_temp
  type (lat_ele_loc_struct), allocatable :: new(:)  ! new(old_ele_index) => new_ele_index
end type
type (ele_index_temp), allocatable :: ibr(:)

integer i, j, ib, ix, i1, i2
integer, allocatable :: ic(:), control(:)

logical, optional :: check_controls

! Allocate

allocate (ibr(0:ubound(lat%branch, 1)) )
do i = 0, ubound(lat%branch, 1)
  allocate (ibr(i)%new(lat%branch(i)%n_ele_max))
enddo

allocate (control(lat%n_control_max))
allocate (ic(lat%n_ic_max))

control = 0
ic = 0

! Mark entries in control and ic arrays for deletion.

do i = 1, lat%n_control_max
  ctl => lat%control(i)
  if (lat%branch(ctl%ix_branch)%ele(ctl%ix_slave)%key == -1) control(i) = -1
  if (lat%ele(ctl%ix_lord)%key == -1) control(i) = -1
enddo

do i = 1, lat%n_ic_max
  if (control(lat%ic(i)) == -1) ic(i) = -1
enddo

! Compress lat%ele(:) array and fill in ibr(:) array.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (ib > 0) branch%ix_from_ele = ibr(branch%ix_from_branch)%new(branch%ix_from_ele)%ix_ele

  i2 = 0
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (ele%key == -1) then
      ibr(ib)%new(i)%ix_ele = -1
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
    call init_ele(branch%ele(i), ix_ele = i)
  enddo

  branch%n_ele_max = i2

enddo

! Compress lat%control() array and correct %ix_lord and %ix_slave pointers.

i2 = 0
do i = 1, lat%n_control_max
  if (control(i) == -1) cycle
  i2 = i2 + 1
  control(i) = i2
  ctl => lat%control(i)
  lat%control(i2) = ctl
  lat%control(i2)%ix_lord  = ibr(0)%new(ctl%ix_lord)%ix_ele
  lat%control(i2)%ix_slave = ibr(ctl%ix_branch)%new(ctl%ix_slave)%ix_ele
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
    i1 = ele%ix1_slave; i2 = ele%ix2_slave
    if (i1 < 1) cycle
    if (control(i1) == i1 .and. control(i2) == i2) cycle

    ele%ix1_slave = 0  ! Start with no slaves
    ele%ix2_slave = -1
    ele%n_slave = 0
    do j = i1, i2
      if (control(j) /= -1 .and. ele%ix1_slave == 0) ele%ix1_slave = control(j)
      if (control(j) /= -1) ele%ix2_slave = control(j)
    enddo
    ele%n_slave = ele%ix2_slave - ele%ix1_slave + 1
  enddo

enddo

! Correct lord info.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do i = 1, branch%n_ele_max
    ele => branch%ele(i)

    i1 = ele%ic1_lord; i2 = ele%ic2_lord
    if (i1 < 1) cycle
    if (ic(i1) == i1 .and. ic(i2) == i2) cycle

    ele%ic1_lord = 0  ! Start with no lords
    ele%ic2_lord = -1
    ele%n_lord = 0
    do j = i1, i2
      if (ic(j) /= -1 .and. ele%ic1_lord == 0) ele%ic1_lord = ic(j)
      if (ic(j) /= -1) ele%ic2_lord = ic(j)
    enddo
    ele%n_lord = ele%ic2_lord - ele%ic1_lord + 1
    if (ele%slave_status == super_slave$ .and. ele%n_lord == 0) ele%slave_status = free$
  enddo

enddo

! deallocate and do a check

deallocate (ibr, control, ic)

if (logic_option(.true., check_controls)) call check_lat_controls (lat, .true.)

end subroutine
          
