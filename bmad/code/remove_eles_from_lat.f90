!+
! Subroutine remove_eles_from_lat (lat, check_sanity)
!
! Subroutine to compress the ele(:), branch(:), control(:), and ic(:) arrays to remove elements no longer used. 
!
! To mark element ie of branch ib for removal set:
!     lat%branch(ib)%ele(ie)%ix_ele = -1
! To mark a super_lord/multipass_lord for removal but *not* remove the slaves set:
!     lat%branch(ib)%ele(ie)%ix_ele = -2
! To remove an entire branch mark the branch index:
!     lat%branch(ib)%ix_branch = -1
! Branches where all the elements with index above 0 are marked for removal will also be removed. 
!
! Individual lat%control(i) elements, along with the corresponding lat%ic(j), can 
! be removed by setting:
!     lat%control(i)%attribute = 'REMOVE'
!
! Important: To save computation time, lattice_bookkeeper is not called by this routine.
! If you use this routine you should call lattice_bookkeeper after all lattice adjustments have been made.
!
! Notes:
!   * A fork and photon_fork element that is retained but whose "to" element is removed will be converted
!       to a marker element.
!   * If branch 0 is to be removed: The lord elements of the branch will be saved unless individually 
!       marked for deletion.
!   * When removing branches or elements, the lat%control and lat%ic arrays will be appropriately adjusted.
!   * Branch and ele arrays will be compressed. So, for example, if branch #0 is removed, then
!       branch #1 will become branch #0, etc.
!   * Whether a super_slave element is removed or not is dependent upon how its super_lord is marked for removal.
!       This is independent of how the super_slave element itself is marked.
!   * If multiple super_lord elements overlap, and if one of them is not marked to be removed, none of the elements
!       will be removed. Exception: With %ix_ele = -2, that super_lord will always be removed.
!
! Input:
!   lat            -- lat_struct: Lattice to compress.
!   check_sanity   -- logical, optional: If True (default) then call lat_sanity_check
!
! Output:
!   lat -- lat_struct: Compressed lattice.
!   ele_all_loc   -- ele_all_location_struct, optional: Used to keep track of which elements move where.
!-

subroutine remove_eles_from_lat (lat, check_sanity)

use bookkeeper_mod, except => remove_eles_from_lat

implicit none
                         
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, ele2, slave
type (branch_struct), pointer :: branch, b0, branch2, new_lord_branch
type (branch_struct), allocatable :: temp_branch(:)
type (control_struct), pointer :: ctl
type (control_struct) :: ctl0
type (lat_ele_loc_struct), pointer :: loc

type this_ele_loc_struct
  integer :: new_ix_ele = -1
  integer :: new_ix_branch = 0
  integer :: old_ix_ele = -1
  integer :: old_ix_branch = 0
end type

type ele_index_temp
  type (this_ele_loc_struct), allocatable :: loc(:)  ! maps old -> new

end type
type (ele_index_temp), allocatable :: ibr(:)

integer i, j, n, ib, ib0, ib2, ie0, ie, ie2, ix, i1, i2, icon, iz
integer :: ic(lat%n_ic_max), control(lat%n_control_max), control_to_ic(lat%n_control_max)

logical, optional :: check_sanity
logical err_flag, found_one

! Some init

lat%ramper_slave_bookkeeping = stale$
control = 0
ic = 0
control_to_ic = 0
iz = 0
allocate (ibr(0:ubound(lat%branch, 1)))

! Convert from old style ele%key = -1 to mark element removal to ele%ix_ele = -1
! Old style was bad since it complicated handling of super_slaves

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 0, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%key == -1) ele%ix_ele = -1
    if (ele%slave_status == super_slave$) ele%ix_ele = ie   ! Removal is determined by the lord's setting.
  enddo
enddo

! Mark super_slave elements according to the super_lord setting.
! First take care of overlapping super_lord_elements.

do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  call retain_overlapping_lords (lord)
enddo

do ie = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ie)
  if (lord%lord_status /= super_lord$) cycle
  do ie2 = 1, lord%n_slave
    slave => pointer_to_slave(lord, ie2)
    ib = slave%ix_branch
    if (lat%branch(ib)%ix_branch == -1) lord%ix_ele = -1
    if (lord%ix_ele == -1) then
      slave%ix_ele = -1
    elseif (slave%ix_ele == -1) then
      slave%ix_ele = int_garbage$   ! Will be corrected below.
    endif
  enddo
enddo

! Correct ele%ix_ele.
! Also a branch will be removed if all elements with index above 0 are marked for removal.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  do ie = 0, branch%n_ele_max
    ele => branch%ele(ie)
    if (ele%ix_ele == int_garbage$) ele%ix_ele = ie
  enddo
  if (all(branch%ele(1:branch%n_ele_track)%ix_ele == -1)) branch%ix_branch = -1
enddo

! If all elements in all branches are to be removed, keep branch 0 and element 0 of branch 0.

if (all(lat%branch%ix_branch == -1)) lat%branch(0)%ix_branch = 0

! Mark elements as dead if in a dead branch.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (branch%ix_branch == -1) then
    n = branch%n_ele_track
    branch%ele(0:n)%ix_ele = -1
  else
    branch%ele(0)%ix_ele = 0   ! Make sure not removed.
    iz = iz + 1
  endif

  n = lat%branch(ib)%n_ele_max
  if (ib /= 0 .and. iz == 1) n = n + lat%n_ele_max - lat%n_ele_track
  allocate (ibr(ib)%loc(0:n))
enddo

! Mark entries in control and ic arrays for deletion.

do i = 1, lat%n_ic_max
  control_to_ic(lat%ic(i)) = i
enddo

do i = 1, lat%n_control_max
  ctl => lat%control(i)
  slave => pointer_to_ele(lat, ctl%slave)
  if (slave%ix_ele == -1) control(i) = -1
  lord => pointer_to_ele(lat, ctl%lord)
  if (lord%ix_ele == -1 .or. lord%ix_ele == -2) control(i) = -1
  if (ctl%ix_attrib == int_garbage$ .or. ctl%attribute == 'REMOVE') control(i) = -1
  if (control(i) == -1 .and. control_to_ic(i) > 0) ic(control_to_ic(i)) = -1
enddo

do i = 1, lat%n_ic_max
  if (control(lat%ic(i)) == -1) ic(i) = -1
  ! Can happen for drifts that have been turned into lord null_ele elements due to superposition.
  ! In this case the lord will point to the element that is put in the space of the drift but,
  ! to simplify matters, the "slave" element will not point back to the null_ele.
  if (lat%ic(i) == 0) ic(i) = -1
enddo

! Compress branch%ele(:) array and fill in ibr(:)%loc(:) which maps old to new positions.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%ix_from_ele > -1) branch%ix_from_ele = ibr(branch%ix_from_branch)%loc(branch%ix_from_ele)%new_ix_ele
  if (branch%ix_from_ele < 0) then
    if (branch%ix_to_ele == 0) branch%ele(0)%value(inherit_from_fork$) = false$
    branch%ix_from_branch = -1
    branch%ix_to_ele = -1
  endif

  i2 = -1
  do i = 0, branch%n_ele_max
    ele => branch%ele(i)

    if (ele%ix_ele == -1 .or. ele%ix_ele == -2) then
      ibr(ib)%loc(i)%new_ix_ele    = -1
      ibr(ib)%loc(i)%new_ix_branch = -1
    else
      i2 = i2 + 1
      ibr(ib)%loc(i)%new_ix_ele    = i2
      ibr(ib)%loc(i)%new_ix_branch = ib
      ibr(ib)%loc(i2)%old_ix_ele    = i
      ibr(ib)%loc(i2)%old_ix_branch = ib
      if (i2 /= i) then
        call ele_equals_ele (branch%ele(i2), ele, .false.) ! Note: Nametable will be updated at end.
        branch%ele(i2)%ix_ele = i2
      endif
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

! Move lords if branch 0 is slated for deletion

if (lat%branch(0)%ix_branch == -1) then
  do ib2 = 1, ubound(lat%branch, 1)
    new_lord_branch => lat%branch(ib2)
    if (new_lord_branch%ix_branch > -1) exit
  enddo

  branch => lat%branch(0)
  call allocate_lat_ele_array(lat, new_lord_branch%n_ele_max + (branch%n_ele_max - branch%n_ele_track), ib2)

  ie2 = new_lord_branch%n_ele_track
  do ie = branch%n_ele_track+1, branch%n_ele_max
    ele => branch%ele(ie)
    ie2 = ie2 + 1
    new_lord_branch%ele(ie2) = ele
    j = ibr(0)%loc(ie)%old_ix_ele
    ibr(0)%loc(j)%new_ix_ele = ie2
    ibr(0)%loc(j)%new_ix_branch = ib2
    ibr(ib2)%loc(ie2)%old_ix_ele = j
    ibr(ib2)%loc(ie2)%old_ix_branch = 0
    call init_ele(branch%ele(ie), ix_ele = ie, branch = branch)
  enddo
  new_lord_branch%n_ele_max = ie2
endif

! Remove unwanted branches

ib2 = -1
do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (branch%ix_branch == -1) then
    if (associated(branch%wall3d)) call unlink_wall3d(branch%wall3d)
    cycle
  endif
  ib2 = ib2 + 1
  if (ib2 == ib) cycle
  branch2 => lat%branch(ib2)

  call deallocate_ele_array_pointers(branch2%ele)
  call transfer_branch (branch, branch2)
  nullify(branch%ele, branch%wall3d)
  branch2%ix_branch = ib2
  branch2%ele%ix_branch = ib2

  do ie = 0, branch%n_ele_max
    ie0 = ibr(ib)%loc(ie)%old_ix_ele
    ib0 = ibr(ib)%loc(ie)%old_ix_branch
    ibr(ib0)%loc(ie0)%new_ix_branch = ib2
  enddo

  do j = 0, ubound(lat%branch, 1)
    if (lat%branch(j)%ix_from_branch == ib) lat%branch(j)%ix_from_branch = ib2
  enddo
enddo

if (ib2 /= ib) call allocate_branch_array (lat, ib2)

! Compress lat%control() array and correct %lord and %slave pointers.

i2 = 0
do i = 1, lat%n_control_max
  if (control(i) == -1) cycle
  i2 = i2 + 1
  control(i) = i2
  ctl0 = lat%control(i)
  lat%control(i2) = ctl0
  lat%control(i2)%lord%ix_ele     = ibr(0)%loc(ctl0%lord%ix_ele)%new_ix_ele
  lat%control(i2)%slave%ix_ele    = ibr(ctl0%slave%ix_branch)%loc(ctl0%slave%ix_ele)%new_ix_ele
  lat%control(i2)%slave%ix_branch = ibr(ctl0%slave%ix_branch)%loc(ctl0%slave%ix_ele)%new_ix_branch
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
    if (control(i1) == i1 .and. control(i2) == i2) cycle  ! Nothing to do if control info has not changed

    ele%ix1_slave = 0  ! Start with no slaves
    ele%n_slave = 0
    ele%n_slave_field = 0
    do j = i1, i2
      if (control(j) == -1) cycle
      if (ele%ix1_slave == 0) ele%ix1_slave = control(j)
      ctl => lat%control(control(j))
      if (ctl%attribute == 'FIELD_OVERLAPS') then
        ele%n_slave_field = ele%n_slave_field + 1
      else
        ele%n_slave = ele%n_slave + 1
      endif
    enddo
  enddo

enddo

! Correct lord info.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  if (branch%ix_to_ele > 0) branch%ix_to_ele = ibr(branch%ix_branch)%loc(branch%ix_to_ele)%new_ix_ele
  if (branch%ix_to_ele < 0) branch%ix_to_ele = -1

  do i = 0, branch%n_ele_max
    ele => branch%ele(i)

    ! Update fork element

    if (ele%key == fork$ .or. ele%key == photon_fork$) then
      ie0 = nint(ele%value(ix_to_element$))
      ib0 = nint(ele%value(ix_to_branch$))
      if (ibr(ib0)%loc(ie0)%new_ix_ele == -1) then
        ele%key = marker$
      else
        ele%value(ix_to_element$) = ibr(ib0)%loc(ie0)%new_ix_ele
        ele%value(ix_to_branch$) = ibr(ib0)%loc(ie0)%new_ix_branch
      endif
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
      if (ic(j) == -1) cycle
      if (ele%ic1_lord == 0) ele%ic1_lord = ic(j)
      icon = lat%ic(ic(j))
      ctl => lat%control(icon)
      if (ctl%attribute == 'FIELD_OVERLAPS') then
        ele%n_lord_field = ele%n_lord_field + 1
      else
        ele%n_lord = ele%n_lord + 1
      endif
    enddo

    ! Correct slave_status

    ele%slave_status = free$
    do j = 1, ele%n_lord
      lord => pointer_to_lord(ele, j)
      select case(lord%lord_status)
      case (group_lord$, overlay_lord$, girder_lord$)
        if (ele%slave_status == free$) ele%slave_status = minor_slave$
      case (multipass_lord$)
        ele%slave_status = multipass_slave$
      case (super_lord$)
        ele%slave_status = super_slave$
      end select
    enddo

  enddo

enddo

! Sanity check

call create_lat_ele_nametable(lat, lat%nametable)
call set_flags_for_changed_attribute(lat)
if (logic_option(.true., check_sanity)) call lat_sanity_check (lat, err_flag)

!---------------------------------------------------------------------------
contains

recursive subroutine retain_overlapping_lords (lord)

type (ele_struct) lord
type (ele_struct), pointer :: slave, lord2

integer ie, ie2

!

if (lord%lord_status /= super_lord$) return
if (lord%ix_ele == -1) return

do ie = 1, lord%n_slave
  slave => pointer_to_slave(lord, ie)
  do ie2 = 1, slave%n_lord
    lord2 => pointer_to_lord(slave, ie2)
    if (lord2%lord_status /= super_lord$) cycle
    if (lord2%ix_ele /= -1) cycle
    lord2%ix_ele = int_garbage$  ! Will recover true %ix_ele later.
    call retain_overlapping_lords (lord2)
  enddo
enddo

end subroutine retain_overlapping_lords

end subroutine
          
