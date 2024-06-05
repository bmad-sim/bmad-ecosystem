!+
! Subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, 
!                               check_sanity, save_null_drift, err_flag, choose_max, ix_insert)
!
! Routine to split an element of the lattice into two to create a lattice that has an element boundary at the point s = s_split. 
! This routine will not split the lattice if the split would create a "runt" element with length less 
! than 5*bmad_com%significant_length.
!
! split_lat will create a super_lord element if needed and will redo the 
! appropriate bookkeeping for lords and slaves.
!
! Note: split_lat does NOT call make_mat6. The Twiss parameters are also not recomputed.
!
! Input:
!   lat             -- lat_struct: Original lat structure.
!   s_split         -- real(rp): Position at which lat%branch(ix_branch) is to be split.
!   ix_branch       -- integer: Index of lat%branch(:) to use.
!   add_suffix      -- logical, optional: If True (default) add '#1' and '#2" suffixes
!                        to the split elements. 
!   check_sanity    -- logical, optional: If True (default) then call lat_sanity_check
!                        after the split to make sure everything is ok.
!   save_null_drift -- logical, optional: Save a copy of a drift to be split as a null_ele?
!                         This is useful when superpositions are done. See add_superimpose for more info.
!                         Default is False.
!   choose_max      -- logical, optional: If no splitting of an element is needed, that is, s_split is at an element 
!                       boundary, there can be multiple possible values for ix_split if there exist zero length elements 
!                       at the split point. If choose_max = True, ix_split will be chosen to be the maximum possible 
!                       index and if choose_max = False ix_split will be chosen to be the minimal possible index.
!                       If s_split is not at an element boundary, the setting of choose_max is immaterial.
!                       If ix_insert is present, the default value of choose_max is set to give the closest element to ix_insert.
!                       If ix_insert is not present, the default value of choose_max is False.
!   ix_insert       -- integer, optional: Element index near the point to be split. ix_insert is useful in the case where
!                       there is a patch with a negative length which can create an ambiguity as to where to do the split
!                       In this case ix_insert will remove the ambiguity. Also useful to ensure where to split if there
!                       are elements with zero length nearby. Ignored if negative.
!
! Output:
!   lat           -- lat_struct: Modified lat structure.
!   ix_split      -- integer: Index of element just before the s = s_split point.
!   split_done    -- logical: True if lat was split.
!   err_flag      -- logical, optional: Set true if there is an error, false otherwise.
!-

subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, &
                            check_sanity, save_null_drift, err_flag, choose_max, ix_insert)

use bmad_interface, except_dummy => split_lat

implicit none

type (lat_struct), target :: lat
type (ele_struct), save :: ele
type (ele_struct), pointer :: ele1, ele2, slave, lord, super_lord
type (branch_struct), pointer :: branch, br
type (control_struct), pointer :: ctl

real(rp) s_split, len_orig, len1, len2, ds_fudge
real(rp) dl, w_inv(3,3)

integer, optional :: ix_insert
integer i, j, k, ix, ix_branch, ib, ie, n_slave
integer ix_split, ixc, ix_attrib, ix_super_lord, ix_start
integer ix2, inc, nr, n_ic2, ct

logical split_done, err, controls_need_removing, found, ix_insert_present
logical, optional :: add_suffix, check_sanity, save_null_drift, err_flag, choose_max

character(*), parameter :: r_name = "split_lat"

! Check for s_split out of bounds.

if (present(err_flag)) err_flag = .true.

lat%ramper_slave_bookkeeping = stale$
branch => lat%branch(ix_branch)
ds_fudge = bmad_com%significant_length

nr = branch%n_ele_track
if (s_split < branch%ele(0)%s - ds_fudge .or. s_split > branch%ele(nr)%s + ds_fudge) then
  call out_io (s_fatal$, r_name, 'POSITION OF SPLIT NOT WITHIN LAT: \es12.3\ ', r_array = [s_split])
  if (global_com%exit_on_error) call err_exit
endif

! Find where to split.

if (integer_option(-1, ix_insert) >= 0) then
  ix_start = ix_insert
  ix_insert_present = .true.
elseif (logic_option(.false., choose_max)) then
  ix_start = branch%n_ele_track
  ix_insert_present = .false.
else
  ix_start = 0
  ix_insert_present = .false.
endif

found = .false.

if (branch%ele(ix_start)%s < s_split + 5*ds_fudge .or. &
          (abs(branch%ele(ix_start)%s - s_split) < 5*ds_fudge .and. logic_option(.false., choose_max))) then
  do ix_split = ix_start, branch%n_ele_track
    if (abs(branch%ele(ix_split)%s - s_split) < 5*ds_fudge) then
      if (.not. logic_option(.false., choose_max) .or. ix_split == branch%n_ele_track .or. ix_insert_present) then
        split_done = .false.
        if (present(err_flag)) err_flag = .false.
        return
      else ! choose_max = True so check if next boundary is the correct choice
        found = .true.
      endif
    ! If we have gone past
    elseif (branch%ele(ix_split)%s >= s_split + 5*ds_fudge) then
      if (found) then
        split_done = .false.
        if (present(err_flag)) err_flag = .false.
        return
      else
        exit ! Element split needed.
      endif
    endif
  enddo

else
  do ix_split = ix_start, 0, -1
    if (abs(branch%ele(ix_split)%s - s_split) < 5*ds_fudge) then
      if (logic_option(.false., choose_max) .or. ix_split == 0 .or. ix_insert_present) then
        split_done = .false.
        if (present(err_flag)) err_flag = .false.
        return
      else ! choose_max = False so check if next boundary is the correct choice
        found = .true.
      endif
    ! If we have gone past 
    elseif (branch%ele(ix_split)%s_start <= s_split - 5*ds_fudge) then
      if (found) then
        split_done = .false.
        if (present(err_flag)) err_flag = .false.
        return
      else
        exit ! Element split needed.
      endif
    endif
  enddo
endif

! Here if split is to be done.

split_done = .true.
ele = branch%ele(ix_split)
ele%branch => branch   ! So we can use pointer_to_lord
len_orig = ele%value(l$)
len2 = branch%ele(ix_split)%s - s_split
len1 = len_orig - len2

! There is a problem with custom elements in that we don't know which attributes (if any) scale with length.
! Also splitting a Taylor element is not currently permitted.
! In the future this might change but the bookkeeping changes are non-negligible.

if (ele%key == custom$ .or. ele%key == match$ .or. ele%key == taylor$) then
  call out_io (s_fatal$, r_name, "I DON'T KNOW HOW TO SPLIT THIS ELEMENT:" // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Save element to be split as a null element if it is a drift that has not been already split.

if (branch%ele(ix_split)%key == drift$ .and. branch%ele(ix_split)%value(split_id$) == 0 .and. &
                                                                 logic_option(.false., save_null_drift)) then
  call new_control (lat, ix, branch%ele(ix_split)%name)
  lord => lat%ele(ix)
  lord = branch%ele(ix_split)
  lord%key = null_ele$
  lord%sub_key = drift$  ! To mark that the element was formally a drift
  lord%value(ix_branch$) = branch%ix_branch  ! Save branch index.
endif

! Insert a new element.
! Note: Any lat%control()%ix_ele pointing to ix_split will now point to ix_split+1.

call insert_element (lat, ele, ix_split, ix_branch)

ele1 => branch%ele(ix_split)
ele2 => branch%ele(ix_split+1)

if (logic_option(.true., add_suffix)) then
  ix = len_trim(ele%name)
  call set_ele_name (ele1, ele%name(:ix) // '#1')
  call set_ele_name (ele2, ele%name(:ix) // '#2')
endif

! The split elements have blank "alias" and "type" values so as to not to pollute the results of
! searching by alias or type

ele1%type = '';  ele1%alias = ''
ele2%type = '';  ele2%alias = ''

! split_id is used by bmad_parser to keep track of which elements where originally one bigger element

if ((ele1%key == drift$ .or. ele1%key == instrument$ .or. ele1%key == monitor$ .or. ele1%key == pipe$) &
                                                                       .and. ele1%value(split_id$) == 0) then
  ele1%value(split_id$) = s_split 
  ele2%value(split_id$) = s_split 
endif

! Kill any talyor series, etc.
! Note that %a_pole and %b_pole components are the exception

call deallocate_ele_pointers (ele1, nullify_branch = .false., dealloc_poles = .false.)
call deallocate_ele_pointers (ele2, nullify_branch = .false., dealloc_poles = .false.)

! put in correct lengths and s positions

ele1%value(l$) = len1
ele1%s = s_split

ele2%value(l$) = len2
ele2%s_start = s_split

!-------------------------------------------------------------
! Now to correct the slave/lord bookkeeping...

ix_super_lord = 0   ! no super lord made yet.
controls_need_removing = .false.

! A free drift needs nothing more.

if (ele%key == drift$ .and. ele%n_lord == 0) goto 8000

! If we have split a super_slave we need to make a 2nd control list for one
! of the split elements (can't have both split elements using the same list).
! Also: Redo the control list for the lord elements.

if (ele%slave_status == super_slave$) then
  if (ele%n_lord == 0) goto 8000  ! nothing to do for free element

  ixc = lat%n_ic_max
  n_ic2 = ixc 

  do j = 1, ele%n_lord
    ! If lord does not overlap ele1 then adjust padding and do not add
    ! as a lord to ele1

    lord => pointer_to_lord(ele, j, ctl)
    ix_attrib = ctl%ix_attrib

    if (.not. has_overlap(ele1, lord, branch)) then
      lord%value(lord_pad1$) = lord%value(lord_pad1$) - ele1%value(l$)
      cycle
    endif

    !

    call add_lattice_control_structs (lord, n_add_slave = 1)

    n_ic2 = n_ic2 + 1

    ix2 = lord%ix1_slave+lord%n_slave-1
    lat%control(ix2)%slave = lat_ele_loc_struct(ix_split, ix_branch)
    lat%control(ix2)%ix_attrib = ix_attrib
    lat%ic(n_ic2) = ix2

    ele1%ic1_lord = ixc + 1
    lat%n_ic_max = n_ic2

    if (lord%lord_status == super_lord$) call order_super_lord_slaves (lat, lord%ix_ele)
  enddo

  ! Remove lord/slave control for ele2 if the lord does not overlap

  do j = 1, ele2%n_lord
    lord => pointer_to_lord(ele2, j, ctl)
    if (has_overlap(ele2, lord, branch)) cycle
    ctl%attribute = 'REMOVE'  ! Mark for deletion
    lord%value(lord_pad2$) = lord%value(lord_pad2$) - ele2%value(l$)
    controls_need_removing = .true.
  enddo

  goto 8000   ! and return
endif   ! split element is a super_slave

! Here if split element is not a super_slave.
! Need to make a super lord to control the split elements.

call new_control (lat, ix_super_lord, ele%name)
ele1 => branch%ele(ix_split)
ele2 => branch%ele(ix_split+1)
super_lord => lat%ele(ix_super_lord)
lat%n_ele_max = ix_super_lord
super_lord = ele
super_lord%lord_status = super_lord$
super_lord%value(l$) = len_orig
ixc = lat%n_control_max
n_slave = 2 + ele%n_slave_field
if (ixc+2+ele%n_slave_field > size(lat%control)) call reallocate_control (lat, ixc+ele%n_slave_field+100)
super_lord%ix1_slave = ixc + 1
super_lord%n_slave = 2
super_lord%n_slave_field = ele%n_slave_field
lat%n_control_max = ixc + n_slave
lat%control(ixc+1)%lord  = lat_ele_loc_struct(ix_super_lord, 0)
lat%control(ixc+1)%slave = lat_ele_loc_struct(ix_split, ix_branch)
lat%control(ixc+2)%lord  = lat_ele_loc_struct(ix_super_lord, 0)
lat%control(ixc+2)%slave = lat_ele_loc_struct(ix_split+1, ix_branch)
do i = 1, super_lord%n_slave_field
  lat%control(ixc+2+i)%lord = lat_ele_loc_struct(ix_super_lord, 0)
  ix = ele%ix1_slave+ele%n_slave+i-1
  lat%control(ixc+2+i)%slave = lat%control(ix)%slave
  lat%control(ix)%attribute = 'REMOVE'    ! Mark for deletion with remove_eles_from_lat
enddo

if (ele2%n_slave_field /= 0) controls_need_removing = .true.
ele2%n_slave_field = 0

! lord elements of the split element must now point towards the super_lord

do i = 1, ele%n_lord+ele%n_lord_field
  lord => pointer_to_lord(ele, i)
  do k = lord%ix1_slave, lord%ix1_slave+lord%n_slave+lord%n_slave_field-1
    if (lat%control(k)%slave%ix_ele == ix_split+1) then
      lat%control(k)%slave = lat_ele_loc_struct(ix_super_lord, 0)
    endif
  enddo
enddo

! point split elements towards their lord

if (lat%n_ic_max+2 > size(lat%ic)) call reallocate_control (lat, lat%n_ic_max+100)

ele1%slave_status = super_slave$
inc = lat%n_ic_max + 1
ele1%ic1_lord = inc
ele1%n_lord = 1
lat%n_ic_max = inc
lat%ic(inc) = ixc + 1

ele2%slave_status = super_slave$
inc = lat%n_ic_max + 1
ele2%ic1_lord = inc
ele2%n_lord = 1
lat%n_ic_max = inc
lat%ic(inc) = ixc + 2

!---------------------------------------
! Last details...

8000  continue

if (controls_need_removing) call remove_eles_from_lat(lat, .false.)

! The length of a Patch element is a dependent attribute so adjust offsets and pitches accordingly.

if (ele1%key == patch$) then
  lord => pointer_to_lord (ele1, 1)
  ele1 => pointer_to_slave (lord, 1)
  call floor_angles_to_w_mat (lord%value(x_pitch$), lord%value(y_pitch$), lord%value(tilt$), w_mat_inv = w_inv)
  dl = lord%value(l$) - ele1%value(l$)
  ele1%value(x_offset$)     = lord%value(x_offset$) - dl * w_inv(3,1)
  ele1%value(y_offset$)     = lord%value(y_offset$) - dl * w_inv(3,2)
  ele1%value(z_offset$)     = lord%value(z_offset$) - dl * w_inv(3,3)
  ele1%value(t_offset$)     = lord%value(t_offset$)
  ele1%value(e_tot_offset$) = lord%value(e_tot_offset$)
  ele1%value(e_tot_set$)    = lord%value(e_tot_set$)
  ele1%value(p0c_set$)      = lord%value(p0c_set$)

  do i = 2, lord%n_slave
    ele2 => pointer_to_slave (lord, i)
    ele2%value(x_pitch$)      = 0
    ele2%value(y_pitch$)      = 0
    ele2%value(tilt$)         = 0
    ele2%value(x_offset$)     = 0
    ele2%value(y_offset$)     = 0
    ele2%value(z_offset$)     = ele2%value(l$)
    ele2%value(t_offset$)     = 0
    ele2%value(e_tot_offset$) = 0
    ele1%value(e_tot_set$)    = 0
    ele1%value(p0c_set$)      = 0
  enddo
endif

!

ele1%bookkeeping_state%attributes = stale$
ele1%bookkeeping_state%floor_position = stale$
ele2%bookkeeping_state%attributes = stale$
ele2%bookkeeping_state%floor_position = stale$

if (ix_super_lord == 0) then
  call control_bookkeeper (lat, ele1)
  call control_bookkeeper (lat, ele2)
else
  super_lord%bookkeeping_state%control = stale$
  call control_bookkeeper (lat, super_lord)
endif

err = .false.  ! In case lat_sanity_check is not called.
if (logic_option(.true., check_sanity)) call lat_sanity_check (lat, err)
if (present(err_flag)) err_flag = err

!--------------------------------------------------------------
contains

function has_overlap (slave, lord, branch) result (overlap)

type (ele_struct) slave, lord
type (branch_struct) branch

real (rp) s0_lord, s0_slave
logical overlap

!

overlap = .true.
if (lord%lord_status /= super_lord$) return

s0_lord = lord%s_start + bmad_com%significant_length
s0_slave = slave%s_start - bmad_com%significant_length

! Case where the lord does not wrap around the lattice

if (s0_lord >= branch%ele(0)%s) then
  if (slave%s < s0_lord .or. s0_slave > lord%s) overlap = .false.

! Case where the lord does wrap

else
  if (slave%s < s0_lord .and. s0_slave > lord%s) overlap = .false.
endif

end function has_overlap

end subroutine
