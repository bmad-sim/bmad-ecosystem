!+
! Subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, check_sanity, save_null_drift, err_flag)
!
! Subroutine to split a lat at a point. Subroutine will not split the lat if the split
! would create a "runt" element with length less than 5*bmad_com%significant_length.
!
! split_lat will create a super_lord element if needed and will redo the 
! appropriate bookkeeping for lords and slaves. 
!
! Note: split_lat does NOT call make_mat6. The Twiss parameters are also not recomputed.
!
! Modules Needed:
!   use bmad
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
!                         This is useful if when superpositions are done. See add_superimpose for more info.
!                         Default is False.
!
! Output:
!   lat        -- lat_struct: Modified lat structure.
!   ix_split   -- integer: Index of element just before the split.
!   split_done -- logical: True if lat was split.
!   err_flag   -- logical, optional: Set true if there is an error, false otherwise.
!-

subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, check_sanity, save_null_drift, err_flag)

use bmad_interface, except_dummy => split_lat
use bookkeeper_mod, only: control_bookkeeper
use geometry_mod, only: floor_angles_to_w_mat

implicit none

type (lat_struct), target :: lat
type (ele_struct), save :: ele
type (ele_struct), pointer :: ele1, ele2, slave, lord, super_lord
type (branch_struct), pointer :: branch, br

real(rp) s_split, len_orig, len1, len2, ds_fudge
real(rp) dl, w_inv(3,3)

integer i, j, k, ix, ix_branch, ib, ie
integer ix_split, ixc, ix_attrib, ix_super_lord
integer icon, ix2, inc, nr, n_ic2, ct

logical split_done, err, controls_need_removing
logical, optional :: add_suffix, check_sanity, save_null_drift, err_flag

character(16) :: r_name = "split_lat"

! Check for s_split out of bounds.

if (present(err_flag)) err_flag = .true.

branch => lat%branch(ix_branch)
ds_fudge = bmad_com%significant_length

nr = branch%n_ele_track
if (s_split < branch%ele(0)%s - ds_fudge .or. s_split > branch%ele(nr)%s + ds_fudge) then
  call out_io (s_fatal$, r_name, 'POSITION OF SPLIT NOT WITHIN LAT: \es12.3\ ',  &
                                  r_array = [s_split] )
  if (global_com%exit_on_error) call err_exit
endif

! Find where to split.

do ix_split = 0, branch%n_ele_track
  if (abs(branch%ele(ix_split)%s - s_split) < 5*ds_fudge) then
    split_done = .false.
    if (present(err_flag)) err_flag = .false.
    return
  endif
  if (branch%ele(ix_split)%s > s_split) exit
enddo

split_done = .true.
ele = branch%ele(ix_split)
ele%branch => branch   ! So we can use pointer_to_lord
len_orig = ele%value(l$)
len2 = branch%ele(ix_split)%s - s_split
len1 = len_orig - len2

! There is a problem with custom elements in that we don't know which
! attributes (if any) scale with length.

if (ele%key == custom$ .or. ele%key == match$) then
  call out_io (s_fatal$, r_name, "I DON'T KNOW HOW TO SPLIT THIS ELEMENT:" // ele%name)
  if (global_com%exit_on_error) call err_exit
endif

! Save element to be split as a null element if needed

if (branch%ele(ix_split)%key == drift$ .and. logic_option(.false., save_null_drift)) then
  call new_control (lat, ixc)
  lat%ele(ixc) = branch%ele(ix_split)
  lat%ele(ixc)%key = null_ele$
endif

! Insert a new element.
! Note: Any lat%control()%ix_ele pointing to ix_split will now 
!  point to ix_split+1.

call insert_element (lat, ele, ix_split, ix_branch)

ele1 => branch%ele(ix_split)
ele2 => branch%ele(ix_split+1)

if (logic_option(.true., add_suffix)) then
  ix = len_trim(ele%name)
  ele1%name = ele%name(:ix) // '#1'
  ele2%name = ele%name(:ix) // '#2'
endif

! Kill any talyor series, etc.
! Note that %a_pole and %b_pole components are the exception

call deallocate_ele_pointers (ele1, nullify_branch = .false., dealloc_poles = .false.)
call deallocate_ele_pointers (ele2, nullify_branch = .false., dealloc_poles = .false.)

! put in correct lengths and s positions

ele1%value(l$) = len1
ele1%s = s_split
ele2%value(l$) = len2

!-------------------------------------------------------------
! Now to correct the slave/lord bookkeeping...

ix_super_lord = 0   ! no super lord made yet.

! A free drift needs nothing more.

if (ele%key == drift$ .and. ele%slave_status == free$) goto 8000

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

    lord => pointer_to_lord(ele, j, icon)
    if (.not. has_overlap(ele1, lord)) then
      lord%value(lord_pad1$) = lord%value(lord_pad1$) - ele1%value(l$)
      cycle
    endif

    !

    n_ic2 = n_ic2 + 1

    ix_attrib = lat%control(icon)%ix_attrib

    lord%n_slave = lord%n_slave + 1
    call add_lattice_control_structs (lat, lord)

    ix2 = lord%ix2_slave
    lat%control(ix2)%slave = lat_ele_loc_struct(ix_split, ix_branch)
    lat%control(ix2)%ix_attrib = ix_attrib
    lat%ic(n_ic2) = ix2

    ele1%ic1_lord = ixc + 1
    ele1%ic2_lord = n_ic2
    lat%n_ic_max = n_ic2

    if (lord%lord_status == super_lord$) call order_super_lord_slaves (lat, lord%ix_ele)

  enddo

  ! Remove lord/slave control for ele2 if the lord does not overlap

  controls_need_removing = .false.
  do j = 1, ele2%n_lord
    lord => pointer_to_lord(ele2, j, ixc)
    if (has_overlap(ele2, lord)) cycle
    lat%control(ixc)%ix_attrib = int_garbage$
    lord%value(lord_pad2$) = lord%value(lord_pad2$) - ele2%value(l$)
    controls_need_removing = .true.
  enddo

  if (controls_need_removing) call remove_eles_from_lat(lat, .false.)

  goto 8000   ! and return

endif   ! split element is a super_slave

! Here if split element is not a super_slave.
! Need to make a super lord to control the split elements.

call new_control (lat, ix_super_lord)
ele1 => branch%ele(ix_split)
ele2 => branch%ele(ix_split+1)
super_lord => lat%ele(ix_super_lord)
lat%n_ele_max = ix_super_lord
super_lord = ele
super_lord%lord_status = super_lord$
super_lord%value(l$) = len_orig
ixc = lat%n_control_max + 1
if (ixc+1 > size(lat%control)) call reallocate_control (lat, ixc+500)
super_lord%ix1_slave = ixc
super_lord%ix2_slave = ixc + 1
super_lord%n_slave = 2
lat%n_control_max = ixc + 1
lat%control(ixc)%ix_lord   = ix_super_lord
lat%control(ixc)%slave = lat_ele_loc_struct(ix_split, ix_branch)
lat%control(ixc+1)%ix_lord   = ix_super_lord
lat%control(ixc+1)%slave = lat_ele_loc_struct(ix_split + 1, ix_branch)

! lord elements of the split element must now point towards the super lord

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  do k = lord%ix1_slave, lord%ix2_slave
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
ele1%ic2_lord = inc
ele1%n_lord = 1
lat%n_ic_max = inc
lat%ic(inc) = ixc

ele2%slave_status = super_slave$
inc = lat%n_ic_max + 1
ele2%ic1_lord = inc
ele2%ic2_lord = inc
ele2%n_lord = 1
lat%n_ic_max = inc
lat%ic(inc) = ixc + 1

8000  continue

! Last details...

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

function has_overlap (slave, lord) result (overlap)

type (ele_struct) slave, lord
real (rp) s0_lord, s0_slave
logical overlap

!

overlap = .true.
if (lord%lord_status /= super_lord$) return

s0_lord = lord%s - lord%value(l$) + bmad_com%significant_length
s0_slave = slave%s - slave%value(l$) - bmad_com%significant_length

! Case where the lord does not wrap around the lattice

if (s0_lord >= branch%ele(0)%s) then
  if (slave%s < s0_lord .or. s0_slave > lord%s) overlap = .false.

! Case where the lord does wrap

else
  if (slave%s < s0_lord .and. s0_slave > lord%s) overlap = .false.
endif

end function has_overlap

end subroutine
