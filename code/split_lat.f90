!+
! Subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, check_controls, save_null_drift, err_flag)
!
! Subroutine to split a lat at a point. Subroutine will not split the lat if the split
! would create a "runt" element with length less than bmad_com%significant_length.
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
!   s_split         -- Real(rp): Position at which lat is to be split.
!   add_suffix      -- Logical, optional: If True (default) add '#1' and '#2" suffixes
!                        to the split elements. 
!   check_controls  -- Logical, optional: If True (default) then call check_lat_controls
!                        after the split to make sure everything is ok.
!   save_null_drift -- Logical, optional: Save a copy of a drift to be split as a null_ele?
!                         This is useful if when superpositions are done. See add_superimpose for more info.
!                         Default is False.
!
! Output:
!   lat        -- lat_struct: Modified lat structure.
!   ix_split   -- Integer: Index of element just before the split.
!   split_done -- Logical: True if lat was split.
!   err_flag   -- Logical, optional: Set true if there is an error, false otherwise.
!-

subroutine split_lat (lat, s_split, ix_branch, ix_split, split_done, add_suffix, check_controls, save_null_drift, err_flag)

use bmad_struct
use bmad_interface, except_dummy => split_lat
use bookkeeper_mod, only: control_bookkeeper

implicit none

type (lat_struct), target :: lat
type (ele_struct), save :: ele
type (ele_struct), pointer :: ele1, ele2, slave, lord
type (branch_struct), pointer :: branch

real(rp) s_split, len_orig, len1, len2, coef1, coef2, coef_old, ds_fudge

integer i, j, k, ix, ix_branch
integer ix_split, ixc, ix_attrib, ix_super_lord
integer icon, ix2, inc, nr, n_ic2, ct

logical split_done, err
logical, optional :: add_suffix, check_controls, save_null_drift, err_flag

character(16) :: r_name = "split_lat"

! Check for s_split out of bounds.

if (present(err_flag)) err_flag = .true.

branch => lat%branch(ix_branch)
ds_fudge = bmad_com%significant_length

nr = branch%n_ele_track
if (s_split < branch%ele(0)%s - ds_fudge .or. s_split > branch%ele(nr)%s + ds_fudge) then
  call out_io (s_fatal$, r_name, 'POSITION OF SPLIT NOT WITHIN LAT: \es12.3\ ',  &
                                  r_array = [s_split] )
  if (bmad_status%exit_on_error) call err_exit
endif

! Find where to split.

do ix_split = 0, branch%n_ele_track
  if (abs(branch%ele(ix_split)%s - s_split) < 10*ds_fudge) then
    split_done = .false.
    if (present(err_flag)) err_flag = .false.
    return
  endif
  if (branch%ele(ix_split)%s > s_split) exit
enddo

split_done = .true.
ele = branch%ele(ix_split)
ele%lat => lat   ! So we can use pointer_to_lord
len_orig = ele%value(l$)
len2 = branch%ele(ix_split)%s - s_split
len1 = len_orig - len2

! there is a problem with custom elements in that we don't know which
! attributes (if any) scale with length.

if (ele%key == custom$ .or. ele%key == match$) then
  call out_io (s_fatal$, r_name, "I DON'T KNOW HOW TO SPLIT THIS ELEMENT:" // ele%name)
  if (bmad_status%exit_on_error) call err_exit
endif

! save element to be split as a null element if needed

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

! kill any talyor series, etc.
! Note that %a_pole and %b_pole components are the exception

call deallocate_ele_pointers (ele1)
call deallocate_ele_pointers (ele2)

ele1%lat => lat ! reinstate
ele2%lat => lat ! reinstate

if (associated(ele%a_pole)) then
  call multipole_init(ele1)
  call multipole_init(ele2)
  ele1%a_pole = ele%a_pole; ele1%b_pole = ele%b_pole
  ele2%a_pole = ele%a_pole; ele2%b_pole = ele%b_pole
endif

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
  n_ic2 = ixc + ele%n_lord
  ele1%ic1_lord = ixc + 1
  ele1%ic2_lord = n_ic2
  lat%n_ic_max = n_ic2

  do j = 1, ele%n_lord

    lord => pointer_to_lord(ele, j, icon)

    coef_old = lat%control(icon)%coef
    ix_attrib = lat%control(icon)%ix_attrib

    if (ele%slave_status == super_slave$ .or.  &
          ix_attrib == hkick$ .or. ix_attrib == vkick$) then
      coef1 = coef_old * len1 / len_orig
      coef2 = coef_old * len2 / len_orig
    else
      coef1 = coef_old
      coef2 = coef_old
    endif

    lat%control(icon)%coef = coef2

    lord%n_slave = lord%n_slave + 1
    call add_lattice_control_structs (lat, lord)

    ix2 = lord%ix2_slave
    lat%control(ix2)%ix_slave  = ix_split
    lat%control(ix2)%ix_branch = ix_branch
    lat%control(ix2)%ix_attrib = ix_attrib
    lat%control(ix2)%coef = coef1
    lat%ic(ixc+j) = ix2

    if (lord%lord_status == super_lord$) call order_super_lord_slaves (lat, lord%ix_ele)

  enddo

  goto 8000   ! and return

endif   ! split element is a super_slave

! Here if a free or overlay element
! Need to make a super lord to control the split elements.

call new_control (lat, ix_super_lord)
ele1 => branch%ele(ix_split)
ele2 => branch%ele(ix_split+1)
lat%n_ele_max = ix_super_lord
lat%ele(ix_super_lord) = ele
lat%ele(ix_super_lord)%lord_status = super_lord$
lat%ele(ix_super_lord)%value(l$) = len_orig
ixc = lat%n_control_max + 1
if (ixc+1 > size(lat%control)) call reallocate_control (lat, ixc+500)
lat%ele(ix_super_lord)%ix1_slave = ixc
lat%ele(ix_super_lord)%ix2_slave = ixc + 1
lat%ele(ix_super_lord)%n_slave = 2
lat%n_control_max = ixc + 1
lat%control(ixc)%ix_lord   = ix_super_lord
lat%control(ixc)%ix_slave  = ix_split
lat%control(ixc)%ix_branch = ix_branch
lat%control(ixc)%coef = len1 / len_orig
lat%control(ixc+1)%ix_lord   = ix_super_lord
lat%control(ixc+1)%ix_slave  = ix_split + 1
lat%control(ixc+1)%ix_branch = ix_branch
lat%control(ixc+1)%coef = len2 / len_orig

! overlay lord elements of the split element must now point towards the
! super lord

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  do k = lord%ix1_slave, lord%ix2_slave
    if (lat%control(k)%ix_slave == ix_split+1) then
      lat%control(k)%ix_slave  = ix_super_lord
      lat%control(k)%ix_branch = 0
    endif
  enddo
enddo

! split elements must now be pointing towards their lord

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

! last details:
!     1) Groups that point to the split element must be redirected to the lord
!     2) Call control_bookkeeper to remake the split elements
!     3) And return

8000  continue

do i = lat%n_ele_track+1, lat%n_ele_max
  ct = lat%ele(i)%lord_status
  if (ct /= group_lord$ .and. ct /= girder_lord$) cycle

  do k = 1, lat%ele(i)%n_slave
    slave => pointer_to_slave(lat%ele(i), k, j) 
    if (slave%ix_ele /= ix_split+1 .or. slave%ix_branch /= ix_branch) cycle
    if (lat%control(j)%ix_attrib == l$) then
      call out_io (s_warn$, r_name, 'GROUP: ' // lat%ele(i)%name, &
                                    'CONTROLS L$ OF SPLIT ELEMENT: ' // ele%name)
    elseif (ix_super_lord /= 0) then
      lat%control(j)%ix_slave  = ix_super_lord
      lat%control(j)%ix_branch = 0
    else
      call out_io (s_warn$, r_name, &
                    'GROUP: ' // lat%ele(i)%name, &
                    'CONTROLS SPLIT ELEMENT: ' // ele%name, &
                    'BUT NO LORD WAS MADE!')
      if (bmad_status%exit_on_error) call err_exit
    endif
  enddo

enddo

call control_bookkeeper (lat, ele1)
call control_bookkeeper (lat, ele2)

err = .false.  ! In case check_lat_controls is not called.
if (logic_option(.true., check_controls)) call check_lat_controls (lat, err)
if (present(err_flag)) err_flag = err

end subroutine
