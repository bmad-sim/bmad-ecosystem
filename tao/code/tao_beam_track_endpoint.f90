!+
! Function tao_beam_track_endpoint (ele_id, lat, branch_str, where, u) result (ele)
!
! Routine to point to the track endpoint element.
!
! Input:
!   ele_id      -- character(*): Name or index of the element.
!   lat         -- lat_struct: Lattice.
!   branch_str  -- integer: Branch where the tracking is done. '' => Branch not specified.
!   where       -- character(*): 'TRACK_END', 'TRACK_START', etc.. Used for error messages.
!   u           -- tao_universe_struct: Universe beam is being tracked in.
!
! Output:
!   ele         -- ele_struct, pointer: Pointer to the track endpoint element. Nullified if error.
!                   Note: blank ele_id is handled if "where" contains 'END' or 'START'
!-

function tao_beam_track_endpoint (ele_id, lat, branch_str, where, u) result (ele)

use tao_interface, dummy => tao_beam_track_endpoint

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele0, ele
type (ele_pointer_struct), allocatable, target :: eles(:)
type (branch_struct), pointer :: branch
type (tao_universe_struct), target :: u

integer ix_branch, n_loc, ie
logical err

character(*) ele_id, where, branch_str
character(*), parameter :: r_name = 'tao_beam_track_endpoint'

!

ele => null()
ix_branch = -1

branch => pointer_to_branch(branch_str, lat)
if (.not. associated(branch) .and. branch_str /= '') then
  call out_io (s_error$, r_name, trim(where) // ' BRANCH NOT FOUND: ' // branch_str)
  return
endif
if (associated(branch)) ix_branch = branch%ix_branch

if (ele_id == '') then
  if (.not. associated(branch)) branch => lat%branch(0)

  if (index(where, 'START') /= 0) then
    ele => branch%ele(0)

  elseif (index(where, 'END') /= 0) then
    if (branch%param%geometry == open$) then
      ele => branch%ele(branch%n_ele_track)
    else
      ie = u%model_branch(branch%ix_branch)%beam%ix_track_start - 1
      if (ie == -1) then
        ele => branch%ele(0)
      else
        ele => branch%ele(ie)
      endif
    endif

  else
    call out_io(s_error$, r_name, 'BLANK STRING FOR ' // trim(where) // ' ELEMENT!')
  endif
  return

else
  call lat_ele_locator (ele_id, lat, eles, n_loc, err, ix_dflt_branch = ix_branch)
  if (err .or. n_loc == 0) then
    if (ix_branch > -1) then
      call out_io (s_error$, r_name, trim(where) // ' ELEMENT NOT FOUND: ' // ele_id, &
                                              'IN BRANCH: ' // int_str(ix_branch))
    else
      call out_io (s_error$, r_name, trim(where) // ' ELEMENT NOT FOUND: ' // ele_id)
    endif
    return
  endif

  if (n_loc > 1) then
    call out_io (s_error$, r_name, 'MULTIPLE ' // trim(where) // ' ELEMENTS FOUND: ' // ele_id)
    return
  endif
endif


ele0 => eles(1)%ele
select case (ele0%lord_status)
case (multipass_lord$) 
  call out_io (s_error$, r_name, 'BEAM ' // trim(where) // ' ELEMENT IS A MULTIPASS LORD: ' // ele_id, &
                   'WHERE TO START IN THE LATTICE IS AMBIGUOUS SINCE WHICH PASS TO USE IS NOT SPECIFIED.')
  return
case (group_lord$, overlay_lord$)
  if (ele0%n_slave == 1) then
    ele0 => pointer_to_slave(ele0, 1)
  else
    call out_io (s_error$, r_name, 'BEAM ' // trim(where) // ' ELEMENT IS A CONTROLLER TYPE ELEMENT: ' // ele_id, &
                                   'THIS DOES NOT MAKE SENSE.')
    return
  endif
case (girder_lord$)
  call out_io (s_error$, r_name, 'BEAM ' // trim(where) // ' ELEMENT IS A GIRDER TYPE ELEMENT: ' // ele_id, &
                                 'THIS DOES NOT MAKE SENSE.')
  return
end select

if (ele0%lord_status == super_lord$) ele0 => pointer_to_slave(ele0, ele0%n_lord)

if (ele0%n_slave /= 0) then
  call out_io (s_error$, r_name, 'UNABLE TO ASSOCIATE BEAM ' // trim(where) // ' ELEMENT: ' // ele_id, &
                                 'TO SOMEPLACE IN THE TRACKING LATTICE.')
  return
endif

ele => ele0

if (ix_branch > -1 .and. ele%ix_branch /= ix_branch) then
  call out_io (s_error$, r_name, 'LATTICE ELEMENT: ' // ele_id, &
                                 'NOT IN CORRECT BRANCH: ' // branch_str)
  return
endif

end function tao_beam_track_endpoint
