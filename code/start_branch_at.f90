!+
! Subroutine start_branch_at (lat, ele_start, move_end_marker, error)
!
! Routine to, while preserving element order, rearrange the elements in a branch such that 
! the first element is ele_start. In the finished lattice, the original first element will
! be just after the original last element. That is, if you think of the branch elements
! as arranged in a circle, this routine shifts where the cut is made to form the branch element array.
!
! Input:
!   lat             -- lat_struct: Lattice to modify.
!   ele_start       -- character(*): Start element. Ele_start will identify the lattice branch to modify.
!   move_end_marker -- logical: If True then the end marker (if it is present) will be shifted like
!                       any other element. False means that the end marker will stay at the end.
!
! Output:
!   lat             -- lat_struct: Modified lattice.
!   error           -- logical: Set True if there is an error Set False if not.
!-

subroutine start_branch_at (lat, ele_start, move_end_marker, error)

use bmad, dummy => start_branch_at

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele0, ele1, e2_arr(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (control_struct), pointer :: ctl

integer ic, ie, n_loc, ie1, iet
logical move_end_marker, error, err, shift_end

character(*) ele_start
character(*), parameter :: r_name = 'start_branch_at'

!

error = .true.

call lat_ele_locator (ele_start, lat, eles, n_loc, err, above_ubound_is_err = .false.)
if (err) return

if (n_loc == 0) then
  call out_io (s_error$, r_name, 'NO LATTICE ELEMENTS FOUND: ' // ele_start, 'LATTICE NOT MODIFIED.')
  return
endif

if (n_loc > 1) then
  call out_io (s_error$, r_name, 'MULTIPLE ELEMENTS FOUND MATCHING: ' // ele_start, 'LATTICE NOT MODIFIED.')
  return
endif

ele1 => eles(1)%ele
if (ele1%n_slave > 0) ele1 => pointer_to_slave(ele1, 1)
branch => pointer_to_branch(ele1)

if (ele1%key == overlay$ .or. ele1%key == group$) then
  call out_io (s_error$, r_name, &
          'IT DOES NOT MAKE SENSE TO START BRANCH AT A GROUP OR OVERLAY ELEMENT: ' // ele_start, & 
          'LATTICE NOT MODIFIED.')
  return
endif

if (ele1%lord_status == multipass_lord$) then
  call out_io (s_error$, r_name, &
          'IT DOES NOT MAKE SENSE TO START BRANCH AT A MULTIPASS_LORD ELEMENT: ' // ele_start, & 
          'LATTICE NOT MODIFIED.')
  return
endif

ie1 = ele1%ix_ele
iet = branch%n_ele_track
shift_end = (move_end_marker .or. branch%ele(iet)%name /= 'END' .or. branch%ele(iet)%key /= marker$)
if (.not. shift_end) iet = iet - 1

!

e2_arr => branch%ele
nullify (branch%ele)
call allocate_lat_ele_array(lat, ubound(e2_arr, 1), branch%ix_branch)

do ie = 0, ubound(e2_arr, 1)
  if (ie == 0 .or. ie > iet) then
    branch%ele(ie) = e2_arr(ie)
  elseif (ie < ie1) then
    branch%ele(ie+1+iet-ie1) = e2_arr(ie)
  else
    branch%ele(ie+1-ie1) = e2_arr(ie)
  endif
enddo

do ic = 1, lat%n_control_max
  ctl => lat%control(ic)
  if (ctl%slave%ix_branch /= branch%ix_branch) cycle
  if (ctl%slave%ix_ele > branch%n_ele_track) cycle  ! slave is a lord so do nothing

  if (ctl%slave%ix_ele < ie1) then
    ctl%slave%ix_ele = ctl%slave%ix_ele + 1 + iet - ie1
  else
    ctl%slave%ix_ele = ctl%slave%ix_ele + 1 - ie1
  endif
enddo

! Transfer Twiss from first non-deleted element back to beginning element.

ele0 => e2_arr(ie1-1)
if (ele1%value(p0c$) > 0 .and. ele0%a%beta /= 0 .and. ele0%b%beta /= 0) then   ! Ref energy has been computed
  call transfer_twiss (ele0, branch%ele(0))
  branch%ele(0)%s_start             = ele0%s_start
  branch%ele(0)%s                   = ele0%s
  branch%ele(0)%floor               = ele0%floor
  branch%ele(0)%ref_time            = ele0%ref_time
  branch%ele(0)%value(e_tot$)       = ele0%value(e_tot$)
  branch%ele(0)%value(e_tot_start$) = branch%ele(0)%value(e_tot$)
  branch%ele(0)%value(p0c$)         = ele0%value(p0c$)
  branch%ele(0)%value(p0c_start$)   = branch%ele(0)%value(p0c$)
  branch%ele(0)%a%phi = 0
  branch%ele(0)%b%phi = 0
  branch%ele(0)%z%phi = 0
  call set_flags_for_changed_attribute(branch%ele(0), branch%ele(0)%value(p0c$))
endif

call deallocate_ele_array_pointers(e2_arr)

! Finish

call create_lat_ele_nametable(lat, lat%nametable)
call set_flags_for_changed_attribute(branch%ele(0))
call lattice_bookkeeper (lat)

error = .false.

end subroutine start_branch_at
