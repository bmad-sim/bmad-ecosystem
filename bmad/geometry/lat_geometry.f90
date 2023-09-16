!+
! Subroutine lat_geometry (lat)
!
! Routine to calculate the physical placement of all the elements in a lattice.
! That is, the layout on the floor. This is the same as the MAD convention.
!
! Note: This routine does NOT update %ele(i)%s. To do this call s_calc.
!
! Input:
!   lat -- lat_struct: The lattice.
!     %ele(0)%floor  -- Floor_position_struct: The starting point for the calculations.
!
! Output:
!   lat -- lat_struct: The lattice.
!     %ele(i)%floor --  floor_position_struct: Floor position.
!-

subroutine lat_geometry (lat)

use bmad_interface, dummy => lat_geometry

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave, b_ele, ele0, ele2
type (branch_struct), pointer :: branch
type (floor_position_struct) :: floor_dum = floor_position_struct(vec3_zero$, mat3_unit$, 0.0_rp, 0.0_rp, 0.0_rp)

integer i, i2, ix2, ie, ib, ie0
logical is_stale

character(16), parameter :: r_name = 'lat_geometry'

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (bmad_com%auto_bookkeeper) then
    is_stale = .true.
  else
    if (branch%param%bookkeeping_state%floor_position /= stale$) cycle
    is_stale = .false.
  endif

  ! If there are fiducial elements then survey the fiducial regions

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    if (ele%key /= fiducial$) cycle
    call ele_geometry (floor_dum, ele)

    ! It is important to do the preceding elements first in case the following elements includes a patch.

    branch%ele(i-1)%floor = ele%floor  ! Save time

    do i2 = i-1, 1, -1
      ele2 => branch%ele(i2)
      if (ele2%key == patch$ .and. is_true(ele2%value(flexible$))) exit
      if (ele2%key == fiducial$) then
        call out_io (s_fatal$, r_name, 'FIDUCIAL ELEMENTS IN A BRANCH MUST BE SEPARATED BY A FLEXIBLE PATCH')
        if (global_com%exit_on_error) call err_exit
        exit
      endif
      call ele_geometry (ele2%floor, ele2, branch%ele(i2-1)%floor, -1.0_rp)
      branch%ele(i2-1)%bookkeeping_state%floor_position = ok$
    enddo

    do i2 = i+1, branch%n_ele_track
      ele2 => branch%ele(i2)
      if (ele2%key == patch$ .and. is_true(ele2%value(flexible$))) exit
      if (ele2%key == fiducial$) then
        call out_io (s_fatal$, r_name, 'FIDUCIAL ELEMENTS IN A BRANCH MUST BE SEPARATED BY A FLEXIBLE PATCH')
        if (global_com%exit_on_error) call err_exit
        exit
      endif
      call ele_geometry (branch%ele(i2-1)%floor, ele2)
    enddo
  enddo

  ! Transfer info from the from_branch element if that element exists.

  if (branch%ix_from_branch > -1 .and. (is_stale .or. branch%ele(0)%bookkeeping_state%floor_position == stale$)) then
    b_ele => pointer_to_ele (lat, branch%ix_from_ele, branch%ix_from_branch)
    ie0 = nint(b_ele%value(ix_to_element$))
    call ele_geometry (b_ele%floor, b_ele, branch%ele(ie0)%floor)
    branch%ele(ie0)%bookkeeping_state%floor_position = ok$
    is_stale = .true.
  else
    ie0 = 0
  endif

  if (branch%ele(ie0)%bookkeeping_state%floor_position == stale$) then
    branch%ele(ie0)%bookkeeping_state%floor_position = ok$
    if (ie0+1 <= branch%n_ele_track) branch%ele(ie0+1)%bookkeeping_state%floor_position = stale$
    if (ie0-1 >= 0)                  branch%ele(ie0-1)%bookkeeping_state%floor_position = stale$
  endif

  do i = ie0+1, branch%n_ele_track
    call propagate_geometry(branch%ele(i), 1, is_stale)
  enddo

  do i = ie0-1, 0, -1
    call propagate_geometry(branch%ele(i), -1, is_stale)
  enddo

  branch%param%bookkeeping_state%floor_position = ok$
enddo

! put info in super_lords and multipass_lords

lat%lord_state%floor_position = ok$
lat%param%bookkeeping_state%floor_position = ok$

if (bmad_com%auto_bookkeeper) lat%ele(lat%n_ele_track+1:lat%n_ele_max)%bookkeeping_state%floor_position = stale$

do i = lat%n_ele_track+1, lat%n_ele_max  
  lord => lat%ele(i)

  if (lord%bookkeeping_state%floor_position /= stale$) cycle
  lord%bookkeeping_state%floor_position = ok$

  if (lord%n_slave == 0) cycle

  select case (lord%lord_status)
  case (super_lord$)
    slave => pointer_to_slave(lord, lord%n_slave) ! Last slave is at exit end.
    lord%floor = slave%floor
  case (multipass_lord$)
    slave => pointer_to_slave(lord, 1)
    lord%floor = slave%floor
  case (girder_lord$)
    call girder_lord_geometry (lord)
  end select

enddo

!---------------------------------------------------
contains

subroutine propagate_geometry (ele, dir, is_stale)

type (ele_struct) ele
type (branch_struct), pointer :: branch
type (floor_position_struct) floor0

integer ie, dir, ix, ix_pass, n_links, k
logical is_stale

!

ie = ele%ix_ele
branch => ele%branch

if (ele%key == patch$ .and. is_true(ele%value(flexible$))) then
  if (dir == -1 .or. ie == branch%n_ele_track) then
    call out_io (s_fatal$, r_name, 'CONFUSION! PLEASE CONTACT DAVID SAGAN!')
    if (global_com%exit_on_error) call err_exit
    return
  endif
  ! If the position of the element just after a flexible patch is is_stale, 
  ! the patch geometry should be recomputed just to be on the safe side.
  if (branch%ele(ie+1)%bookkeeping_state%floor_position /= ok$) is_stale = .true.  
endif

!

if (.not. is_stale .and. ele%bookkeeping_state%floor_position /= stale$) return

floor0 = ele%floor

if (dir == 1) then
  call ele_geometry (branch%ele(ie-1)%floor, ele)
else
  call ele_geometry (branch%ele(ie+1)%floor, branch%ele(ie+1), ele%floor, -1.0_rp)
  ele%bookkeeping_state%floor_position = ok$
endif

is_stale = (.not. (ele%floor == floor0))

end subroutine propagate_geometry

!---------------------------------------------------
! contains

! When girders themselves have girder lords, must do computation in order: Lord before slave.

recursive subroutine girder_lord_geometry (ele)
type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (floor_position_struct) :: floor_dum = floor_position_struct(vec3_zero$, mat3_unit$, 0.0_rp, 0.0_rp, 0.0_rp)
integer i

!

do i = 1, ele%n_lord
  lord => pointer_to_lord(ele, i)
  if (lord%key /= girder$) cycle
  if (lord%bookkeeping_state%floor_position /= stale$) cycle
  call girder_lord_geometry (lord)
enddo

call ele_geometry (floor_dum, ele)

end subroutine girder_lord_geometry

end subroutine lat_geometry
