module reverse_mod

use bookkeeper_mod

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lat_reverse (lat_in, lat_rev)
!
! Subroutine to construct a lat structure with the elements in reversed order.
! This may be used for backward tracking through the lat. 
!
! The correspondence between elements in the two lattices is as follows:
!     lat_rev%ele(lat%n_ele_track+1-i) = lat_in%ele(i)  For 0 < i <= lat%n_ele_track
!     lat_rev%ele(i)                   = lat_in%ele(i)  For lat%n_ele_track < i 
!
! The transformation from particle coordinates
! in one lat to the corrisponding coordinates in the other is:
!     (x, P_x, y, P_y, z, P_z) -> (x, -P_x, y, -P_y, -z, P_z)
!
! Note: The Twiss parameters will not be correct for the reversed lat.
! You will need to compute them.
!
! Modules needed:
!   use reverse_mod
!
! Input:
!   lat_in -- lat_struct: Input lat.
!
! Output:
!   lat_rev -- lat_struct: Lat with the elements in reversed order.
!               The lat_rev actual argument may not be the same as the lat_in actual argument.
!-

subroutine lat_reverse (lat_in, lat_rev)

implicit none

type (lat_struct), target :: lat_in, lat_rev
type (ele_struct), pointer :: lord, ele
type (control_struct), pointer :: con
type (branch_struct), pointer :: branch, branch_in

integer i, n, i1, i2, nr, n_con, ib
integer :: ix_con(size(lat_in%control))

logical err_flag

! Correct control information

lat_rev = lat_in

do i = 1, lat_rev%n_control_max
  con => lat_rev%control(i)
  if (con%ix_slave <= lat_rev%n_ele_track) con%ix_slave = lat_rev%n_ele_track+1-con%ix_slave
  if (con%ix_lord <= lat_rev%n_ele_track)  con%ix_lord  = lat_rev%n_ele_track+1-con%ix_lord
enddo

! Slaves of a super lord must be in assending sequence.
! ix_con keeps track of the switching.
! Also: adjust s-position of lords.

n_con = size(lat_in%control)
forall (i = 1:n_con) ix_con(i) = i 

do i = lat_rev%n_ele_track+1, lat_rev%n_ele_max
  lord => lat_rev%ele(i)
  if (lord%lord_status /= super_lord$) cycle
  i1 = lord%ix1_slave
  i2 = lord%ix2_slave
  lat_rev%control(i1:i2) = lat_rev%control(i2:i1:-1)
  ix_con(i1:i2) = ix_con(i2:i1:-1)
  if (lord%s > lord%value(l$)) then
    lord%s = lat_rev%param%total_length - (lord%s - lord%value(l$))
  else  ! Lord wraps around zero case
    lord%s = lord%value(l$) - lord%s
  endif
enddo

n = lat_rev%n_ic_max
lat_rev%ic(1:n) = ix_con(lat_rev%ic(1:n))

! Transfer info from lat_in to lat_rev.
! the lat lattice is used since the actual arguments of lat_in and lat_rev
! may be the same

do ib = 0, ubound(lat_in%branch, 1)

  branch => lat_rev%branch(ib)
  branch_in => lat_in%branch(ib)

  nr = branch%n_ele_track
  branch%ele(1:nr) = branch_in%ele(nr:1:-1)

  ! Flip longitudinal stuff, maps

  do i = 0, branch%n_ele_max
    ele => branch%ele(i)
    ele%orientation = -ele%orientation
    if (i <= nr) then
      ele%s = branch%param%total_length - (ele%s - ele%value(l$))
      ele%ref_time = branch_in%ele(nr)%ref_time - branch_in%ele(nr-i)%ref_time
    endif
  enddo

  ! Cleanup

  branch%param%t1_with_RF = 0  ! Init
  branch%param%t1_no_RF = 0    ! Init

enddo

! Finish

call lat_sanity_check (lat_rev, err_flag)
call set_ele_status_stale (lat_rev%ele(0), floor_position_group$)
call lattice_bookkeeper (lat_rev)

end subroutine lat_reverse

end module


