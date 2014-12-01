module reverse_mod

use bookkeeper_mod

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine lat_reverse (lat_in, lat_rev, track_antiparticle, make_mats6)
!
! Subroutine to construct a lat structure with the elements in reversed order
! and with the elements with reversed orientaiton.
! This may be used for backward tracking through the lat. 
!
! The correspondence between elements in the two lattices is as follows:
!     lat_rev%ele(lat%n_ele_track+1-i) = lat_in%ele(i)  For 0 < i <= lat%n_ele_track
!     lat_rev%ele(i)                   = lat_in%ele(i)  For lat%n_ele_track < i 
!
! Note: The Twiss parameters will not be correct for the reversed lat.
! You will need to compute them.
!
! Modules needed:
!   use reverse_mod
!
! Input:
!   lat_in     -- lat_struct: Input lat.
!   track_antiparticle 
!              -- logical, optional: Set the default type of particle to track in lat_rev
!                   to the antiparticle of the default in lat? Default is True.
!   make_mats6 -- logical, optional: Make the matrices for lat_rev? Default is True.
!
! Output:
!   lat_rev -- lat_struct: Lat with the elements in reversed order.
!-

subroutine lat_reverse (lat_in, lat_rev, track_antiparticle, make_mats6)

implicit none

type (lat_struct), target :: lat_in, lat_rev, rlat
type (ele_struct), pointer :: lord, ele
type (control_struct), pointer :: con
type (branch_struct), pointer :: branch, branch_in

integer i, n, i1, i2, nr, n_con, ib
integer :: ix_con(size(lat_in%control))

logical, optional :: track_antiparticle, make_mats6
logical err_flag

! Correct control information

rlat = lat_in  ! Use temp lattice in case lat_in and lat_rev are the same actual argument.

do i = 1, rlat%n_control_max
  con => rlat%control(i)
  if (con%ix_slave <= rlat%n_ele_track) con%ix_slave = rlat%n_ele_track+1-con%ix_slave
  if (con%ix_lord <= rlat%n_ele_track)  con%ix_lord  = rlat%n_ele_track+1-con%ix_lord
enddo

! Slaves of a super lord must be in assending sequence.
! ix_con keeps track of the switching.
! Also: adjust s-position of lords.

n_con = size(lat_in%control)
forall (i = 1:n_con) ix_con(i) = i 

do i = rlat%n_ele_track+1, rlat%n_ele_max
  lord => rlat%ele(i)
  if (lord%lord_status /= super_lord$) cycle
  i1 = lord%ix1_slave
  i2 = lord%ix2_slave
  rlat%control(i1:i2) = rlat%control(i2:i1:-1)
  ix_con(i1:i2) = ix_con(i2:i1:-1)
  if (lord%s > lord%value(l$)) then
    lord%s = rlat%param%total_length - (lord%s - lord%value(l$))
  else  ! Lord wraps around zero case
    lord%s = lord%value(l$) - lord%s
  endif
enddo

n = rlat%n_ic_max
rlat%ic(1:n) = ix_con(rlat%ic(1:n))

! Transfer info from lat_in to rlat.
! the lat lattice is used since the actual arguments of lat_in and rlat
! may be the same

do ib = 0, ubound(lat_in%branch, 1)

  branch => rlat%branch(ib)
  branch_in => lat_in%branch(ib)

  nr = branch%n_ele_track
  branch%ele(1:nr) = branch%ele(nr:1:-1)

  if (logic_option(.true., track_antiparticle)) branch%param%default_tracking_species = &
                                              antiparticle(branch%param%default_tracking_species)

  ! Flip longitudinal stuff, maps

  do i = 0, branch%n_ele_max
    ele => branch%ele(i)
    ele%ix_ele = i
    ele%orientation = -ele%orientation
    if (associated(ele%taylor(1)%term)) call kill_taylor(ele%taylor)
    if (associated(ele%ptc_genfield)) call kill_ptc_genfield(ele%ptc_genfield)
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

lat_rev = rlat
call deallocate_lat_pointers (rlat) 

call lat_sanity_check (lat_rev, err_flag)
call set_flags_for_changed_attribute(lat_rev)
call lattice_bookkeeper (lat_rev)

if (logic_option(.true., make_mats6)) call lat_make_mat6(lat_rev, ix_branch = -1)

end subroutine lat_reverse

end module


