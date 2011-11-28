!+
! Subroutine s_calc (lat)
!
! Subroutine to calculate the longitudinal distance S for the elements
! in a lattice.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lat -- lat_struct:
!
! Output:
!   lat -- lat_struct:
!-

subroutine s_calc (lat)

use bmad_struct
use bmad_interface, except_dummy => s_calc
use lat_ele_loc_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave
type (branch_struct), pointer :: branch

integer i, j, n, ic, icon, ix2
real(8) ss, s_end

! Just go through all the elements and add up the lengths.

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (.not. bmad_com%auto_bookkeeper .and. branch%param%status%length /= stale$) cycle
  ! Branches that branch from another branch start from zero
  if (branch%ix_from_branch > -1) branch%ele(0)%s = 0  
  ss = branch%ele(0)%s
  do n = 1, branch%n_ele_track
    ele => branch%ele(n)
    ss = ss + ele%value(l$)
    ele%s = ss
    if (ele%status%length == stale$) ele%status%length = ok$
  enddo
  branch%param%total_length = ss - branch%ele(0)%s
  branch%param%status%length = ok$
enddo


! Now fill in the s positions of the super_lords and zero everyone else.
! Exception: A null_ele lord element is the result of a superposition on a multipass section.
! We need to preserve the s value of this element.

do n = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(n)
  if (.not. bmad_com%auto_bookkeeper .and. lord%status%length /= stale$) cycle
  if (lord%key == null_ele$) cycle
  if (lord%n_slave == 0) cycle  ! Can happen when manipulating a lattice.
  if (lord%lord_status == super_lord$) then
    slave => pointer_to_slave(lord, lord%n_slave)
    lord%s = slave%s
  else
    lord%s = 0
  endif
  lord%status%length = ok$
enddo

end subroutine
