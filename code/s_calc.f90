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

integer i, j, n, ic, icon, ix2, nt
real(8) ss, s_end

! Just go through all the elements and add up the lengths.

do i = 0, ubound(lat%branch, 1)
  branch => lat%branch(i)
  if (i > 0) branch%ele(0)%s = 0  ! Branches start from zero
  ss = branch%ele(0)%s
  nt = branch%n_ele_track
  do n = 1, nt
    ele => branch%ele(n)
    ss = ss + ele%value(l$)
    branch%ele(n)%s = ss
  enddo
enddo

lat%param%total_length = ss - lat%ele(0)%s

! Now fill in the s positions of the super_lords and zero everyone else.
! Exception: A null_ele lord element is the result of a superposition on a multipass section.
! We need to preserve the s value of this element.

do n = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(n)
  if (lord%key == null_ele$) cycle
  if (lord%n_slave == 0) cycle  ! Can happen when manipulating a lattice.
  if (lord%lord_status == super_lord$) then
    slave => pointer_to_slave (lat, lord, lord%n_slave)
    lord%s = slave%s
  else
    lord%s = 0
  endif
enddo

end subroutine
