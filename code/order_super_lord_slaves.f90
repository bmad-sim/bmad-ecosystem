!+
! Subroutine order_super_lord_slaves (lat, ix_lord)
!
! Subroutine to make the slave elements of a super_lord in order.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat    -- lat_struct: Lat.
!   ix_lord -- Integer: Index of lord element.
!
! Output
!   lat -- lat_struct: Lat with fixed controls.
!-

subroutine order_super_lord_slaves (lat, ix_lord)

use bmad_struct
use bmad_interface, except_dummy => order_super_lord_slaves
use nr

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: lord, slave
type (control_struct), allocatable :: cs(:)

integer i, ix_lord, ix1, ix2, ns
integer, allocatable :: ixx(:), iyy(:)

real(rp) ds
real(rp), allocatable :: s_rel(:)

! Init setup.

lord => lat%ele(ix_lord)

if (lord%lord_status /= super_lord$) then
  print *, 'ERROR IN ORDER_SUPER_LORD_SLAVES: ELEMENT NOT A SUPER_LORD'
  call err_exit
endif

! Make an array of distances between the slave elements and the lord element.
! Note that all distances should be negative but since the lord can wrap
! around the ends of the lattice, we need to correct for that if that is the situation.


ns = lord%n_slave
allocate (s_rel(ns), ixx(ns), iyy(ns), cs(ns))

do i = 1, lord%n_slave
  slave => pointer_to_slave (lat, lord, i)
  ds = slave%s - lord%s
  if (ds > bmad_com%significant_longitudinal_length) ds = ds - lat%branch(slave%ix_branch)%param%total_length
  if (-ds > lord%value(l$)) then
    print *, 'ERROR IN ORDER_SUPER_LORD_SLAVES: INTERNAL ERROR! ', trim(slave%name)
    print *, slave%s, lord%s, ds, lord%value(l$)
    call err_exit
  endif
  s_rel(i) = ds
enddo

! Sort slaves by distance.

call indexx (s_rel, ixx)

ix1 = lord%ix1_slave; ix2 = lord%ix2_slave
cs = lat%control(ix1:ix2) 

do i = 1, ns
  lat%control(i+ix1-1) = cs(ixx(i))
  iyy(ixx(i)) = i
enddo

do i = 1, lat%n_ic_max
  if (lat%ic(i) >= ix1 .and. lat%ic(i) <= ix2) then
    lat%ic(i) = iyy(lat%ic(i)+1-ix1) + ix1 - 1
  endif
enddo

deallocate (s_rel, ixx, iyy, cs)

end subroutine
