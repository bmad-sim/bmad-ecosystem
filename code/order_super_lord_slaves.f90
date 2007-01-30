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

#include "CESR_platform.inc"

subroutine order_super_lord_slaves (lat, ix_lord)

  use bmad_struct
  use bmad_interface, except => order_super_lord_slaves
  use nr

  implicit none

  type (lat_struct), target :: lat
  type (ele_struct), pointer :: ele
  type (control_struct), allocatable :: cs(:)

  integer i, ix, ix_lord, ix1, ix2, ns
  integer, allocatable :: ixx(:), iyy(:)

  real(rp) ds
  real(rp), allocatable :: s_rel(:)

! Init setup.

  ele => lat%ele(ix_lord)
  ix1 = ele%ix1_slave; ix2 = ele%ix2_slave

  if (ele%control_type /= super_lord$) then
    print *, 'ERROR IN ORDER_SUPER_LORD_SLAVES: ELEMENT NOT A SUPER_LORD'
    call err_exit
  endif

! Make an array of distances between the slave elements and the lord element.
! Note that all distances are negative.

  ns = ele%n_slave
  allocate (s_rel(ns), ixx(ns), iyy(ns), cs(ns))

  do i = ix1, ix2
    ix = lat%control(i)%ix_slave
    ds = lat%ele(ix)%s - ele%s
    if (ds > 0) ds = ds - lat%param%total_length
    if (-ds > ele%value(l$)) then
      print *, 'ERROR IN ORDER_SUPER_LORD_SLAVES: INTERNAL ERROR!'
      call err_exit
    endif
    s_rel(i+1-ix1) = ds
  enddo

! Sort slaves by distance.

  call indexx (s_rel, ixx)
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
