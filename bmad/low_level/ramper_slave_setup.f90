!+
! Subroutine ramper_slave_setup (lat)
!
! Routine to setup slave%controller%ix_ramper_lord(:) array.
!
! Input:
!   lat     -- lat_struct: Lattice to be setup.
!
! Output:
!   lat     -- lat_struct: Lattice with ramper slaves setup.
!-

subroutine ramper_slave_setup (lat)

use bmad_routine_interface, dummy => ramper_slave_setup

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele, lord
type (control_ramp1_struct), pointer :: r1
type (ele_pointer_struct), allocatable :: eles(:)

integer ib, ie, ir, iv, n_loc
logical err
character(*), parameter :: r_name = 'ramper_slave_setup'

! Clean 

err = .false.
if (lat%ramper_slave_bookkeeping_done) return

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch%ele%n_lord_ramper = 0
enddo

! 

do ir = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ir)
  if (lord%key /= ramper$) cycle
  do iv = 1, size(lord%control%ramp)
    r1 => lord%control%ramp(iv)
    call lat_ele_locator(r1%slave_name, lat, eles, n_loc, err)
    if (err) return
    do ie = 1, n_loc
      call set_this_slave(eles(ie)%ele, lord, iv, lat, err); if (err) return
    enddo
  enddo
enddo

lat%ramper_slave_bookkeeping_done = .true.

!----------------------------------------------------------------------------------
contains

recursive subroutine set_this_slave (slave, lord, ix_control, lat, err)

type (ele_struct), target :: slave, lord
type (ele_struct), pointer :: slave2
type (lat_struct) lat
type (controller_struct), pointer :: ctl
type (ramper_lord_struct), pointer :: r0

integer ix_control, is
logical err

!

if (.not. associated(slave%control)) allocate(slave%control)
ctl => slave%control

if (slave%n_lord_ramper == 0) then
  ctl%ramper_lord = [ramper_lord_struct(lord%ix_ele, ix_control)]
  slave%n_lord_ramper = slave%n_lord_ramper + 1
else
  ! Avoid duplicate entries. 
  ! EG: "*" element match will match to all elements and super slaves will be duplicated.
  r0 => ctl%ramper_lord(slave%n_lord_ramper)
  if (r0%ix_ele /= lord%ix_ele .or. r0%ix_con /= ix_control) then
    ctl%ramper_lord = [ctl%ramper_lord, ramper_lord_struct(lord%ix_ele, ix_control)]
    slave%n_lord_ramper = slave%n_lord_ramper + 1
  endif
endif

! And mark slaves of slave as well

do is = 1, slave%n_slave
  slave2 => pointer_to_slave(slave, is)
  call set_this_slave (slave2, lord, ix_control, lat, err)
enddo

end subroutine set_this_slave

end subroutine
