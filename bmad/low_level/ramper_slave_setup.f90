!+
! Subroutine ramper_slave_setup (lat, force_setup)
!
! Routine to setup slave%controller%ix_ramper_lord(:) array.
!
! Input:
!   lat         -- lat_struct: Lattice to be setup.
!   force_setup -- logical, optional: Default False. 
!                   If True, do the setup even if lat%ramper_slave_bookkeeping = ok$.
!                   But the setup will never be done if lat%ramper_slave_bookkeeping = super_ok$.
!
! Output:
!   lat         -- lat_struct: Lattice with ramper slaves setup.
!-

subroutine ramper_slave_setup (lat, force_setup)

use bmad_routine_interface, dummy => ramper_slave_setup
use attribute_mod, only: attribute_free

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: lord
type (control_ramp1_struct), pointer :: r1
type (ele_pointer_struct), allocatable :: eles(:)

integer ib, ie, ir, iv, n_loc, n_slave
logical, optional :: force_setup
logical err
character(*), parameter :: r_name = 'ramper_slave_setup'

! Clean 

err = .false.
if (lat%ramper_slave_bookkeeping == ok$ .and. .not. logic_option(.false., force_setup)) return
if (lat%ramper_slave_bookkeeping == super_ok$) return 

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch%ele%n_lord_ramper = 0
enddo

! 

do ir = lat%n_ele_track+1, lat%n_ele_max
  lord => lat%ele(ir)
  if (lord%key /= ramper$) cycle
  n_slave = 0
  do iv = 1, size(lord%control%ramp)
    r1 => lord%control%ramp(iv)
    call lat_ele_locator(r1%slave_name, lat, eles, n_loc, err)
    if (err) return
    n_slave = n_slave + n_loc
    do ie = 1, n_loc
      call set_this_slave(eles(ie)%ele, lord, iv, r1, lat, n_slave, err)
      if (err) return
    enddo
  enddo

  if (n_slave == 0) then
    call out_io (s_warn$, r_name, 'Ramper: ' // lord%name, &
                                  'Does not control any lattice element attributes!')
  endif
enddo

lat%ramper_slave_bookkeeping = ok$

!----------------------------------------------------------------------------------
contains

recursive subroutine set_this_slave (slave, lord, ix_control, r1, lat, n_slave, err_flag)

type (ele_struct), target :: slave, lord
type (ele_struct), pointer :: slave2
type (lat_struct) lat
type (controller_struct), pointer :: ctl
type (control_ramp1_struct), pointer :: r1
type (ramper_lord_struct), pointer :: r0
type (all_pointer_struct) a_ptr

integer ix_control, is, n_slave
logical err_flag, err, has_wild, free, energy_ramp

! If the slave name has wild card characters, do not match to controllers.

err_flag = .false.
has_wild = index(r1%slave_name, '*') /= 0 .or. index(r1%slave_name, '%') /= 0
energy_ramp = (r1%attribute == 'P0C' .or. r1%attribute == 'E_TOT')
r1%is_controller = (slave%key == overlay$ .or. slave%key == group$ .or. slave%key == girder$)
if (has_wild .and. r1%is_controller) return

if (slave%key == ramper$ .and. .not. has_wild) then
  call out_io (s_error$, r_name, 'RAMPER: ' // lord%name, &
                                 'MAY NOT CONTROL ANOTHER RAMPER: ' // r1%slave_name)
  err_flag = .true.
  return
endif

! Check attribute.

call pointer_to_attribute (slave, r1%attribute, .true., a_ptr, err, .false.)
if (err .or. .not. associated(a_ptr%r)) then
  if (has_wild) return
  call out_io (s_error$, r_name, 'BAD SLAVE ATTRIBUTE FOR RAMPER LORD: ' // lord%name, &
                                 'ATTRIBUTE: ' // r1%attribute, &
                                 'CONTROLLING SLAVE: ' // slave%name)
  err_flag = .true.
  return
endif

! Do not set super_slave elements. EG: If attribute is hkick it must be divided between the slaves.
! Exception: If p0c or E_tot is being set

if (slave%slave_status == super_slave$ .and. .not. energy_ramp) then
  if (.not. has_wild) then
    call out_io (s_error$, r_name, 'RAMPER: ' // lord%name, &
                                   'MAY NOT CONTROL A SUPER SLAVE ELEMENT: ' // r1%slave_name)
    err_flag = .true.
  endif
  return
endif

! Is free


free = attribute_free(slave, r1%attribute, .false.)

if (.not. free .and. .not. energy_ramp) then
  if (.not. has_wild) then
    call out_io (s_error$, r_name, 'RAMPER: ' // lord%name, &
                                   'IS TRYING TO CONTROL A NON-FREE ATTRIBUTE: ' // r1%attribute, &
                                   'OF SLAVE ELEMENT: ' // r1%slave_name)
    err_flag = .true.
  endif
  return
endif

!

if (.not. associated(slave%control)) allocate(slave%control)
ctl => slave%control

if (slave%n_lord_ramper == 0) then
  ctl%ramper_lord = [ramper_lord_struct(lord%ix_ele, ix_control, a_ptr%r)]
  slave%n_lord_ramper = slave%n_lord_ramper + 1
  n_slave = n_slave + 1
else
  ! Avoid duplicate entries. 
  ! EG: "*" element match will match to all elements and super slaves will be duplicated.
  r0 => ctl%ramper_lord(slave%n_lord_ramper)
  if (r0%ix_ele /= lord%ix_ele .or. r0%ix_con /= ix_control) then
    ctl%ramper_lord = [ctl%ramper_lord, ramper_lord_struct(lord%ix_ele, ix_control, a_ptr%r)]
    slave%n_lord_ramper = slave%n_lord_ramper + 1
    n_slave = n_slave + 1
  endif
endif

end subroutine set_this_slave

end subroutine
