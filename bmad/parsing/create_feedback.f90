!+
! Subroutine create_feedback (lord, input, output, err_flag)
!
! Subroutine to add the lord/slave bookkeeping information for a feedback lord element.
!
! Input:
!   lord          -- ele_struct: Feedback element.
!   input(:)      -- character(*): Names of input slaves.
!   output(:)     -- character(*): Names of output slaves.
!   err_flag      -- Logical: Set True if there is a problem.
!
! Output:
!   lord          -- ele_struct: Modified feedback elment.
!-

subroutine create_feedback (lord, input, output, err_flag)

use bmad_parser_mod, except_dummy => create_feedback

implicit none

type (ele_struct), target :: lord
type (ele_struct), pointer :: slave
type (lat_struct), pointer :: lat

integer nc0, is, ix_slave
logical err_flag, err

character(*) input(:), output(:)
character(*), parameter :: r_name = 'create_feedback'

! Error check

if (size(input) == 0) then
  call out_io(s_error$, r_name, 'NUMBER OF INPUT_ELE SLAVES IS ZERO FOR FEEDBACK ELEMENT: ' // lord%name)
endif

if (size(output) == 0) then
  call out_io(s_error$, r_name, 'NUMBER OF OUTPUT_ELE SLAVES IS ZERO FOR FEEDBACK ELEMENT: ' // lord%name)
endif

lat => lord%branch%lat
err_flag = .true.

nc0 = lat%n_control_max
ix_slave = 0

lord%lord_status = control_lord$
lord%ix1_slave = nc0 + 1

do is = 1, size(input)
  call add_this_slave(lat, lord, input(is), 'INPUT', nc0, ix_slave, err); if (err) return
enddo

do is = 1, size(output)
  call add_this_slave(lat, lord, output(is), 'OUTPUT', nc0, ix_slave, err); if (err) return
enddo

err_flag = .false.

!------------------------------------------------------------------------
contains

subroutine add_this_slave (lat, lord, ele_id, who, nc0, ix_slave, err)
type (lat_struct), target :: lat
type (ele_struct), target :: lord
type (ele_struct), pointer :: slave
type (ele_pointer_struct), allocatable :: eles(:)
type (control_struct), pointer :: c

integer is, nc0, ix_slave, n_loc
logical err
character(*) ele_id, who

! Locate element

call lat_ele_locator (ele_id, lat, eles, n_loc, err)
if (n_loc == 0 .or. err) then
  call parser_error ('CANNOT FIND ' // ele_id // ' WHICH IS A ' // who // ' OF: ' // lord%name)
  err = .true.
  return
endif

!

do is = 1, n_loc
  slave => eles(is)%ele

  ix_slave = ix_slave + 1
  lord%n_slave = ix_slave

  lat%n_control_max = nc0 + ix_slave
  call reallocate_control(lat, lat%n_control_max)

  c => lat%control(nc0+ix_slave)
  c%lord  = lat_ele_loc_struct(lord%ix_ele, 0)
  c%slave = lat_ele_loc_struct(slave%ix_ele, slave%ix_branch)
  c%slave_name = slave%name
  c%attribute = who

  if (slave%slave_status == free$) slave%slave_status = minor_slave$

  ! You cannot control super_slaves 

  if (slave%slave_status == super_slave$ .or. slave%slave_status == multipass_slave$) then
    call parser_error ('FEEDBACK ELEMENT ' // trim(lord%name) // ' CANNOT CONTROL SUPER_SLAVE NOR MULTIPASS_SLAVE: ' // slave%name)
    return
  endif

  call add_lattice_control_structs(slave, n_add_lord = 1)
  lat%ic(slave%ic1_lord+slave%n_lord-1) = lord%ix1_slave + ix_slave - 1
enddo

end subroutine add_this_slave

end subroutine


