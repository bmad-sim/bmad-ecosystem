!+
! Subroutine create_feedback (lord, pickup, kicker, err_flag)
!
! Subroutine to add the lord/slave bookkeeping information for a feedback lord element.
!
! Input:
!   lord          -- ele_struct: Feedback element.
!   pickup        -- character(*): Name of pickup slave.
!   kicker        -- character(*): Name of kicker slave.
!   err_flag      -- Logical: Set True if there is a problem.
!
! Output:
!   lord          -- ele_struct: Modified feedback elment.
!-

subroutine create_feedback (lord, pickup, kicker, err_flag)

use bmad_parser_mod, except_dummy => create_feedback

implicit none

type (ele_struct), target :: lord
type (ele_struct), pointer :: slave
type (lat_struct), pointer :: lat

integer nc0
character(*) pickup, kicker

logical err_flag, err

! Error check

lat => lord%branch%lat
err_flag = .true.

nc0 = lat%n_control_max
lat%n_control_max = nc0 + 2
call reallocate_control(lat, lat%n_control_max)

lord%lord_status = control_lord$
lord%n_slave = 2
lord%ix1_slave = nc0 + 1

call add_this_slave(lat, lord, pickup, 'PICKUP', nc0, 1, err); if (err) return
call add_this_slave(lat, lord, kicker, 'KICKER', nc0, 2, err); if (err) return

err_flag = .false.

!------------------------------------------------------------------------
contains

subroutine add_this_slave (lat, lord, name, who, nc0, ix_slave, err)
type (lat_struct), target :: lat
type (ele_struct), target :: lord
type (ele_struct), pointer :: slave
type (ele_pointer_struct), allocatable :: eles(:)
type (control_struct), pointer :: c

integer nc0, ix_slave, n_loc
logical err
character(*) name, who

! Locate element

call lat_ele_locator (name, lat, eles, n_loc, err)
if (n_loc == 0 .or. err) then
  call parser_error ('CANNOT FIND ' // name // ' WHICH IS A ' // who // ' OF: ' // lord%name)
  err = .true.
  return
endif

if (n_loc > 1) then
  call parser_error ('MULTIPLE ELEMENTS MATCH NAME: ' // name // 'WHICH IS A ' // who // ' OF: ' // lord%name)
  err = .true.
  return
endif

slave => eles(1)%ele

!

c => lat%control(nc0+ix_slave)
c%lord  = lat_ele_loc_struct(lord%ix_ele, 0)
c%slave = lat_ele_loc_struct(slave%ix_ele, slave%ix_branch)
c%slave_name = slave%name

if (slave%slave_status == free$) slave%slave_status = minor_slave$

! You cannot control super_slaves 

if (slave%slave_status == super_slave$ .or. slave%slave_status == multipass_slave$) then
  call parser_error ('FEEDBACK ELEMENT ' // trim(lord%name) // ' CANNOT CONTROL SUPER_SLAVE NOR MULTIPASS_SLAVE: ' // slave%name)
  return
endif

call add_lattice_control_structs(slave, n_add_lord = 1)
lat%ic(slave%ic1_lord+slave%n_lord-1) = lord%ix1_slave + ix_slave - 1

end subroutine add_this_slave

end subroutine


