!+
! Subroutine create_ramper (lord, contrl, err)
!
! Subroutine to add the controller information to slave elements of an ramper_lord.
!
! Note: If the stack (in contrl(i)%stack(:)) array has a single numeric term,
! the arithmatic expression is modified so that the controlled attribute is linear
! in lord%control%var(1) with a coefficient given by the single numeric term.
!
! Input:
!   lord              -- ele_struct: Ramper element.
!     %control%type
!   contrl(:)         -- Control_struct: control info. 1 element for each slave.
!     %stack(:)         -- Arithmetic expression stack for evaluating the controlled parameter value.
!     %y_knot(:)        -- Knot points for spline or linear interpolation.
!     %attribute        -- name of attribute to be controlled.
!   err               -- Logical: Set True if an attribute is not free to be controlled.
!
! Output:
!   lord              -- ele_struct: Modified ramper elment
!-

subroutine create_ramper (lord, contrl, err)

use bmad_parser_mod, except_dummy => create_ramper
use expression_mod

implicit none

type (ele_struct), target :: lord
type (ele_struct), pointer :: slave
type (lat_struct), pointer :: lat
type (control_struct), target :: contrl(:), con1
type (control_struct), pointer :: con0
type (control_ramp1_struct), pointer :: r1

integer i, j, is, n, iv
integer n_slave, ix_attrib, key

character(40) attrib_name

logical err, err2, free

! Error check

lat => lord%branch%lat
n_slave = size (contrl)
err = .true.

do iv = 1, size(lord%control%var)
  call upcase_string(lord%control%var(iv)%name)
enddo

! Mark element as an ramper lord

lord%lord_status = ramper_lord$
lord%key = ramper$

if (n_slave == 0) return ! If no slaves then nothing to do.

! Loop over all slaves.

if (allocated(lord%control%ramp)) deallocate (lord%control%ramp)
allocate(lord%control%ramp(n_slave))

do j = 1, n_slave
  con0 => contrl(j)
  call parser_transfer_control_struct(con0, con1, lord, 1)

  r1 => lord%control%ramp(j)
  r1%attribute = con0%attribute
  r1%slave_name = con0%slave_name
  if (allocated(con1%y_knot)) r1%y_knot = con1%y_knot
  if (allocated(con1%stack))  r1%stack  = con1%stack

  ! If slave is an overlay, group, or girder, need to mark it as such
  r1%is_controller = .false.
  slave => pointer_to_ele(lat, r1%slave_name)
  if (associated(slave)) then
    key = slave%key
    r1%is_controller = (key == overlay$ .or. key == group$ .or. key == girder$)
    if (r1%is_controller) r1%slave = ele_loc(slave)
  endif
enddo

end subroutine
