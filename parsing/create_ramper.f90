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
!   lord           -- ele_struct: Ramper element.
!     %control%type
!   contrl(:)      -- Control_struct: control info. 1 element for each slave.
!     %stack(:)      -- Arithmetic expression stack for evaluating the controlled parameter value.
!     %y_knot(:)     -- Knot points for spline or linear interpolation.
!     %attribute     -- name of attribute to be controlled
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!
! Output:
!   lord          -- ele_struct: Modified ramper elment
!-

subroutine create_ramper (lord, contrl, err)

use bmad_parser_mod, except_dummy => create_ramper
use expression_mod

implicit none

type (ele_struct), target :: lord
type (lat_struct), pointer :: lat
type (control_struct)  contrl(:)
type (control_struct), pointer :: c

integer i, j, is, n, iv
integer n_slave, ix_attrib

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

call check_controller_controls (ramper$, contrl, lord%name, err)
if (err) return

lord%lord_status = ramper_lord$
lord%key = ramper$

if (n_slave == 0) return ! If no slaves then nothing to do.

! Loop over all slaves.

if (allocated(lord%control%ramp)) deallocate (lord%control%ramp)
allocate(lord%control%ramp(n_slave))

do j = 1, n_slave
  ix_attrib = contrl(j)%ix_attrib
  attrib_name = contrl(j)%attribute

  !

  c => lord%control%ramp(j)
  c%ix_attrib = contrl(j)%ix_attrib
  c%attribute = contrl(j)%attribute
  c%slave_name = contrl(j)%slave_name

  call parser_transfer_control_struct(contrl(j), c, lord, 1)
enddo

end subroutine


