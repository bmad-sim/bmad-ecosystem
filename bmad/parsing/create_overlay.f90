!+
! Subroutine create_overlay (lord, contrl, err)
!
! Subroutine to add the controller information to slave elements of an overlay_lord.
!
! Note: If the stack (in contrl(i)%stack(:)) array has a single numeric term,
! the arithmatic expression is modified so that the controlled attribute is linear
! in lord%control%var(1) with a coefficient given by the single numeric term.
!
! Input:
!   lord           -- ele_struct: Overlay element.
!     %control%type
!   contrl(:)      -- Control_struct: control info. 1 element for each slave.
!     %stack(:)      -- Arithmetic expression stack for evaluating the controlled parameter value.
!     %y_knot(:)     -- Knot points for spline interpolation.
!     slave%ix_ele   -- Index of element to control
!     %ix_branch     -- Index of branch element belongs to.
!     %attribute     -- name of attribute to be controlled
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!
! Output:
!   lord          -- ele_struct: Modified overlay elment
!-

subroutine create_overlay (lord, contrl, err)

use bmad_parser_mod, except_dummy => create_overlay
use expression_mod

implicit none

type (ele_struct), target :: lord
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: slave
type (control_struct)  contrl(:)
type (control_struct), pointer :: c

integer i, j, nc0, nc2, is, n, iv
integer ix_slave, n_slave, ix_attrib, ix_branch

character(40) attrib_name

logical err, err2, free

! Error check

lat => lord%branch%lat
n_slave = size (contrl)
err = .true.

do iv = 1, size(lord%control%var)
  call upcase_string(lord%control%var(iv)%name)
enddo

do j = 1, n_slave
  ix_slave  = contrl(j)%slave%ix_ele
  ix_branch = contrl(j)%slave%ix_branch

  if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) then
    call parser_error ('BRANCH INDEX OUT OF BOUNDS. \i0\ ', &
                       'CONSTRUCTING OVERLAY: ' // lord%name, i_array = [ix_branch])
    return
  endif

  if (ix_slave < 0 .or. ix_slave > ubound(lat%branch(ix_branch)%ele, 1)) then
    call parser_error ('LATTICE ELEMENT INDEX OUT OF BOUNDS. \i0\ ', &
                       'CONSTRUCTING OVERLAY: ' // lord%name, i_array = [ix_slave])
    return
  endif
enddo

! Mark element as an overlay lord

call check_controller_controls (overlay$, contrl, lord%name, err)
if (err) return

lord%lord_status = overlay_lord$
lord%key = overlay$

if (n_slave == 0) return ! If no slaves then nothing to do.

! Loop over all slaves.

nc0 = lat%n_control_max
nc2 = nc0

do j = 1, n_slave

  ix_attrib = contrl(j)%ix_attrib
  ix_slave = contrl(j)%slave%ix_ele
  attrib_name = contrl(j)%attribute
  slave => lat%branch(contrl(j)%slave%ix_branch)%ele(ix_slave)

  if (nc2+4 > size(lat%control)) call reallocate_control (lat, nc2+100)

  ! If the slave attribute is a multipole component, make sure it exists.

  if (is_attribute(ix_attrib, multipole$) .and. .not. associated (slave%a_pole)) then
    call multipole_init(slave, magnetic$)
  endif

  if (is_attribute(ix_attrib, elec_multipole$) .and. .not. associated (slave%a_pole_elec)) then
    call multipole_init(slave, electric$)
  endif

  !

  free = attribute_free (slave, attrib_name, .true., .true., .true.)
  err = (err .or. .not. free)
  if (err) then
    call parser_error ('ATTRIBUTE NOT FREE TO VARY: ' // attrib_name, &
                       'FOR ELEMENT: ' // lord%name)

    return
  endif
  c => lat%control(nc2+1)

  c%ix_attrib = contrl(j)%ix_attrib
  c%slave     = contrl(j)%slave
  c%attribute = contrl(j)%attribute
  c%lord      = lat_ele_loc_struct(lord%ix_ele, 0)
  nc2 = nc2 + 1

  call parser_transfer_control_struct(contrl(j), c, lord, 1)
enddo

lord%n_slave = n_slave
lord%ix1_slave = nc0 + 1
lat%n_control_max = nc2

! Loop over all slaves
! Free elements convert to minor slaves.

do i = 1, lord%n_slave

  slave => pointer_to_slave(lord, i)
  if (slave%slave_status == free$) slave%slave_status = minor_slave$

  ! You cannot overlay super_slaves 

  if (slave%slave_status == super_slave$) then
    call parser_error ('ILLEGAL OVERLAY ON ' // slave%name, ' BY: ' // lord%name)
    return
  endif

  ! update controller info for the slave ele

  call add_lattice_control_structs (slave, n_add_lord = 1)
  lat%ic(slave%ic1_lord+slave%n_lord-1) = lord%ix1_slave + i - 1
enddo

end subroutine


