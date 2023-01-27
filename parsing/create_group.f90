!+
! Subroutine create_group (lord, contrl, err)
!
! Subroutine to add the controller information to slave elements of a group_lord.
!
! Note: If the stack (in contrl(i)%stack(:)) array has a single numeric term, and
! if there is only one control variable, then
! the arithmatic expression is modified so that the controlled attribute is linear
! in lord%control%var(1) with a coefficient given by the single numeric term.
!
! Note: See the Bmad manual for directions as to how to use this routine.
!
! Input:
!   lord           -- ele_struct: Group element.
!     %control%type
!   contrl(:)      -- Control_struct: control info. 1 element for each slave.
!     %stack         -- Arithmetic expression stack for evaluating the controlled parameter value.
!     %y_knot(:)     -- Knot points for spline interpolation.
!     %slave         -- Integer: Index to lat%branch()%ele() of element controlled.
!     %attribute     -- character(40): Attribute name.
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!
! Output:
!   lord          -- ele_struct: Modified group elment
!-

subroutine create_group (lord, contrl, err)

use bmad_parser_mod, except_dummy => create_group
use expression_mod

implicit none

type (ele_struct), target :: lord
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: slave
type (control_struct)  contrl(:)
type (branch_struct), pointer :: branch
type (control_struct), pointer :: c

integer i, j, n_control, n_con, is, iv, n, ix_attrib
integer ix1, ix2, ix_min, ix_max, ix_slave, ix_branch

logical err, err2, free

character(40) attrib_name

! Error check

n_control = size(contrl)
lat => lord%branch%lat

do i = 1, n_control
  ix_slave  = contrl(i)%slave%ix_ele
  ix_branch = contrl(i)%slave%ix_branch

  if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) then
    call parser_error ( 'BRANCH INDEX OUT OF BOUNDS. \i0\ ', &
                                    'CONSTRUCTING GROUP: ' // lord%name, i_array = [ix_branch])
    return
  endif

  if (ix_slave < 0 .or. ix_slave > ubound(lat%branch(ix_branch)%ele, 1)) then
    call parser_error ('LATTICE ELEMENT INDEX OUT OF BOUNDS. \i0\ ', &
                       'CONSTRUCTING GROUP: ' // lord%name, i_array = [ix_slave])
    return
  endif
enddo

! init

call check_controller_controls (group$, contrl, lord%name, err)
if (err) return

lord%lord_status = group_lord$
lord%key = group$

if (n_control == 0) return ! If no slaves then nothing to do.

err = .true.
n_con = lat%n_control_max
lord%ix1_slave = n_con + 1

do iv = 1, size(lord%control%var)
  call upcase_string(lord%control%var(iv)%name)
enddo

! loop over all controlled attributes

do i = 1, n_control

  ! For position control: We need to figure out the elements that
  ! need to be controlled.
  ! Find beginning and ending positions of element
  ! if a super_lord then we must go to the slave elements to find the ends
  ! else not a super lord so finding the ends is simple

  ix_slave = contrl(i)%slave%ix_ele
  ix_branch = contrl(i)%slave%ix_branch
  attrib_name = contrl(i)%attribute
  branch => lat%branch(ix_branch)

  slave => branch%ele(ix_slave)
  if (slave%slave_status == free$) slave%slave_status = minor_slave$

  ! If the slave attribute is a multipole component, make sure it exists.

  if (is_attribute(ix_attrib, multipole$) .and. .not. associated (slave%a_pole)) then
    call multipole_init(slave, magnetic$)
  endif

  if (is_attribute(ix_attrib, elec_multipole$) .and. .not. associated (slave%a_pole_elec)) then
    call multipole_init(slave, electric$)
  endif

  ! Varying the length of a super_slave is permitted so do not check in this case.

  select case (attrib_name)
  case ('START_EDGE', 'END_EDGE', 'ACCORDION_EDGE', 'S_POSITION')
    free = attribute_free (slave, 'L', .false., .false., .true.)
  case default
    free = attribute_free (slave, attrib_name, .false., .false., .true.)
  end select

  if (.not. free) then
    call parser_error ('SLAVE ATTRIBUTE NOT FREE TO VARY FOR GROUP LORD: ' // lord%name)
    err = .true.
    return
  endif

  !

  n_con = n_con + 1
  if (n_con > size(lat%control)) call reallocate_control (lat, n_con+100)
  c => lat%control(n_con)
  c%ix_attrib = contrl(i)%ix_attrib
  c%slave     = contrl(i)%slave
  c%attribute = contrl(i)%attribute
  c%lord      = lat_ele_loc_struct(lord%ix_ele, 0)

  call add_lattice_control_structs (slave, n_add_lord = 1)
  lat%ic(slave%ic1_lord+slave%n_lord-1) = n_con

  call parser_transfer_control_struct(contrl(i), c, lord, 1)

enddo

! End stuff

lord%n_slave = n_con - lord%ix1_slave + 1
lat%n_control_max = n_con

err = .false.

end subroutine
