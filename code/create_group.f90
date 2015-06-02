!+
! Subroutine create_group (lord, contrl, err, err_print_flag)
!
! Subroutine to add the controller information to slave elements of a group_lord.
!
! Note: If the stack (in contrl(i)%stack(:)) array has a single numeric term,
! the arithmatic expression is modified so that the controlled attribute is linear
! in lord%control_var(1) with a coefficient given by the single numeric term.
!
! Note: See the Bmad manual for directions as to how to use this routine.
!
! Modules Needed:
!   use bmad
!
! Input:
!   lord           -- ele_struct: Overlay element.
!   contrl(:)      -- Control_struct: control info. 1 element for each slave.
!     %stack         -- Arithmetic expression stack for evaluating the controlled parameter value.
!     %slave         -- Integer: Index to lat%branch()%ele() of element controlled.
!     %ix_attrib     -- Integer: Index in %value() array of attribute controlled.
!   err            -- Logical: Set True if an attribute is not free to be controlled.
!   err_print_flag -- Logical, optional: If present and False then suppress.
!                       printing of an error message if attribute is not free.  
!
! Output:
!   lord          -- ele_struct: Modified overlay elment
!-

subroutine create_group (lord, contrl, err, err_print_flag)

use bmad_interface, except_dummy => create_group
use expression_mod

implicit none

type (ele_struct), target :: lord
type (lat_struct), pointer :: lat
type (ele_struct), pointer :: slave
type (control_struct)  contrl(:)
type (branch_struct), pointer :: branch
type (control_struct), pointer :: c

integer i, j, ix_attrib, n_control, n_con, is, iv, n
integer ix1, ix2, ix_min, ix_max, ix_slave, ix_branch

logical err, free, var_found
logical, optional :: err_print_flag

character(16) :: r_name = 'create_group'

! Error check

n_control = size(contrl)
lat => lord%branch%lat

do i = 1, n_control
  ix_slave  = contrl(i)%slave%ix_ele
  ix_branch = contrl(i)%slave%ix_branch

  if (ix_branch < 0 .or. ix_branch > ubound(lat%branch, 1)) then
    call out_io (s_fatal$, r_name, 'BRANCH INDEX OUT OF BOUNDS. \i0\ ', ix_branch)
    if (global_com%exit_on_error) call err_exit
  endif

  if (ix_slave <= 0 .or. ix_slave > ubound(lat%branch(ix_branch)%ele, 1)) then
    call out_io (s_fatal$, r_name, 'INDEX OUT OF BOUNDS. \i0\ ', ix_slave)
    if (global_com%exit_on_error) call err_exit
  endif
enddo

! init

call check_controller_controls (contrl, lord%name, err)
if (err) return

lord%lord_status = group_lord$
lord%key = group$
call set_ele_defaults (lord)

if (n_control == 0) return ! If no slaves then nothing to do.

err = .true.
n_con = lat%n_control_max
lord%ix1_slave = n_con + 1

do iv = 1, size(lord%control_var)
  call upcase_string(lord%control_var(iv)%name)
enddo

! loop over all controlled elements

do i = 1, n_control

  ! For position control: We need to figure out the elements that
  ! need to be controlled.
  ! Find beginning and ending positions of element
  ! if a super_lord then we must go to the slave elements to find the ends
  ! else not a super lord so finding the ends is simple

  ix_slave = contrl(i)%slave%ix_ele
  ix_branch = contrl(i)%slave%ix_branch
  ix_attrib = contrl(i)%ix_attrib
  branch => lat%branch(ix_branch)
  slave => branch%ele(ix_slave)

  ! If the slave attribute is a multipole component, make sure it exists.
  if (is_attribute(ix_attrib, multipole$) .and. .not. associated (slave%a_pole)) then
    call multipole_init(slave)
  endif

  ! Varying the length of a super_slave is permitted so do not check in this case.

  select case (ix_attrib)
  case (start_edge$, end_edge$, accordion_edge$, s_position$)
    free = attribute_free (slave, 'L', err_print_flag)
  case default
    free = attribute_free (slave, attribute_name(slave, ix_attrib), err_print_flag)
  end select

  if (.not. free) then
    if (logic_option(.true., err_print_flag)) call out_io (s_error$, r_name, &
          'SLAVE ATTRIBUTE NOT FREE TO VARY FOR GROUP LORD: ' // lord%name)
    err = .true.
  endif

  !

  n_con = n_con + 1
  if (n_con > size(lat%control)) call reallocate_control (lat, n_con+100)
  c => lat%control(n_con)

  do is = 1, size(contrl(i)%stack)
    if (contrl(i)%stack(is)%type == end_stack$) exit
  enddo
  call reallocate_expression_stack(c%stack, is-1)

  c%stack     = contrl(i)%stack(1:is-1)
  c%ix_attrib = contrl(i)%ix_attrib
  c%slave     = contrl(i)%slave
  c%ix_lord   = lord%ix_ele

  ! Convert numeric$ type to variable index if name matches a variable
  do is = 1, size(c%stack)
    if (c%stack(is)%type == end_stack$) exit
    if (c%stack(is)%type /= numeric$) cycle
    do iv = 1, size(lord%control_var)
      if (upcase(c%stack(is)%name) /= lord%control_var(iv)%name) cycle
      c%stack(is)%type = iv + var_offset$
      exit
    enddo
  enddo

  ! Convert a stack of a single constant "const" to "const * control_var(1)"
  var_found = .false.
  do is = 1, size(c%stack)
    if (c%stack(is)%type < old_var_offset$) cycle
    if (c%stack(is)%type == end_stack$) exit
    var_found = .true.
    exit
  enddo

  if (.not. var_found) then
    if (size(c%stack) == 1 .and. c%stack(1)%name == '1' .or. c%stack(1)%name == '1.0') then
      c%stack(1) = expression_atom_struct(lord%control_var(1)%name, 1+var_offset$, 0.0_rp)
    else
      n = size(c%stack)
      call reallocate_expression_stack(c%stack, n+2)
      c%stack(n+1) = expression_atom_struct(lord%control_var(1)%name, 1+var_offset$, 0.0_rp)
      c%stack(n+2) = expression_atom_struct('', times$, 0.0_rp)
    endif
  endif

  ! Update controller info for the slave element

  slave%n_lord = slave%n_lord + 1
  call add_lattice_control_structs (lat, slave)
  lat%ic(slave%ic2_lord) = n_con
  if (slave%slave_status == free$) slave%slave_status = control_slave$

enddo

! End stuff

lord%ix2_slave = n_con
lord%n_slave = n_con - lord%ix1_slave + 1
lat%n_control_max = n_con

err = .false.

end subroutine
