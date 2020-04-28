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

logical err, err2, free, var_found

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

  if (ix_slave <= 0 .or. ix_slave > ubound(lat%branch(ix_branch)%ele, 1)) then
    call parser_error ('LATTICE ELEMENT INDEX OUT OF BOUNDS. \i0\ ', &
                       'CONSTRUCTING GROUP: ' // lord%name, i_array = [ix_slave])
    return
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

do iv = 1, size(lord%control%var)
  call upcase_string(lord%control%var(iv)%name)
enddo

if (lord%control%type /= expression$ .and. size(lord%control%var) /= 1) then
  call parser_error ('A SPLINE BASED CONTROLLER MAY ONLY HAVE ONE CONTROL VARIABLE: ' // lord%name)
  return
endif  

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


  if (lord%control%type == expression$) then

    do is = 1, size(contrl(i)%stack)
      if (contrl(i)%stack(is)%type == end_stack$) exit
    enddo
    call reallocate_expression_stack(c%stack, is-1)

    c%stack = contrl(i)%stack(1:is-1)

    ! Convert variable$ type to group variable index if name matches a group variable name
    do is = 1, size(c%stack)
      if (c%stack(is)%type == end_stack$) exit
      if (c%stack(is)%type /= variable$) cycle
      do iv = 1, size(lord%control%var)
        if (upcase(c%stack(is)%name) /= lord%control%var(iv)%name) cycle
        c%stack(is)%type = iv + var_offset$
        exit
      enddo
    enddo

    ! Convert a stack of a single constant "const" to "const * control_var(1)"
    var_found = .false.
    do is = 1, size(c%stack)
      if (.not. is_attribute(c%stack(is)%type, all_control_var$)) cycle
      if (c%stack(is)%type == end_stack$) exit
      var_found = .true.
      exit
    enddo

    if (.not. var_found) then
      if (size(c%stack) == 1 .and. c%stack(1)%name == '1' .or. c%stack(1)%name == '1.0') then
        c%stack(1) = expression_atom_struct(lord%control%var(1)%name, 1+var_offset$, 0.0_rp)
      else
        n = size(c%stack)
        call reallocate_expression_stack(c%stack, n+2)
        c%stack(n+1) = expression_atom_struct(lord%control%var(1)%name, 1+var_offset$, 0.0_rp)
        c%stack(n+2) = expression_atom_struct('', times$, 0.0_rp)
      endif
    endif

    ! Update controller info for the slave element

    call add_lattice_control_structs (slave, n_add_lord = 1)
    lat%ic(slave%ic1_lord+slave%n_lord-1) = n_con

    ! Evaluate any variable values.

    do is = 1, size(c%stack)
      select case (c%stack(is)%type)
      case (ran$, ran_gauss$)
        call parser_error ('RANDOM NUMBER FUNCITON MAY NOT BE USED WITH A GROUP', &
                           'FOR ELEMENT: ' // lord%name)
        if (global_com%exit_on_error) call err_exit
        return
      case (variable$)
        call word_to_value (c%stack(is)%name, lat, c%stack(is)%value, err); if (err) return
        ! Variables in the arithmetic expression are immediately evaluated and never reevaluated.
        ! If the variable is an element attribute (looks like: "ele_name[attrib_name]") then this may
        ! be confusing if the attribute value changes later. To avoid some (but not all) confusion, 
        ! turn the variable into a numeric$ so the output from the type_ele routine looks "sane".
        if (index(c%stack(is)%name, '[') /= 0) then
          c%stack(is)%type = numeric$
          c%stack(is)%name = ''
        endif
      end select
    enddo

  else
    c%y_knot = contrl(i)%y_knot
    if (size(c%y_knot) /= size(lord%control%x_knot)) then
      call parser_error ('NUMBER OF Y_SPLINE POINTS FOR SLAVE: ' // slave%name, &
                         'IS NOT THE SAME AS THE NUMBER OF X_SPLINE POINTS FOR ELEMENT: ' // lord%name)
      if (global_com%exit_on_error) call err_exit
      return
    endif
  endif

enddo

! End stuff

lord%n_slave = n_con - lord%ix1_slave + 1
lat%n_control_max = n_con

err = .false.

end subroutine
