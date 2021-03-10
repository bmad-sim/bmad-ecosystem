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

logical err, err2, free, var_found

! Error check

lat => lord%branch%lat
n_slave = size (contrl)
err = .true.

do iv = 1, size(lord%control%var)
  call upcase_string(lord%control%var(iv)%name)
enddo

if (lord%control%type /= expression$ .and. size(lord%control%var) /= 1) then
  call parser_error ('A SPLINE BASED CONTROLLER MAY ONLY HAVE ONE CONTROL VARIABLE: ' // lord%name)
  return
endif  

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

  !

  if (lord%control%type == expression$) then

    do is = 1, size(contrl(j)%stack)
      if (contrl(j)%stack(is)%type == end_stack$) exit
    enddo
    call reallocate_expression_stack(c%stack, is-1)

    c%stack = contrl(j)%stack(1:is-1)

    ! Convert variable$ type to ramper variable index if name matches an ramper variable name.

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
      if (.not. is_attribute (c%stack(is)%type, all_control_var$)) cycle
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

    ! Evaluate any variable values.

    do is = 1, size(c%stack)
      select case (c%stack(is)%type)
      case (ran$, ran_gauss$)
        call parser_error ('RANDOM NUMBER FUNCITON MAY NOT BE USED WITH AN RAMPER OR GROUP', &
                           'FOR ELEMENT: ' // lord%name)
        return
      case (variable$)
        call word_to_value (c%stack(is)%name, lat, c%stack(is)%value, err2)
        err = (err .or. err2)
        if (err) then
          call parser_error ('ERROR CONVERTING WORD TO VALUE: ' // c%stack(is)%name, &
                             'FOR ELEMENT: ' // lord%name)
          return
        endif
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
    c%y_knot = contrl(j)%y_knot
    if (size(c%y_knot) /= size(lord%control%x_knot)) then
      call parser_error ('NUMBER OF Y_SPLINE POINTS CONTROLLING PARAMETER: ' // c%attribute, &
                         'IS NOT THE SAME AS THE NUMBER OF X_SPLINE POINTS FOR ELEMENT: ' // lord%name)
      return
    endif
  endif

enddo

end subroutine


