!+
! Subroutine tao_evaluate_stack_new (stack, n_size_in, use_good_user, value, info, err_flag, print_err, expression)
!
! Routine to evaluate an expression stack.
!
! Input:
!   stack(:)      -- tao_eval_stack1_struct: Expression stack
!   n_size_in     -- integer: Desired array size. If the expression evaluates to a
!                      a scalar, each value in the value array will get this value.
!                      If n_size = 0, the natural size is determined by the expression itself.
!   use_good_user -- logical: Use the good_user logical in evaluating good(:)
!   print_err     -- logical: If False then supress evaluation error messages.
!                      This does not affect syntax error messages. Default is True.
!   expression    -- character(*): Original expression. Used for error messages.
!
! Output:
!   value(:)      -- Real(rp), allocatable: Value of arithmetic expression.
!   info(:)       -- tao_expression_info_struct, allocatable: Is the value valid? 
!                      Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                      orbit.x[23]|good_user is False.
!   err_flag      -- Logical: True on error. False otherwise
!-

subroutine tao_evaluate_stack_new (stack, n_size_in, use_good_user, value, err_flag, print_err, expression, info_in)

use tao_interface, dummy => tao_evaluate_stack_new
use expression_mod

type (tao_eval_stack1_struct), target :: stack(:)
type (tao_eval_stack1_struct), pointer :: ss
type (tao_eval_stack1_struct), pointer :: s2(:)
type (tao_eval_stack1_struct) stk2(20)
type (tao_expression_info_struct), allocatable, optional :: info_in(:)
type (tao_expression_info_struct), allocatable :: info(:)

real(rp), allocatable :: value(:)

integer n_size_in, species
integer i, i2, j, n, ns, ni, n_size

logical err_flag, use_good_user, print_err, info_allocated

character(*) expression
character(*), parameter :: r_name = 'tao_evaluate_stack_new'

! Calculate good

s2 => stack   ! For debugging purposes
err_flag = .true.
n_size = max(1, n_size_in)

do i = 1, size(stack)
  ss => stack(i)

  select case (ss%name)
  case ('average', 'sum', 'rms', 'min', 'max'); n_size = 1
  end select

  if (allocated(ss%value)) then
    if (size(ss%value) > 1 .and. n_size > 1 .and. size(ss%value) /= n_size) then
      if (print_err) call out_io (s_error$, r_name, 'Array sizes mismatch in expression')
      err_flag = .true.
      return
    endif
    n_size = max(n_size, size(ss%value))
  endif

  if (allocated(ss%value_ptr)) then
    if (associated(ss%value_ptr(1)%good_value)) then    
      do j = 1, size(ss%value_ptr)
        if (use_good_user) then
          ss%info(j)%good = ss%value_ptr(j)%good_value .and. ss%value_ptr(j)%good_user
        else
          ss%info(j)%good = ss%value_ptr(j)%good_value
        endif
      enddo
    endif
  endif
enddo

call tao_re_allocate_expression_info(info, n_size)

! Go through the stack and perform the operations...

i2 = 0  ! stack pointer
do i = 1, size(stack)

  if (allocated(stack(i)%info)) then
    ns = size(stack(i)%info)
    if (.not. allocated(info)) then
      call tao_re_allocate_expression_info(info, ns)
      info = tao_expression_info_struct()
    endif
    ni = size(info)

    if (.not. this_size_check(ns, ni)) return
    if (ns == ni) then
      info%good = info%good .and. stack(i)%info%good
      do j = 1, size(info)
        if (stack(i)%info(j)%s /= real_garbage$) info(j)%s = stack(i)%info(j)%s
        info(j)%ele => stack(i)%info(j)%ele
      enddo 
    elseif (ns == 1) then
      info%good = info%good .and. stack(i)%info(1)%good
    elseif (ni == 1) then
      call tao_re_allocate_expression_info(info, ns)
      info%good = info(1)%good .and. stack(i)%info(1)%good
    endif
  endif

  !

  select case (stack(i)%type)
  case (arg_count$)
    cycle

  case (numeric$)
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (species_const$) ! Something like "electron". Just push on stack.
    i2 = i2 + 1
    stk2(i2)%name = stack(i)%name
    call re_allocate(stk2(i2)%value, 1)

  case (species$)
    stk2(i2)%value = species_id(stk2(i2)%name)

  case (lat_num$, ele_num$)
    !!! This needs to be fixed to include default stuff
    !!! call tao_param_value_routine (stack(i)%name, '', stack(i), err_flag, print_err)
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (var_num$, data_num$)
    do j = 1, size(stack(i)%value)
      stack(i)%value(j) = stack(i)%value_ptr(j)%r
    enddo
    stack(i)%value = stack(i)%value * stack(i)%scale
    i2 = i2 + 1
    call value_transfer (stk2(i2)%value, stack(i)%value)

  case (unary_minus$) 
    stk2(i2)%value = -stk2(i2)%value

  case (unary_plus$) 
    ! Nothing to do

  case (plus$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) + stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value + stk2(i2)%value
    endif
    i2 = i2 - 1

  case (minus$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) - stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value - stk2(i2)%value
    endif
    i2 = i2 - 1

  case (times$)
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) * stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value * stk2(i2)%value
    endif
    i2 = i2 - 1

  case (divide$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    n = size(stk2(i2)%value)
    do j = 1, n
      if (stk2(i2)%value(j) /= 0) cycle
      if (print_err) call out_io (s_error$, r_name, 'Divide by zero!')
      err_flag = .true.
      return
      ! Propably can get rid of this stuff...
      stk2(i2)%value(j) = 1
      if (n == 1) then
        info%good = .false.  ! All are false
      else
        info(j)%good = .false.
      endif
    enddo

    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) / stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value / stk2(i2)%value
    endif
    i2 = i2 - 1

  case (power$)
    if (.not. this_size_check(size(stk2(i2)%value), size(stk2(i2-1)%value))) return
    if (size(stk2(i2)%value) < size(stk2(i2-1)%value)) then
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value(1)
    elseif (size(stk2(i2)%value) > size(stk2(i2-1)%value)) then
      call value_transfer (stk2(i2-1)%value, stk2(i2-1)%value(1) ** stk2(i2)%value)
    else
      stk2(i2-1)%value = stk2(i2-1)%value ** stk2(i2)%value
    endif
    i2 = i2 - 1

  case (function$)
    select case (stack(i)%name)
    case ('cot');       stk2(i2)%value = 1.0_rp / tan(stk2(i2)%value)
    case ('csc');       stk2(i2)%value = 1.0_rp / sin(stk2(i2)%value)
    case ('sec');       stk2(i2)%value = 1.0_rp / cos(stk2(i2)%value)
    case ('sin');       stk2(i2)%value = sin(stk2(i2)%value)
    case ('sinc');      stk2(i2)%value = sinc(stk2(i2)%value)
    case ('cos');       stk2(i2)%value = cos(stk2(i2)%value)
    case ('tan');       stk2(i2)%value = tan(stk2(i2)%value)
    case ('asin');      stk2(i2)%value = asin(stk2(i2)%value)
    case ('acos');      stk2(i2)%value = acos(stk2(i2)%value)
    case ('atan');      stk2(i2)%value = atan(stk2(i2)%value)
    case ('atan2');     stk2(i2-1)%value = atan2(stk2(i2-1)%value, stk2(i2)%value)
      i2 = i2 - 1
    case ('modulo');    stk2(i2-1)%value = modulo(stk2(i2-1)%value, stk2(i2)%value)
      i2 = i2 - 1
    case ('sinh');      stk2(i2)%value = sinh(stk2(i2)%value)
    case ('cosh');      stk2(i2)%value = cosh(stk2(i2)%value)
    case ('tanh');      stk2(i2)%value = tanh(stk2(i2)%value)
    case ('coth');      stk2(i2)%value = 1.0_rp / tanh(stk2(i2)%value)
    case ('asinh');     stk2(i2)%value = asinh(stk2(i2)%value)
    case ('acosh');     stk2(i2)%value = acosh(stk2(i2)%value)
    case ('atanh');     stk2(i2)%value = atanh(stk2(i2)%value)
    case ('acoth');     stk2(i2)%value = 1.0_rp / atanh(stk2(i2)%value)
    case ('abs');       stk2(i2)%value = abs(stk2(i2)%value)
    case ('sqrt');      stk2(i2)%value = sqrt(stk2(i2)%value)
    case ('log');       stk2(i2)%value = log(stk2(i2)%value)
    case ('exp');       stk2(i2)%value = exp(stk2(i2)%value)
    case ('int');       stk2(i2)%value = int(stk2(i2)%value)
    case ('sign');      stk2(i2)%value = sign_of(stk2(i2)%value)
    case ('nint');      stk2(i2)%value = nint(stk2(i2)%value)
    case ('floor');     stk2(i2)%value = floor(stk2(i2)%value)
    case ('ceiling');   stk2(i2)%value = ceiling(stk2(i2)%value)

    case ('rms');       stk2(i2)%value(1) = rms_value(stk2(i2)%value, info%good)
      call re_allocate(stk2(i2)%value, 1)
      info(1)%good = any(info%good)
      call tao_re_allocate_expression_info(info, 1)

    case ('average', 'mean')
      if (any(info%good)) then
        stk2(i2)%value(1) = sum(stk2(i2)%value, mask = info%good) / count(info%good)
      endif
      call re_allocate(stk2(i2)%value, 1)
      info(1)%good = any(info%good)
      call tao_re_allocate_expression_info(info, 1)

    case ('sum', 'min', 'max')
      select case (stack(i)%name)
      case ('sum'); stk2(i2)%value(1) = sum(stk2(i2)%value, mask = info%good)
      case ('min'); stk2(i2)%value(1) = minval(stk2(i2)%value, mask = info%good)
      case ('max'); stk2(i2)%value(1) = maxval(stk2(i2)%value, mask = info%good)
      end select
      call re_allocate(stk2(i2)%value, 1)
      info(1)%good = .true.
      call tao_re_allocate_expression_info(info, 1)

    case ('factorial');
    do n = 1, size(stk2(i2)%value)
      stk2(i2)%value(n) = factorial(nint(stk2(i2)%value(n)))
    enddo

    case ('ran');
    i2 = i2 + 1
    call re_allocate(stk2(i2)%value, n_size)
    call ran_uniform(stk2(i2)%value)

    if (size(info) == 1) then
      call tao_re_allocate_expression_info(info, n_size)
      info%good = info(1)%good
    endif

    case ('ran_gauss')
      if (nint(stack(i-1)%value(1)) == 0) then
        i2 = i2 + 1
        call re_allocate(stk2(i2)%value, n_size)
        call ran_gauss(stk2(i2)%value)
      else
        call re_allocate(value, n_size)
        call ran_gauss(value, sigma_cut = stk2(i2)%value(1))
        call re_allocate(stk2(i2)%value, n_size)
        stk2(i2)%value = value
      endif

      if (size(info) == 1) then
        call tao_re_allocate_expression_info(info, n_size)
        info%good = info(1)%good
      endif

    case ('mass_of', 'charge_of', 'anomalous_moment_of')
      species = species_id(stk2(i2)%name)
      if (species == invalid$) then
        if (print_err) call out_io (s_error$, r_name, 'Not a valid species name: ' // stk2(i2)%name)
        err_flag = .true.
        return
      endif
      select case (stack(i)%name)
      case ('mass_of');              stk2(i2)%value = mass_of(species)
      case ('charge_of');            stk2(i2)%value = charge_of(species)
      case ('anomalous_moment_of');  stk2(i2)%value = anomalous_moment_of(species)
      end select
    end select

  case default
    call out_io (s_error$, r_name, 'INTERNAL ERROR')
    call err_exit
  end select
enddo

!

if (i2 /= 1) then
  call out_io (s_error$, r_name, 'INTERNAL ERROR')
  call err_exit
endif

if (size(stk2(1)%value) == 1 .and. n_size_in > 1) then
  call re_allocate(value, n_size_in)
  value = stk2(1)%value(1)
  if (.not. info(1)%good) value = 0
elseif (size(stk2(1)%value) > 1 .and. size(info) == 1) then
  call value_transfer (value, stk2(1)%value)
  if (.not. info(1)%good) value = 0
else
  call value_transfer (value, stk2(1)%value)
  where (.not. info%good) value = 0
endif

if (present(info_in)) then
  if (allocated(info_in)) deallocate(info_in)
  info_in = info
endif

n_size = size(value)
if (n_size_in /= 0) then
  if (n_size /= 1 .and. n_size_in /= n_size) then
    call out_io (s_error$, r_name, 'ARRAY SIZE MISMATCH FROM WHAT IS EXPECTED IN EXPRESSION: ' // expression)
    return
  endif
endif

err_flag = .false.

!-------------------------------------------------------------------------
contains

subroutine value_transfer (to_array, from_array)

real(rp), allocatable :: to_array(:)
real(rp) from_array(:)

!

call re_allocate (to_array, size(from_array))
to_array = from_array

end subroutine value_transfer

!-------------------------------------------------------------------------
! contains

function this_size_check (isize1, isize2) result (ok)
integer isize1, isize2
logical ok
ok = (isize1 == 1 .or. isize2 == 1 .or. isize1 == isize2)
if (.not. ok) then
  call out_io (s_error$, r_name, 'ARRAY SIZE MISMATCH IN EXPRESSION: ' // expression)
endif
end function this_size_check

end subroutine tao_evaluate_stack_new 

