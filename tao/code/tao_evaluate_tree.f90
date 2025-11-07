!+
! Subroutine tao_evaluate_tree (tao_tree, n_size_in, use_good_user, value, info, err_flag, print_err, expression)
!
! Routine to evaluate an expression tao_tree.
!
! Input:
!   tao_tree      -- tao_eval_node_struct: Expression tree
!   n_size_in     -- integer: Desired array size. If the expression evaluates to a
!                      a scalar, each value in the value array will get this value.
!                      If n_size = 0, the natural size is determined by the expression itself.
!   use_good_user -- logical: Use the good_user logical in evaluating good(:)
!   print_err     -- logical: If False then supress evaluation error messages.
!                      This does not affect syntax error messages. Default is True.
!   expression    -- character(*): Original expression. Used for error messages.
!
! Output:
!   value(:)      -- Real(rp), allocatable: Value(s) of the arithmetic expression.
!   info(:)       -- tao_expression_info_struct, allocatable: Are the returned values valid? 
!                      Example: 'orbit.x[23]|meas' is not good if orbit.x[23]|good_meas or
!                      orbit.x[23]|good_user is False.
!   err_flag      -- Logical: True on error. False otherwise
!-

recursive &
subroutine tao_evaluate_tree (tao_tree, n_size_in, use_good_user, value, err_flag, print_err, expression, info_in)

use tao_interface, dummy => tao_evaluate_tree
use expression_mod

type (tao_eval_node_struct), target :: tao_tree
type (tao_eval_node_struct), pointer :: node1
type (tao_eval_node_struct) stk2(20)
type (tao_expression_info_struct), allocatable, optional :: info_in(:)
type (tao_expression_info_struct), allocatable :: info_loc(:), info2(:)

real(rp) val
real(rp), allocatable :: value(:), val2(:)

integer n_size_in, n_size, id_species
integer i, i2, j, n, nn, nj, ns, ni

logical err_flag, use_good_user, print_err, info_allocated, err

character(*) expression
character(*), parameter :: r_name = 'tao_evaluate_tree'

! Handle root

if (tao_tree%type == root$) then  ! Must have only one child
  call tao_evaluate_tree (tao_tree%node(1), n_size_in, use_good_user, value, err_flag, print_err, expression, info_in)
  return
endif

!

if (.not. associated(tao_tree%node)) return

err_flag = .true.
nn = size(tao_tree%node)

! Calculate good

n_size = max(1, n_size_in)

do i = 1, nn
  node1 => tao_tree%node(i)

  select case (node1%name)
  case ('average', 'sum', 'rms', 'min', 'max'); n_size = 1
  end select

  if (allocated(node1%value)) then
    if (size(node1%value) > 1 .and. n_size > 1 .and. size(node1%value) /= n_size) then
      if (print_err) call out_io (s_error$, r_name, 'Array sizes mismatch in expression')
      err_flag = .true.
      return
    endif
    n_size = max(n_size, size(node1%value))
  endif

  if (allocated(node1%value_ptr)) then
    if (associated(node1%value_ptr(1)%good_value)) then    
      do j = 1, size(node1%value_ptr)
        if (use_good_user) then
          node1%info(j)%good = node1%value_ptr(j)%good_value .and. node1%value_ptr(j)%good_user
        else
          node1%info(j)%good = node1%value_ptr(j)%good_value
        endif
      enddo
    endif
  endif
enddo

call tao_re_allocate_expression_info(info_loc, n_size)

! Go through the tao_tree%node array and perform the operations...

i2 = 0  
do i = 1, nn
  node1 => tao_tree%node(i)

  !

  select case (node1%type)
  case (func_parens$)
    ! Can be something like "max([1,2,3])" or "max(1,2,3)"
    nj = size(node1%node)

    if (nj == 1) then
      call tao_evaluate_tree (node1, 0, use_good_user, val2, err, print_err, expression, info2)
      if (err) return
      nj = size(val2)    
      call re_allocate(node1%value, nj)
      call tao_re_allocate_expression_info(node1%info, nj)
      node1%value = val2
      node1%info = info2

    else
      call re_allocate(node1%value, nj)
      call tao_re_allocate_expression_info(node1%info, nj)

      do j = 1, nj
        call tao_evaluate_tree (node1%node(j), 1, use_good_user, val2, err, print_err, expression, info2)
        if (err) return
        node1%value(j) = val2(1)
        node1%info(j) = info2(1)
      enddo
    endif

  case (square_brackets$)
    nj = size(node1%node)
    call re_allocate(node1%value, nj)
    call tao_re_allocate_expression_info(node1%info, nj)

    do j = 1, nj
      call tao_evaluate_tree (node1%node(j), 1, use_good_user, val2, err, print_err, expression, info2)
      if (err) return
      node1%value(j) = val2(1)
      node1%info(j) = info2(1)
    enddo

  case (compound$, parens$, comma$)
    call tao_evaluate_tree (node1, 0, use_good_user, val2, err, print_err, expression, info2)
    if (err) return
    nj = size(val2)    
    call re_allocate(node1%value, nj)
    call tao_re_allocate_expression_info(node1%info, nj)
    node1%value = val2
    node1%info = info2
  end select

  !

  if (allocated(node1%info) .and. node1%type /= func_parens$) then
    ns = size(node1%info)
    if (.not. allocated(info_loc)) then
      call tao_re_allocate_expression_info(info_loc, ns)
      info_loc = tao_expression_info_struct()
    endif
    ni = size(info_loc)

    if (.not. this_size_check(ns, ni)) return
    if (ns == ni) then
      info_loc%good = info_loc%good .and. node1%info%good
      do j = 1, size(info_loc)
        if (node1%info(j)%s /= real_garbage$) info_loc(j)%s = node1%info(j)%s
        info_loc(j)%ele => node1%info(j)%ele
      enddo 
    elseif (ns == 1) then
      info_loc%good = info_loc%good .and. node1%info(1)%good
    elseif (ni == 1) then
      call tao_re_allocate_expression_info(info_loc, ns)
      info_loc%good = info_loc(1)%good .and. node1%info(1)%good
    endif
  endif

  !

  select case (node1%type)

  case (numeric$, constant$, square_brackets$, compound$, func_parens$, parens$, comma$, variable$)
    i2 = i2 + 1
    call node_transfer (stk2(i2), node1)

  case (species_const$) ! Something like "electron". Just push on stack.
    i2 = i2 + 1
    call node_transfer (stk2(i2), node1)
    call re_allocate(stk2(i2)%value, 1)

  case (lat_num$, ele_num$)
    !!! This needs to be fixed to include default stuff
    !!! call tao_param_value_routine (node1%name, '', node1, err_flag, print_err)
    i2 = i2 + 1
    call node_transfer (stk2(i2), node1)

  case (var_num$, data_num$)
    do j = 1, size(node1%value)
      node1%value(j) = node1%value_ptr(j)%r
    enddo
    node1%value = node1%value * node1%scale
    i2 = i2 + 1
    call node_transfer (stk2(i2), node1)

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
        info_loc%good = .false.  ! All are false
      else
        info_loc(j)%good = .false.
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
    select case (node1%name)
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
    case ('atan2')
      val = atan2(stk2(i2)%value(1), stk2(i2)%value(2))
      call re_allocate(stk2(i2)%value, 1)
      stk2(i2)%value = val
    case ('modulo')
      val = modulo(stk2(i2)%value(1), stk2(i2)%value(2))
      call re_allocate(stk2(i2)%value, 1)
      stk2(i2)%value = val
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

    case ('rms');       stk2(i2)%value(1) = rms_value(stk2(i2)%value, info_loc%good)
      call re_allocate(stk2(i2)%value, 1)
      info_loc(1)%good = any(info_loc%good)
      call tao_re_allocate_expression_info(info_loc, 1)

    case ('average', 'mean')
      if (any(info_loc%good)) then
        stk2(i2)%value(1) = sum(stk2(i2)%value, mask = info_loc%good) / count(info_loc%good)
      endif
      call re_allocate(stk2(i2)%value, 1)
      info_loc(1)%good = any(info_loc%good)
      call tao_re_allocate_expression_info(info_loc, 1)

    case ('sum', 'min', 'max')
      select case (node1%name)
      case ('sum'); stk2(i2)%value(1) = sum(stk2(i2)%value, mask = info_loc%good)
      case ('min'); stk2(i2)%value(1) = minval(stk2(i2)%value, mask = info_loc%good)
      case ('max'); stk2(i2)%value(1) = maxval(stk2(i2)%value, mask = info_loc%good)
      end select
      call re_allocate(stk2(i2)%value, 1)
      info_loc(1)%good = .true.
      call tao_re_allocate_expression_info(info_loc, 1)

    case ('factorial');
      do n = 1, size(stk2(i2)%value)
        stk2(i2)%value(n) = factorial(nint(stk2(i2)%value(n)))
      enddo

      case ('ran');
      call re_allocate(stk2(i2)%value, n_size)
      call ran_uniform(stk2(i2)%value)

      if (size(info_loc) == 1) then
        call tao_re_allocate_expression_info(info_loc, n_size)
        info_loc%good = info_loc(1)%good
      endif

    case ('ran_gauss')
      select case (size(tao_tree%node(i-1)%node))
      case (0)
        call re_allocate(stk2(i2)%value, n_size)
        call ran_gauss(stk2(i2)%value)
      case (1)
        call re_allocate(stk2(i2)%value, n_size)
        call ran_gauss(stk2(i2)%value, sigma_cut = stk2(i2)%value(1))
      case default
        call out_io (s_error$, r_name, 'ran_gauss has more than one argument!')
        return
      end select

      if (size(info_loc) == 1) then
        call tao_re_allocate_expression_info(info_loc, n_size)
        info_loc%good = info_loc(1)%good
      endif

    case ('mass_of', 'charge_of', 'anomalous_moment_of', 'species')
      id_species = nint(stk2(i2)%value(1))
      select case (node1%name)
      case ('mass_of');              stk2(i2)%value = mass_of(id_species)
      case ('charge_of');            stk2(i2)%value = charge_of(id_species)
      case ('anomalous_moment_of');  stk2(i2)%value = anomalous_moment_of(id_species)
      case ('species');              stk2(i2)%value = id_species
      end select
    end select

  case default
    call out_io (s_error$, r_name, 'INTERNAL ERROR')
    call err_exit
  end select
enddo

!

if (i2 /= 1) then
  call out_io (s_error$, r_name, 'Malformed expression: ' // expression)
  return
endif

if (size(stk2(1)%value) == 1 .and. n_size_in > 1) then
  call re_allocate(value, n_size_in)
  value = stk2(1)%value(1)
  if (.not. info_loc(1)%good) value = 0
elseif (size(stk2(1)%value) > 1 .and. size(info_loc) == 1) then
  call value_transfer (value, stk2(1)%value)
  if (.not. info_loc(1)%good) value = 0
else
  call value_transfer (value, stk2(1)%value)
  where (.not. info_loc%good) value = 0
endif

if (present(info_in)) then
  if (allocated(info_in)) deallocate(info_in)
  info_in = info_loc
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

subroutine node_transfer (to_node, from_node)

type (tao_eval_node_struct) to_node, from_node

!

call re_allocate (to_node%value, size(from_node%value))
to_node%value = from_node%value
to_node%type = from_node%type
to_node%name = from_node%name

end subroutine node_transfer

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

end subroutine tao_evaluate_tree 

