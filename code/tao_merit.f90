!+
! function tao_merit (calc_ok) result (this_merit)
! 
! function to calculate the merit.
!
! Input:
!
! Output:
!   calc_ok    -- Logical, optional: Set False if there was an error in the 
!                   calculation like a particle was lost or a lat is unstable.
!   this_merit -- Real(rp): Merit value.
!-

function tao_merit (calc_ok) result (this_merit)

use tao_mod
use tao_lattice_calc_mod

implicit none

type (tao_var_struct), pointer :: var
type (tao_data_struct), pointer :: data(:)
type (tao_d1_data_struct), pointer :: d1

real(rp) this_merit, ave, value, model_value

integer i, j, n, iu0

character(16) :: r_name = "tao_merit"

logical, optional :: calc_ok
logical err, ok, opt_with_ref, opt_with_base


! make sure all calculations are up to date.

call tao_lattice_calc (ok)
if (present(calc_ok)) calc_ok = ok

!----------------------------------------
! Merit contribution from the variables.

this_merit = 0

do j = 1, size(s%var)

  var => s%var(j)
  var%delta_merit = 0
  var%merit = 0

  if (.not. var%useit_opt) cycle

  select case (var%merit_type)
  case ('target', 'match')
    if (s%global%opt_with_ref .and. s%global%opt_with_base) then
      if (tao_com%unified_lattices) then
        var%delta_merit = (var%model_value - var%common%model_value) - &
                                                (var%meas_value - var%ref_value)
      else
        var%delta_merit = (var%model_value - var%base_value) - &
                                                (var%meas_value - var%ref_value)
      endif
    elseif (s%global%opt_with_ref) then
      var%delta_merit = (var%model_value - var%design_value) - &
                                                (var%meas_value - var%ref_value)
    elseif (s%global%opt_with_base) then
      if (tao_com%unified_lattices) then
        var%delta_merit = (var%model_value - var%common%model_value) - var%meas_value
      else
        var%delta_merit = (var%model_value - var%base_value) - var%meas_value
      endif
    else
      var%delta_merit = var%model_value - var%meas_value
    endif
  case ('limit')
    if (var%model_value > var%high_lim) then
      var%delta_merit = var%model_value - var%high_lim
    elseif (var%model_value < var%low_lim) then
      var%delta_merit = var%model_value - var%low_lim
    endif
  case default
    call tao_hook_merit_var (i, j, var)
  end select

  var%merit = var%weight * var%delta_merit**2
  this_merit = this_merit + var%merit

enddo

!----------------------------------------
! Merit contribution from the data:

if (tao_com%unified_lattices) iu0 = tao_com%u_common%ix_uni

do i = lbound(s%u, 1), ubound(s%u, 1)

  data => s%u(i)%data
  data%merit = 0
  data%delta_merit = 0

  opt_with_ref = s%global%opt_with_ref
  opt_with_base = s%global%opt_with_base
  if (s%u(i)%common_uni) then
    opt_with_ref = .false.
    opt_with_base = .false.
  endif

! check if universe is turned off

  if (.not. s%u(i)%is_on) cycle

! First compute the delta for the merit function
  
  do j = 1, size(data)
    if (.not. data(j)%useit_opt) cycle
    if (.not. data(j)%good_model) cycle

    if (data(j)%merit_type(1:4) == 'abs_') then
      model_value = data(j)%model_value
    else
      model_value = abs(data(j)%model_value)
    endif

    if (opt_with_ref .and. opt_with_base) then
      if (tao_com%unified_lattices) then
        data(j)%delta_merit = model_value - &
            data(j)%meas_value + data(j)%ref_value - s%u(iu0)%data(j)%model_value
      else
        data(j)%delta_merit = model_value - &
            data(j)%meas_value + data(j)%ref_value - data(j)%base_value
      endif
    elseif (opt_with_ref) then
      data(j)%delta_merit = model_value - &
            data(j)%meas_value + data(j)%ref_value - data(j)%design_value
    elseif (opt_with_base) then
      if (tao_com%unified_lattices) then
        data(j)%delta_merit = model_value - &
                                data(j)%meas_value -  s%u(iu0)%data(j)%model_value
      else
        data(j)%delta_merit = model_value - &
                                data(j)%meas_value - data(j)%base_value
      endif
    else
      if (data(j)%merit_type(1:3) == 'int') then
        data(j)%delta_merit = data(j)%model_value
      else
        data(j)%delta_merit = model_value - data(j)%meas_value 
      endif
    endif
  enddo

! For phase data, since there is an arbitrary overall phase,
! we choose to make the average delta zero.

  call tao_find_data (err, 'phase.x', d1_ptr = d1, ix_uni = i, print_err = .false.)
  if (.not. err) then
    n = count(d1%d%useit_opt .and. d1%d%good_model)
    if (n /= 0 .and. all(d1%d%ele0_name == ' ')) then
      ave = sum(d1%d%delta_merit, mask = d1%d%useit_opt .and. d1%d%good_model) / n
      d1%d%delta_merit = d1%d%delta_merit - ave
    endif
  endif

  call tao_find_data (err, 'phase.y', d1_ptr = d1, ix_uni = i, print_err = .false.)
  if (.not. err) then
    n = count(d1%d%useit_opt .and. d1%d%good_model)
    if (n /= 0 .and. all(d1%d%ele0_name == ' ')) then
      ave = sum(d1%d%delta_merit, mask = d1%d%useit_opt) / n
      d1%d%delta_merit = d1%d%delta_merit - ave
    endif
  endif

! for max or min merit_types the delta might be modified.

  do j = 1, size(data)
    select case (data(j)%merit_type)
    case ('target', 'match', 'int_max', 'int_min')  ! Nothing to be done
    case ('max', 'abs_max')
      if (data(j)%delta_merit < 0) data(j)%delta_merit = 0  ! it's OK to be less
    case ('min', 'abs_min')
      if (data(j)%delta_merit > 0) data(j)%delta_merit = 0  ! it's OK to be more
    case default
      call tao_hook_merit_data (i, j, data(j))
    end select
  enddo

  where (data%useit_opt .and. data%good_model) data%merit = data%weight * data%delta_merit**2
  this_merit = this_merit + sum (data%merit, mask = data%useit_opt .and. data%good_model)

enddo

end function
