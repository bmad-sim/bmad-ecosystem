!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine tao_set_var_model_value (var, value, print_limit_warning)
!
! Subroutine to set the value for a model variable and do the necessary bookkeeping.
! If value is past the variable's limit, and s%global%var_limits_on = True, the
! variable will be set to the limit value.
!
! Input:
!   var   -- Tao_var_struct: Variable to set
!   value -- Real(rp): Value to set to
!   print_limit_warning
!         -- Logical, optional: Print a warning if the value is past the variable's limits.
!             Default is True.
!-

subroutine tao_set_var_model_value (var, value, print_limit_warning)

use tao_interface, dummy => tao_set_var_model_value
use bookkeeper_mod, only: set_flags_for_changed_attribute

implicit none

type (tao_var_struct), target :: var
type (tao_universe_struct), pointer :: u
type (lat_struct), pointer :: lat
type (tao_var_slave_struct), pointer :: var_slave
type (ele_struct), pointer :: ele

real(rp) value
integer i
logical, optional :: print_limit_warning

!

if (.not. var%exists) return

! check if hit variable limit
if (s%global%var_limits_on .and. (.not. s%global%only_limit_opt_vars .or. var%useit_opt)) then
  if (value < var%low_lim) then
    if (logic_option (.true., print_limit_warning)) &
          call out_io (s_blank$, ' ', "Hit lower limit of variable: " // tao_var1_name(var))
    value = var%low_lim
  elseif (value > var%high_lim) then
    if (logic_option (.true., print_limit_warning)) &
          call out_io (s_blank$, ' ', "Hit upper limit of variable: " // tao_var1_name(var))
    value = var%high_lim
  endif
endif

var%model_value = value
do i = 1, size(var%slave)
  var_slave => var%slave(i)
  var_slave%model_value = value
  u => s%u(var_slave%ix_uni)

  u%calc%lattice = .true.

  lat => u%model%lat
  if (var%ele_name /= 'PARTICLE_START') then
    ele => lat%branch(var_slave%ix_branch)%ele(var_slave%ix_ele)
    call set_flags_for_changed_attribute (ele, var_slave%model_value)
  endif
enddo

end subroutine tao_set_var_model_value
