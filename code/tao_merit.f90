!+
! function tao_merit (s) result (this_merit)
! 
! function to calculate the merit.
!
! Input:
!   s       -- Tao_super_universe_struct:
!
! Output:
!   this_merit -- Real(rp): Merit value.
!-

function tao_merit (s) result (this_merit)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_var_struct), pointer :: var
type (tao_data_struct), pointer :: data(:)
type (tao_d1_data_struct), pointer :: d1

real(rp) this_merit, ave, value

integer i, j, n
logical err

! make sure all calculations are up to date.

call tao_lattice_calc (s)

!----------------------------------------
! Merit contribution from the variables.

this_merit = 0

do i = 1, size(s%u)
  do j = 1, size(s%u(i)%var)

    var => s%u(i)%var(j)
    if (.not. var%useit_opt) cycle

    if (s%global%opt_with_ref .and. s%global%opt_with_base) then
      value = var%model_value - var%base_value + var%ref_value
    elseif (s%global%opt_with_ref) then
      value = var%model_value - var%design_value + var%ref_value
    elseif (s%global%opt_with_base) then
      value = var%model_value - var%base_value
    else
      value = var%model_value 
    endif
    
    var%merit = 0
    select case (var%merit_type)
    case ('target')
      var%merit = var%weight * (value - var%data_value)**2
    case ('limit')
      if (var%model_value > var%high_lim_value) then
        var%merit = var%weight * (value - var%high_lim_value)**2
      elseif (var%model_value < var%low_lim_value) then
        var%merit = var%weight * (value - var%low_lim_value)**2
      endif
    case default
      call tao_hook_merit_var (s, i, j, var)
    end select

    this_merit = this_merit + var%merit

  enddo
enddo

!----------------------------------------
! Merit contribution from the data:

do i = 1, size(s%u)

! First compute the delta for the merit function
  
  data => s%u(i)%data
  data%merit = 0
  data%delta = 0

  if (s%global%opt_with_ref .and. s%global%opt_with_base) then
    where (data%useit_opt) data%delta = data%model_value - &
            data%data_value + data%ref_value - data%base_value
  elseif (s%global%opt_with_ref) then
    where (data%useit_opt) data%delta = data%model_value - &
            data%data_value + data%ref_value - data%design_value
  elseif (s%global%opt_with_base) then
    where (data%useit_opt) data%delta = data%model_value - &
                                data%data_value - data%base_value
  else
    where (data%useit_opt) data%delta = data%model_value - data%data_value 
  endif

! For phase data, since there is an arbitrary overall phase,
! we choose to make the average delta zero.

  call tao_find_data (err, s%u(i), 'phase:x', d1_ptr = d1)
  if (.not. err) then
    n = count(d1%d%useit_opt)
    if (n /= 0) then
      ave = sum(d1%d%delta, mask = d1%d%useit_opt) / n
      d1%d%delta = d1%d%delta - ave
    endif
  endif

  call tao_find_data (err, s%u(i), 'phase:y', d1_ptr = d1)
  if (.not. err) then
    n = count(d1%d%useit_opt)
    if (n /= 0) then
      ave = sum(d1%d%delta, mask = d1%d%useit_opt) / count(d1%d%useit_opt)
      d1%d%delta = d1%d%delta - ave
    endif
  endif

! for max or min merit_types the delta might be modified.

  do j = 1, size(data)
    select case (data(j)%merit_type)
    case ('target')
    case ('max')
      if (data(j)%delta < 0) data(j)%delta = 0  ! it's OK to be less
    case ('min')
      if (data(j)%delta > 0) data(j)%delta = 0  ! it's OK to be more
    case default
      call tao_hook_merit_data (s, i, j, data(j))
    end select
  enddo

  where (data%useit_opt) data%merit = data%weight * data%delta**2
  this_merit = this_merit + sum (data%merit, mask = data%useit_opt)

enddo

end function
