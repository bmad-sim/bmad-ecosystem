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

use tao_interface, dummy => tao_merit

implicit none

type (tao_var_struct), pointer :: var
type (tao_data_struct), pointer :: data(:)
type (tao_universe_struct), pointer :: u

real(rp) this_merit, value, model_value, meas_value, ref_value

integer i, j, k, n, iu0

character(*), parameter :: r_name = "tao_merit"

logical, optional :: calc_ok
logical ok


! make sure all calculations are up to date.

this_merit = 0

if (.not. s%global%lattice_calc_on) return

call tao_lattice_calc (ok)
if (present(calc_ok)) calc_ok = ok

!----------------------------------------
! Merit contribution from the variables.

do j = 1, s%n_var_used
  var => s%var(j)
  var%delta_merit = 0
  var%merit = 0

  if (.not. var%useit_opt) cycle

  ! Calculate the merit delta

  select case (var%merit_type)
  case ('target')
    if (s%global%opt_with_ref .and. s%global%opt_with_base) then
      var%delta_merit = (var%model_value - var%base_value) - (var%meas_value - var%ref_value)
    elseif (s%global%opt_with_ref) then
      var%delta_merit = (var%model_value - var%design_value) - (var%meas_value - var%ref_value)
    elseif (s%global%opt_with_base) then
      var%delta_merit = (var%model_value - var%base_value) - var%meas_value
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

!----------------------------------------------------------------------
! Merit contribution from the data:

do i = lbound(s%u, 1), ubound(s%u, 1)
  u => s%u(i)
  data => u%data
  data%merit = 0
  data%delta_merit = 0

  ! Scale ping data

  call tao_scale_ping_data(u)

  ! Check if universe is turned off

  if (.not. u%is_on) cycle

  ! First compute the delta for the merit function
  
  do j = 1, size(data)
    if (.not. data(j)%useit_opt) cycle

    if (.not. data(j)%good_model .or. (s%global%opt_with_ref .and. .not. data(j)%good_design) .or. &
                                      (s%global%opt_with_base .and. .not. data(j)%good_base)) then
      data(j)%delta_merit = data(j)%invalid_value
      cycle
    endif

    if (data(j)%merit_type(1:4) == 'abs_') then
      model_value = abs(data(j)%model_value)
    else
      model_value = data(j)%model_value
    endif

    meas_value = data(j)%meas_value
    ref_value = data(j)%ref_value

    if (data(j)%d1%d2%name == 'ping_a' .and. data(j)%d1%name(1:5) /= 'phase') then
      meas_value = meas_value * u%ping_scale%a_mode_meas
      ref_value  = ref_value  * u%ping_scale%a_mode_ref
    endif

    if (data(j)%d1%d2%name == 'ping_b' .and. data(j)%d1%name(1:5) /= 'phase') then
      meas_value = meas_value * u%ping_scale%b_mode_meas
      ref_value  = ref_value  * u%ping_scale%b_mode_ref
    endif

    if (s%global%opt_with_ref .and. s%global%opt_with_base) then
      data(j)%delta_merit = (model_value - data(j)%base_value) - (meas_value - ref_value) 
    elseif (s%global%opt_with_ref) then
      data(j)%delta_merit = (model_value - data(j)%design_value) - (meas_value - ref_value)
    elseif (s%global%opt_with_base) then
      data(j)%delta_merit = (model_value - data(j)%base_value) - meas_value
    else
      data(j)%delta_merit = model_value - meas_value 
    endif
  enddo

  ! For phase data, since there is an arbitrary overall phase,
  ! we choose to make the average delta zero.

  call zero_phase_merit ('phase.a', i)
  call zero_phase_merit ('phase.b', i)
  call zero_phase_merit ('bpm_phase.a', i)
  call zero_phase_merit ('bpm_phase.b', i)
  call zero_phase_merit ('ping_a.phase_x', i)
  call zero_phase_merit ('ping_a.phase_y', i)
  call zero_phase_merit ('ping_b.phase_x', i)
  call zero_phase_merit ('ping_b.phase_y', i)

  ! for max or min merit_types the delta might be modified.

  do j = 1, size(data)
    if (.not. data(j)%exists) cycle
    select case (data(j)%merit_type)
    case ('target', 'average', 'integral', 'max-min', 'rms')
      ! Nothing to be done
    case ('max', 'abs_max')
      if (data(j)%delta_merit < 0) data(j)%delta_merit = 0  ! it's OK to be less
    case ('min', 'abs_min')
      if (data(j)%delta_merit > 0) data(j)%delta_merit = 0  ! it's OK to be more
    case default
      call tao_hook_merit_data (i, j, data(j), ok)
      if (.not. ok) then
        call out_io (s_error$, r_name, 'MERIT_TYPE NOT RECOGNIZED FOR DATA: ' // tao_datum_name(data(j)), &
                                       'MERIT_TYPE: ' // data(j)%merit_type, &
                                       'WILL MARK DATUM AS DOES NOT EXIST.')
        data%exists = .false.
      endif
    end select
  enddo

  where (data%useit_opt) data%merit = data%weight * data%delta_merit**2
  this_merit = this_merit + sum (data%merit, mask = data%useit_opt)
enddo

!-------------------------------------------------------------------------------------------
contains

subroutine zero_phase_merit (who, ix_uni)

type (tao_d1_data_array_struct), allocatable :: d1_array(:)
character(*) who
real(rp) ave
integer ix_uni, n
logical err

!

call tao_find_data (err, who, d1_array = d1_array, ix_uni = ix_uni, print_err = .false.)
if (err) return
if (size(d1_array) == 0) return
n = count(d1_array(1)%d1%d%useit_opt .and. d1_array(1)%d1%d%good_model)
if (n /= 0 .and. all(d1_array(1)%d1%d%ele_ref_name == '')) then
  ave = sum(d1_array(1)%d1%d%delta_merit, &
                  mask = d1_array(1)%d1%d%useit_opt .and. d1_array(1)%d1%d%good_model) / n
  d1_array(1)%d1%d%delta_merit = d1_array(1)%d1%d%delta_merit - ave
endif

end subroutine

end function
