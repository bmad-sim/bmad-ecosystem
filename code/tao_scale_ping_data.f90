subroutine tao_scale_ping_data(u)

use tao_mod, except_dummy => tao_scale_ping_data

implicit none

type (tao_universe_struct) u
type (tao_d1_data_array_struct), allocatable, target :: d1_arr(:)
type (tao_d1_data_struct), pointer :: d1
type (tao_data_struct), pointer :: d

real(rp) model_m_sum, model_r_sum, meas_sum, ref_sum

integer i, j

logical err

!

call tao_find_data (err, 'ping_a.amp_x', d1_array = d1_arr, ix_uni = u%ix_uni, print_err = .false.)
call scale_this_ping(u%ping_scale%a_mode_meas, u%ping_scale%a_mode_ref)

call tao_find_data (err, 'ping_b.amp_y', d1_array = d1_arr, ix_uni = u%ix_uni, print_err = .false.)
call scale_this_ping(u%ping_scale%b_mode_meas, u%ping_scale%b_mode_ref)

!-------------------------------------
contains

subroutine scale_this_ping (meas_scale, ref_scale)

real(rp) meas_scale, ref_scale

!

do i = 1, size(d1_arr)
  model_m_sum = 0
  model_r_sum = 0
  meas_sum = 0
  ref_sum = 0

  d1 => d1_arr(i)%d1
  do j = lbound(d1%d, 1), ubound(d1%d, 1)
    d => d1%d(j)

    if (d%exists .and. d%good_user .and. d%good_model .and. d%good_meas) then
      model_m_sum = model_m_sum + d%model_value
      meas_sum  = meas_sum + d%meas_value
    endif

    if (d%exists .and. d%good_user .and. d%good_model .and. d%good_ref) then
      model_r_sum = model_r_sum + d%model_value
      ref_sum   = ref_sum + d%ref_value
    endif
  enddo

  if (meas_sum /= 0) meas_scale = model_m_sum / meas_sum
  if (ref_sum /= 0)  ref_scale  = model_r_sum / ref_sum
enddo
  
end subroutine scale_this_ping

end subroutine
