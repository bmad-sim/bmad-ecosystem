subroutine tao_scale_ping_data(u)

use tao_interface, except_dummy => tao_scale_ping_data

implicit none

type (tao_universe_struct) u

character(40), parameter :: ping_a_name(7) = [character(40):: &
            'ping_a.phase_x', 'ping_a.amp_y', 'ping_a.phase_y', 'ping_a.amp_sin_y', &
            'ping_a.amp_cos_y', 'ping_a.amp_sin_rel_y', 'ping_a.amp_cos_rel_y']

character(40), parameter :: ping_b_name(7) = [character(40):: &
            'ping_b.phase_y', 'ping_b.amp_x', 'ping_b.phase_x', 'ping_b.amp_sin_x', &
            'ping_b.amp_cos_x', 'ping_b.amp_sin_rel_x', 'ping_b.amp_cos_rel_x']

character(*), parameter :: r_name = 'tao_scale_ping_data'

!

call scale_this_ping(u, 'ping_a.amp_x', ping_a_name, u%ping_scale%a_mode_meas, u%ping_scale%a_mode_ref)
call scale_this_ping(u, 'ping_b.amp_y', ping_b_name, u%ping_scale%b_mode_meas, u%ping_scale%b_mode_ref)

!-------------------------------------
contains

subroutine scale_this_ping (u, ping_main, ping_name, meas_scale, ref_scale)

type (tao_universe_struct) u
type (tao_d1_data_array_struct), allocatable, target :: d1_arr(:), d1_arr2(:)
type (tao_d1_data_struct), pointer :: d1, d1b
type (tao_data_struct), pointer :: d

real(rp) meas_scale, ref_scale
real(rp) model_m_sum, model_r_sum, meas_sum, ref_sum

integer i, j, id

logical, allocatable :: good_user(:), good_user2(:)
logical err

character(*) ping_main, ping_name(:)
character(40) ping_here_name

! Potential problem is that the ping_a.amp_x or ping_b.amp_y data may be %good_user vetoed for optimization 
! but %good_user needs to be considered in case bad data points have been vetoed. 
! The solution is to consider other ping data types in this case.

call tao_find_data (err, ping_main, d1_array = d1_arr, ix_uni = u%ix_uni, print_err = .false.)
if (size(d1_arr) == 0) return
if (size(d1_arr) > 1) then
  call out_io (s_error$, r_name, 'MORE THAN ONE ' // trim(ping_main) // ' DATA ARRAY. I CANNOT HANDLE THIS.')
  return
endif

d1 => d1_arr(1)%d1
allocate (good_user(lbound(d1%d, 1):ubound(d1%d, 1)), good_user2(lbound(d1%d, 1):ubound(d1%d, 1)))
forall (id = lbound(d1%d, 1): ubound(d1%d, 1)) good_user(id) = d1%d(id)%good_user

! If all ping_main data is vetoed then look to another ping array for good_user values.

if (.not. any(good_user)) then
  do i = 1, size(ping_name)
    call tao_find_data (err, ping_name(i), d1_array = d1_arr2, ix_uni = u%ix_uni, print_err = .false.)
    do j = 1, size(d1_arr2)
      d1b => d1_arr2(j)%d1
      if (lbound(d1b%d, 1) /= lbound(d1%d, 1) .or. ubound(d1b%d, 1) /= ubound(d1%d, 1)) then
        call out_io(s_error$, r_name, 'DATA ARRAY ' // trim(ping_name(i)) // ' ARRAY BOUNDS DO NOT MATCH ' // ping_main)
        cycle
      endif
      forall (id = lbound(d1%d, 1): ubound(d1%d, 1)) good_user2(id) = d1b%d(id)%good_user
      if (.not. any(good_user2)) cycle
      if (any(good_user) .and. any(good_user .neqv. good_user2)) then
        call out_io(s_error$, r_name, 'DATA ARRAY ' // ping_name(i) // ' HAS GOOD_USER MISMATCH WITH ' // trim(ping_here_name))
      endif          
      good_user = good_user2
      ping_here_name = ping_name(i)
    enddo
  enddo
endif

! Now compute the scale factor

model_m_sum = 0
model_r_sum = 0
meas_sum = 0
ref_sum = 0

do j = lbound(d1%d, 1), ubound(d1%d, 1)
  d => d1%d(j)

  if (d%exists .and. good_user(j) .and. d%good_model .and. d%good_meas) then
    model_m_sum = model_m_sum + d%model_value
    meas_sum  = meas_sum + d%meas_value
  endif

  if (d%exists .and. good_user(j) .and. d%good_model .and. d%good_ref) then
    model_r_sum = model_r_sum + d%model_value
    ref_sum   = ref_sum + d%ref_value
  endif
enddo

if (meas_sum /= 0) meas_scale = model_m_sum / meas_sum
if (ref_sum /= 0)  ref_scale  = model_r_sum / ref_sum
  
end subroutine scale_this_ping

end subroutine
