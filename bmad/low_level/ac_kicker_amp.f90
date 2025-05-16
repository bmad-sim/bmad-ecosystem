!+
! Function ac_kicker_amp(ele, orbit, true_time) result (ac_amp)
!
! Routine to calculate the amplitude of the field for an ac_kicker element
!
! Input:
!   ele       -- ele_struct: ac_kicker element.
!   orbit     -- coord_struct: Contains the time to evaluate the amplitude at.
!   true_time -- real(rp), optional: The actual time. Normally this time is calculated using
!                     orbit%t or orbit%vec(5) but sometimes it is convenient to be able to override this.
!                     For example, time_runge_kutta uses this.
!
! Output:
!   ac_amp  -- real(rp): Amplitude. Will be set to 1 if the element is not an ac_kicker.
!-

function ac_kicker_amp(ele, orbit, true_time) result (ac_amp)

use bmad_interface, dummy => ac_kicker_amp

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (ele_struct), pointer :: lord
type (ac_kicker_struct), pointer :: ac
type (ele_struct), pointer :: ref_ele
type (ele_pointer_struct), allocatable :: chain(:)

real(rp), optional :: true_time
real(rp) t, ac_amp, dt_ds0
integer i, ix_pass, n_links
logical err_flag

character(*), parameter :: r_name = 'ac_kicker_amp'

!

ref_ele => ele
if (ref_ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) ref_ele => pointer_to_super_lord (ref_ele)

ac_amp = 1
if (ele%key /= ac_kicker$) return

! Slice slaves and super slaves have their associated %ac_kick components stored in the lord

if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
  lord => pointer_to_lord(ele, 1)
else
  lord => ele
endif

ac => lord%ac_kick

if (allocated(ac%frequency)) then
  ac_amp = 0
  do i = 1, size(ac%frequency)
    t = real_option(particle_rf_time(orbit, ele, rf_clock_harmonic = ac%frequency(i)%rf_clock_harmonic), true_time)
    ac_amp = ac_amp + ac%frequency(i)%amp * cos(twopi*(ac%frequency(i)%f * t + ac%frequency(i)%phi))
  enddo

else
  t = real_option(particle_rf_time(orbit, ele), true_time)
  ac_amp = knot_interpolate(ac%amp_vs_time%time, ac%amp_vs_time%amp, t, nint(ele%value(interpolation$)), err_flag)
  if (err_flag) then
    call out_io (s_fatal$, r_name, 'INTERPOLATION PROBLEM FOR AC_KICKER: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

end function
