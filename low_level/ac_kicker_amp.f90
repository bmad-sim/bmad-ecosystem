!+
! Function ac_kicker_amp(ele, time) result (ac_amp)
!
! Routine to calculate the amplitude of the field for an ac_kicker element
!
! Input:
!   ele     -- ele_struct: ac_kicker element.
!   orbit   -- coord_struct: Contains the time to evaluate the amplitude at.
!
! Output:
!   ac_amp  -- real(rp): Amplitude. Will be set to 1 if the element is not an ac_kicker.
!-

function ac_kicker_amp(ele, orbit) result (ac_amp)

use bmad_interface, dummy => ac_kicker_amp

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (ele_struct), pointer :: lord
type (ac_kicker_struct), pointer :: ac
type (ele_struct), pointer :: ref_ele
type (ele_pointer_struct), allocatable :: chain(:)

real(rp) t, time, ac_amp, dt_ds0
integer i, ix_pass, n_links
logical err_flag

character(*), parameter :: r_name = 'ac_kicker_amp'

!

ref_ele => ele
if (ref_ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) ref_ele => pointer_to_lord (ref_ele, 1)

if (is_true(ele%value(ref_time_offset$))) then
  time = orbit%t - ref_ele%value(ref_time_start$)
else
  time = orbit%t
endif

ac_amp = 1
if (ele%key /= ac_kicker$) return

! Slice slaves and super slaves have their associated %ac_kick components stored in the lord

if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
  lord => pointer_to_lord(ele, 1)
else
  lord => ele
endif

ac => lord%ac_kick
t = time - ele%value(t_offset$)

if (allocated(ac%frequencies)) then
  ac_amp = 0
  do i = 1, size(ac%frequencies)
    ac_amp = ac_amp + ac%frequencies(i)%amp * cos(twopi*(ac%frequencies(i)%f * t + ac%frequencies(i)%phi))
  enddo

else
  ac_amp = knot_interpolate(ac%amp_vs_time%time, ac%amp_vs_time%amp, t, nint(ele%value(interpolation$)), err_flag)
  if (err_flag) then
    call out_io (s_fatal$, r_name, 'INTERPOLATION PROBLEM FOR AC_KICKER: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

end function
