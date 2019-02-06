!+
! Function ac_kicker_amp(ele, time) result (ac_amp)
!
! Routine to calculate the amplitude of the field for an ac_kicker element
!
! Input:
!   ele     -- ele_struct: ac_kicker element.
!   time    -- real(rp): Time to evaluate the amplitude at.
!
! Output:
!   ac_amp  -- real(rp): Amplitude. Will be set to 1 if the element is not an ac_kicker.
!-

function ac_kicker_amp(ele, time) result (ac_amp)

use bmad_interface, dummy => ac_kicker_amp

implicit none

type (ele_struct), target :: ele
type (ele_struct), pointer :: lord
type (ac_kicker_struct), pointer :: ac
real(rp) t, time, ac_amp, f
integer i, n, ix

character(*), parameter :: r_name = 'ac_kicker_amp'
!

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
  n = size(ac%amp_vs_time)
  call bracket_index(ac%amp_vs_time%time, 1, n, t, ix)
  if (ix < 1) then
    ac_amp = ac%amp_vs_time(1)%amp
  elseif (ix == n) then
    ac_amp = ac%amp_vs_time(n)%amp
  elseif (nint(ele%value(interpolation$)) == linear$) then
    f = (t - ac%amp_vs_time(ix)%time) / (ac%amp_vs_time(ix+1)%time - ac%amp_vs_time(ix)%time)
    ac_amp = ac%amp_vs_time(ix)%amp * (1 - f) + ac%amp_vs_time(ix+1)%amp * f
  elseif (nint(ele%value(interpolation$)) == spline$) then
    ac_amp = spline1 (ac%amp_vs_time(ix)%spline, t)
  else
    call out_io (s_fatal$, r_name, 'UNKNOWN INTERPOLATION TYPE FOR AC_KICKER: ' // ele%name)
    if (global_com%exit_on_error) call err_exit
    return
  endif
endif

end function
