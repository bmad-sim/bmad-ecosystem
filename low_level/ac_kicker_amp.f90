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

type (ele_struct) ele
type (ac_kicker_struct), pointer :: ac
real(rp) t, time, ac_amp
integer i, n, ix

!

ac_amp = 1
if (ele%key /= ac_kicker$) return

ac => ele%ac_kick
t = time - ele%value(t_offset$)

if (allocated(ac%frequencies)) then
  ac_amp = 0
  do i = 1, size(ac%frequencies)
    ac_amp = ac_amp + ac%frequencies(i)%amp * cos(twopi*(ac%frequencies(i)%f * t + ac%frequencies(i)%phi))
  enddo

else
  n = size(ac%amp_vs_time)
  call bracket_index(ac%amp_vs_time%time, 1, n, t, ix)
  if (ix < 1 .or. ix == n) then
    ac_amp = 0
  else
    ac_amp = spline1 (ac%amp_vs_time(ix)%spline, t)
  endif
endif

end function
