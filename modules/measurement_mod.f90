module measurement_mod

use bmad_struct
use bmad_interface

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_orbit_reading (orb, ele, axis, reading)
!
! Calculate the measured reading on a bpm given the actual orbit and the 
! BPM's offsets, noise, etc.
!
! This routine will only give a nonzero reading for Bmad markers, 
! monitors, and instruments.
!
! Modules needed:
!   use measurement_mod
!
! Input: 
!  orb        -- Coord_struct: Orbit position at BPM.
!  ele        -- Ele_struct: Element where the orbit is measured.
!    %value(noise$)         -- relative bpm resolution RMS
!    %value(tilt_tot$)      -- angle error in radians rms.
!    %value(x_gain_calib$)  -- Horizontal gain correction.
!    %value(y_gain_err$)    -- Horizontal gain error.
!    ... etc ...
!  axis       -- Integer: x_plane$ or y_plane$
!
! Output:
!  reading  -- Real(rp): BPM reading
!-

subroutine to_orbit_reading (orb, ele, axis, reading)

use random_mod

implicit none

type (coord_struct) orb
type (ele_struct) ele

real(rp) reading, noise_factor
real(rp) ran_num, angle, x, y, x_gain, y_gain

integer axis

character(20) :: r_name = "to_orbit_reading"

logical err

!

if (.not. ele%is_on) then
  reading = 0.0
  return
endif

if (ele%key /= monitor$ .and. ele%key /= instrument$ .and. ele%key /= marker$) then
  reading = 0.0
  return
endif

if (ele%value(noise$) /= 0) then
  call ran_gauss (ran_num)
  noise_factor = ele%value(noise$) * ran_num
else
  noise_factor = 0
endif

x_gain = 1 + ele%value(x_gain_err$) - ele%value(x_gain_calib$)
y_gain = 1 + ele%value(y_gain_err$) - ele%value(y_gain_calib$)
x = orb%vec(1) - ele%value(x_offset_tot$) + ele%value(x_offset_calib$)
y = orb%vec(3) - ele%value(y_offset_tot$) + ele%value(y_offset_calib$)

if (axis == x_plane$) then
  angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) + &
                (ele%value(crunch$) - ele%value(crunch_calib$))
  reading = noise_factor + x_gain * (x * cos(angle) + y * sin(angle))
elseif (axis == y_plane$) then
  angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) - &
                (ele%value(crunch$) - ele%value(crunch_calib$))
  reading = noise_factor + y_gain * (-x * sin(angle) + y * cos(angle))
else
  reading = 0.0
  call out_io (s_warn$, r_name, "This axis not supported for BPM reading!")
  call err_exit
endif

end subroutine to_orbit_reading

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_eta_reading (eta, ele, axis, reading)
!
! Compute the measured dispersion reading given the true dispersion and the 
! monitor offsets, noise, etc.
!
! This routine will only give a nonzero reading for Bmad markers, 
! monitors, and instruments.
!
! Modules needed:
!   use measurement_mod
!
! Input: 
!  eta_actual(2) -- Real(rp): Actual (eta_x, eta_y) dispersion.
!  ele           -- Ele_struct: Element where the orbit is measured.
!    %value(dE_eta_meas$)   -- Percent energy change used in dispersion measurement.
!    %value(noise$)         -- relative bpm resolution RMS
!    %value(tilt_tot$)      -- angle error in radians rms.
!    %value(x_gain_calib$)  -- Horizontal gain correction.
!    %value(y_gain_err$)    -- Horizontal gain error.
!    ... etc ...
!  axis       -- Integer: x_plane$ or y_plane$
!
! Output:
!  reading  -- Real(rp): BPM reading
!-

subroutine to_eta_reading (eta_actual, ele, axis, reading)

use random_mod

implicit none

type (ele_struct) ele

real(rp) eta_actual(:)
real(rp) reading, noise_factor
real(rp) ran_num, angle, x, y, x_gain, y_gain

integer axis

character(20) :: r_name = "to_eta_reading"

logical err

!

reading = 0.0

if (.not. ele%is_on) return
if (ele%key /= monitor$ .and. ele%key /= instrument$ .and. &
                                      ele%key /= marker$) return

if (ele%value(noise$) /= 0) then
  if (ele%value(de_eta_meas$) == 0) then
    call out_io (s_warn$, r_name, "dE_eta_meas not set for: " // ele%name)
    return
  endif
  call ran_gauss (ran_num)
  noise_factor = sqrt_2 * ele%value(noise$) * ran_num / ele%value(de_eta_meas$)
else
  noise_factor = 0
endif

x_gain = 1 + ele%value(x_gain_err$) - ele%value(x_gain_calib$)
y_gain = 1 + ele%value(y_gain_err$) - ele%value(y_gain_calib$)
x = eta_actual(1)
y = eta_actual(2)

if (axis == x_plane$) then
  angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) + &
                (ele%value(crunch$) - ele%value(crunch_calib$))
  reading = noise_factor + x_gain * (x * cos(angle) + y * sin(angle))
elseif (axis == y_plane$) then
  angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) - &
                (ele%value(crunch$) - ele%value(crunch_calib$))
  reading = noise_factor + y_gain * (-x * sin(angle) + y * cos(angle))
else
  reading = 0.0
  call out_io (s_warn$, r_name, "This axis not supported for BPM reading!")
  call err_exit
endif
  
end subroutine to_eta_reading

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_betatron_phase_reading (actual_phase, ele, reading)
!
! Find the measured betatron phase given the actual phase
!
! This routine will only give a nonzero reading for Bmad markers, 
! monitors, and instruments.
!
! Modules needed:
!   use measurement_mod
!
! Input: 
!  actual_phase -- Real(rp): Actual phase reading.
!  ele          -- Ele_struct: Element where phase is measured.
!    %value(phase_noise$) -- RMS Noise in radians.
!
! Output:
!  reading  -- Real(rp): Betatron phase reading
!-

subroutine to_betatron_phase_reading (actual_phase, ele, reading)

use random_mod

implicit none

type (ele_struct) ele

real(rp) reading, actual_phase
real(rp) ran_num

integer axis

character(40) :: r_name = "to_betatron_phase_reading"

logical err

if (.not. ele%is_on) then
  reading = 0.0
  return
endif

if (ele%key /= monitor$ .and. ele%key /= instrument$ .and. ele%key /= marker$) then
  reading = 0.0
  return
endif

reading = actual_phase

if (ele%value(noise$) /= 0) then
  call ran_gauss (ran_num)
  reading = reading + ele%value(noise$) * ran_num
endif

end subroutine to_betatron_phase_reading

end module
