module measurement_mod

use bmad_routine_interface

private compute_measurement_distortion_mat

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Function ele_is_monitor (ele, print_warning) result (is_monitor)
!
! Routine to check that an element is either a detector, instrument, monitor, or marker. 
! These are the elements where measurement errors can be defined.
!
! Input:
!   ele           -- ele_struct: Lattice element.
!   print_warning -- logical, optional: If True print a warning message if the
!                     element not a monitor like element. Default is True.
!
! Output:
!   is_monitor    -- logical: Set True if the element is a monitor like element.
!-

function ele_is_monitor (ele, print_warning) result (is_monitor)

implicit none

type (ele_struct) ele
logical, optional :: print_warning
logical is_monitor
character(*), parameter :: r_name = 'ele_is_monitor'

!

select case (ele%key)
case (monitor$, instrument$, marker$, detector$) 
  is_monitor = .true.
case default
  if (logic_option(.true., print_warning)) then
    call out_io (s_error$, r_name, &
                'MONITOR CALCULATION CALLED FOR ELEMENT THAT IS NEITHER', &
                'A MONITOR, INSTRUMENT, MARKER, NOR DETECTOR: ' // ele%name)
  endif
  is_monitor = .false.
end select

end function ele_is_monitor

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine compute_measurement_distortion_mat (ele, d_mat)
!
! Routine to compute the matrix associated with the transformation between
! the actual orbit, phase, eta, and coupling, and what the measured values are.
!
! This routine is private and not meant for general use.
!
! Input:
!   ele        -- Ele_struct: An element at which measurements are made.
!
! Output:
!   d_mat(2,2) -- real(rp): Distortion matrix.

subroutine compute_measurement_distortion_mat (ele, d_mat)

implicit none

type (ele_struct) ele

real(rp) d_mat(2,2), x_gain, y_gain, x_angle, y_angle

!

x_gain = 1 + ele%value(x_gain_err$) - ele%value(x_gain_calib$)
y_gain = 1 + ele%value(y_gain_err$) - ele%value(y_gain_calib$)
x_angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) + (ele%value(crunch$) - ele%value(crunch_calib$))
y_angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) - (ele%value(crunch$) - ele%value(crunch_calib$))

if (x_angle == 0 .and. y_angle == 0) then
  d_mat(1,1) = x_gain
  d_mat(1,2) = 0
  d_mat(2,1) = 0
  d_mat(2,2) = y_gain
else
  d_mat(1,1) =  x_gain * cos(x_angle)
  d_mat(1,2) =  x_gain * sin(x_angle)
  d_mat(2,1) = -y_gain * sin(y_angle)
  d_mat(2,2) =  y_gain * cos(y_angle)
endif

end subroutine compute_measurement_distortion_mat

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_orbit_reading (orb, ele, axis, add_noise, reading, err)
!
! Calculate the measured reading on a bpm given the actual orbit and the 
! BPM's offsets, noise, etc.
!
! This routine will only give a nonzero reading for Bmad markers, 
! monitors, and instruments.
!
! Input: 
!   orb        -- Coord_struct: Orbit position at BPM.
!   ele        -- Ele_struct: Element where the orbit is measured.
!    %value(noise$)         -- relative bpm resolution RMS
!    %value(tilt_tot$)      -- angle error in radians rms.
!    %value(x_gain_calib$)  -- Horizontal gain correction.
!    %value(y_gain_err$)    -- Horizontal gain error.
!    ... etc ...
!   axis       -- Integer: x_plane$ or y_plane$
!   add_noise  -- logical: If True add noise to the reading
!
! Output:
!  reading  -- Real(rp): BPM reading
!  err      -- Logical: Set True if there is no valid reading.
!                For example, if ele%is_on = False.
!-

subroutine to_orbit_reading (orb, ele, axis, add_noise, reading, err)

use random_mod

implicit none

type (coord_struct) orb
type (ele_struct) ele

real(rp) d_mat(2,2), reading, noise_factor
real(rp) ran_num, x, y

integer axis

character(20) :: r_name = "to_orbit_reading"

logical add_noise, err, error

!

reading = 0.0
err = .true.

if (.not. ele_is_monitor(ele)) return
if (.not. ele%is_on) return

call compute_measurement_distortion_mat (ele, d_mat)

if (add_noise .and. ele%value(noise$) /= 0) then
  call ran_gauss (ran_num)
  noise_factor = ele%value(noise$) * ran_num
else
  noise_factor = 0
endif

x = orb%vec(1) - (ele%value(x_offset_tot$) + ele%value(x_offset_calib$))
y = orb%vec(3) - (ele%value(y_offset_tot$) + ele%value(y_offset_calib$))

select case (axis)
case (x_plane$)
  reading = noise_factor + (x * d_mat(1,1) + y * d_mat(1,2))
case (y_plane$)
  reading = noise_factor + (x * d_mat(2,1) + y * d_mat(2,2))
case default
  reading = 0.0
  call out_io (s_warn$, r_name, "This axis not supported for BPM reading!")
  if (global_com%exit_on_error) call err_exit
end select

err = .false.

end subroutine to_orbit_reading

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_eta_reading (eta, ele, axis, add_noise, reading, err)
!
! Compute the measured dispersion reading given the true dispersion and the 
! monitor offsets, noise, etc.
!
! This routine will only give a nonzero reading for Bmad markers, 
! monitors, and instruments.
!
! Input: 
!   eta_actual(2) -- Real(rp): Actual (eta_x, eta_y) dispersion.
!   ele           -- Ele_struct: Element where the orbit is measured.
!     %value(dE_eta_meas$)   -- Percent energy change used in dispersion measurement.
!     %value(noise$)         -- relative bpm resolution RMS
!     %value(tilt_tot$)      -- angle error in radians rms.
!     %value(x_gain_calib$)  -- Horizontal gain correction.
!     %value(y_gain_err$)    -- Horizontal gain error.
!     ... etc ...
!   axis          -- Integer: x_plane$ or y_plane$
!   add_noise     -- logical: If True add noise to the reading
!
! Output:
!  reading  -- Real(rp): BPM reading
!  err      -- Logical: Set True if there is an error. False otherwise.
!-

subroutine to_eta_reading (eta_actual, ele, axis, add_noise, reading, err)

use random_mod

implicit none

type (ele_struct) ele

real(rp) eta_actual(:)
real(rp) d_mat(2,2), reading, noise_factor
real(rp) ran_num, x, y

integer axis

character(20) :: r_name = "to_eta_reading"

logical add_noise, err, error

reading = 0.0

!

err = .true.

if (.not. ele%is_on) return
if (.not. ele_is_monitor(ele)) return

if (add_noise .and. ele%value(noise$) /= 0) then
  if (ele%value(de_eta_meas$) == 0) then
    call out_io (s_warn$, r_name, "dE_eta_meas not set for: " // ele%name)
    return
  endif
  call ran_gauss (ran_num)
  noise_factor = sqrt_2 * ele%value(noise$) * ran_num / ele%value(de_eta_meas$)
else
  noise_factor = 0
endif

call compute_measurement_distortion_mat (ele, d_mat)

x = eta_actual(1) - (ele%value(x_dispersion_err$) + ele%value(x_dispersion_calib$))
y = eta_actual(2) - (ele%value(y_dispersion_err$) + ele%value(y_dispersion_calib$))

if (axis == x_plane$) then
  reading = noise_factor + (x * d_mat(1,1) + y * d_mat(1,2))
elseif (axis == y_plane$) then
  reading = noise_factor + (x * d_mat(2,1) + y * d_mat(2,2))
else
  reading = 0.0
  call out_io (s_warn$, r_name, "This axis not supported for BPM reading!")
  if (global_com%exit_on_error) call err_exit
endif

err = .false.
  
end subroutine to_eta_reading

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_phase_and_coupling_reading (ele, add_noise, reading, err)
!
! Find the measured coupling values given the actual ones
!
! This routine will only give a nonzero reading for Bmad markers, 
! monitors, and instruments.
!
! Input: 
!   actual_phase  -- Real(rp): Actual phase reading.
!   ele           -- Ele_struct: Element where phase is measured.
!     %value(phase_noise$) -- RMS Noise in radians.
!   add_noise     -- logical: If True add noise to the reading
!
! Output:
!   reading   -- bpm_phase_coupling_struct: K and Cbar coupling parameters
!   err       -- Logical: Set True if there is an error. False otherwise.
!-

subroutine to_phase_and_coupling_reading (ele, add_noise, reading, err)

use random_mod

implicit none

type (ele_struct) ele
type (bpm_phase_coupling_struct) reading, lab

real(rp) d_mat(2,2), ran_num(6), cbar_lab(2,2), denom, ratio
real(rp) q_a1_lab(2), q_a2_lab(2), q_b1_lab(2), q_b2_lab(2)
real(rp) q_a1_mon(2), q_a2_mon(2), q_b1_mon(2), q_b2_mon(2)

integer axis

character(40) :: r_name = "to_phase_and_coupling_reading"

logical add_noise, err, error

! 

reading%K_22a = 0
reading%K_12a = 0
reading%K_11b = 0
reading%K_12b = 0

!

err = .true.

if (.not. ele%is_on) return
if (.not. ele_is_monitor(ele)) return
if (ele%a%beta == 0) return   ! Can happen if lattice is unstable

call compute_measurement_distortion_mat (ele, d_mat)

if (add_noise .and. ele%value(noise$) /= 0) then
  if (ele%value(n_sample$) == 0 .or. ele%value(osc_amplitude$) == 0) then
    call out_io (s_error$, r_name, 'N_SAMPLE OR OSC_AMPLITUDE NOT SET!')
    return
  endif
  call ran_gauss (ran_num)
  ran_num = ran_num * ele%value(noise$) / (ele%value(n_sample$) * ele%value(osc_amplitude$))
endif

!

call c_to_cbar (ele, cbar_lab)
ratio = sqrt(ele%b%beta / ele%a%beta)

! cbar_lab to k_lab

lab%K_22a = -ratio * cbar_lab(2,2) / ele%gamma_c
lab%K_12a = -ratio * cbar_lab(1,2) / ele%gamma_c
lab%K_11b =  cbar_lab(1,1) / (ratio * ele%gamma_c)
lab%K_12b = -cbar_lab(1,2) / (ratio * ele%gamma_c)

! k_lab to q_lab

q_a1_lab(1) = 1
q_a1_lab(2) = lab%K_22a

q_a2_lab(1) = 0
q_a2_lab(2) = lab%K_12a

q_b1_lab(1) = lab%K_11b
q_b1_lab(2) = 1

q_b2_lab(1) = lab%K_12b
q_b2_lab(2) = 0

! q_lab to q_mon

q_a1_mon = matmul(d_mat, q_a1_lab)
q_a2_mon = matmul(d_mat, q_a2_lab)
q_b1_mon = matmul(d_mat, q_b1_lab)
q_b2_mon = matmul(d_mat, q_b2_lab)

! q_mon to k_mon

denom = q_a1_mon(1)**2 + q_a2_mon(1)**2
reading%K_22a = (q_a1_mon(1) * q_a1_mon(2) + q_a2_mon(1) * q_a2_mon(2)) / denom
reading%K_12a = (q_a1_mon(1) * q_a2_mon(2) - q_a2_mon(1) * q_a1_mon(2)) / denom

denom = q_b1_mon(2)**2 + q_b2_mon(2)**2
reading%K_11b = (q_b1_mon(2) * q_b1_mon(1) + q_b2_mon(2) * q_b2_mon(1)) / denom
reading%K_12b = (q_b1_mon(2) * q_b2_mon(1) - q_b2_mon(2) * q_b1_mon(1)) / denom

! Add random terms to k_mon

if (ele%value(noise$) /= 0) then
  reading%K_22a = reading%K_22a + ran_num(1)
  reading%K_12a = reading%K_12a + ran_num(2)
  reading%K_11b = reading%K_11b + ran_num(3)
  reading%K_12b = reading%K_12b + ran_num(4)
endif

! k_mon to Cbar_mon

reading%cbar22_a = -reading%K_22a * ele%gamma_c / ratio 
reading%cbar12_a = -reading%K_12a * ele%gamma_c / ratio 
reading%cbar11_b =  reading%K_11b * ele%gamma_c * ratio
reading%cbar12_b = -reading%K_12b * ele%gamma_c * ratio

! Phase calc

reading%phi_a = ele%a%phi - atan2(q_a2_mon(1), q_a1_mon(1)) + ran_num(5)
reading%phi_b = ele%b%phi - atan2(q_b2_mon(2), q_b1_mon(2)) + ran_num(6)

err = .false.

end subroutine to_phase_and_coupling_reading

end module
