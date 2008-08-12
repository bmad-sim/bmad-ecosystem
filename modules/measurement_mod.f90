module measurement_mod

use bmad_struct
use bmad_interface

type measurement_common_struct
  real(rp) :: x_off, y_off
  real(rp) :: M_m(2,2)
  real(rp) :: sqrt_beta_a, sqrt_beta_b
  real(rp) :: beta_a = 0, beta_b = 0
  real(rp) :: value(n_attrib_maxx) = real_garbage$
end type

type (measurement_common_struct), private, save :: m_com

private compute_bpm_transformation_numbers

contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine check_if_ele_is_monitor (ele, err)
!
! Routine to check that the element is either an instrument, monitor, or marker.
! This routine is private and not meant for general use.
!-

subroutine check_if_ele_is_monitor (ele, err)

implicit none

type (ele_struct) ele
logical err
character(40) :: r_name = 'check_if_ele_is_monitor'

!

if (ele%key /= monitor$ .and. ele%key /= instrument$ .and. ele%key /= marker$) then
  call out_io (s_error$, r_name, &
                'MONITOR CALCULATION CALLED FOR ELEMENT THAT IS NEITHER', &
                'A MONITOR, INSTRUMENT OR MARKER')
  err = .true.
  return
else
  err = .false.
endif

end subroutine

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine compute_bpm_transformation_numbers (ele)
!
! Routine to compute the numbers associated with the transformation between
! the actual orbit, phase, eta, and coupling, and what the measured values.
!
! This routine is private and not meant for general use.
!
! Input:
!   ele -- Ele_struct: An element at which measurements are made.
!     %value(x_offset_calib$) -- Offset calibration value.
!     ... etc ...
!
! Output:
!   m_com -- Private common block of measurement mod.

subroutine compute_bpm_transformation_numbers (ele)

implicit none

type (ele_struct) ele

real(rp) x_gain, y_gain, x_angle, y_angle
integer, parameter :: ix_attribs(12) = (/ x_gain_err$, x_gain_calib$, &
          y_gain_err$, y_gain_calib$, tilt_tot$, tilt_calib$, crunch$, crunch_calib$, &
          x_offset_tot$, x_offset_calib$, y_offset_tot$, y_offset_calib$ /)

!

if (m_com%beta_a /= ele%a%beta .or. m_com%beta_b /= ele%b%beta) then
  m_com%sqrt_beta_a = sqrt(ele%a%beta)
  m_com%sqrt_beta_b = sqrt(ele%b%beta)
  m_com%beta_a = ele%a%beta
  m_com%beta_b = ele%b%beta
endif

if (any(m_com%value(ix_attribs) /= ele%value(ix_attribs))) then
  x_gain = 1 + ele%value(x_gain_err$) - ele%value(x_gain_calib$)
  y_gain = 1 + ele%value(y_gain_err$) - ele%value(y_gain_calib$)
  x_angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) + &
                        (ele%value(crunch$) - ele%value(crunch_calib$))
  y_angle = (ele%value(tilt_tot$) - ele%value(tilt_calib$)) - &
                        (ele%value(crunch$) - ele%value(crunch_calib$))
  if (x_angle == 0 .and. y_angle == 0) then
    m_com%M_m(1,1) = x_gain
    m_com%M_m(1,2) = 0
    m_com%M_m(2,1) = 0
    m_com%M_m(2,2) = y_gain
  else
    m_com%M_m(1,1) =  x_gain * cos(x_angle)
    m_com%M_m(1,2) =  x_gain * sin(x_angle)
    m_com%M_m(2,1) = -y_gain * sin(y_angle)
    m_com%M_m(2,2) =  y_gain * cos(y_angle)
  endif
  m_com%x_off = ele%value(x_offset_tot$) + ele%value(x_offset_calib$)
  m_com%y_off = ele%value(y_offset_tot$) + ele%value(y_offset_calib$)
  m_com%value = ele%value
endif

end subroutine compute_bpm_transformation_numbers

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_orbit_reading (orb, ele, axis, reading, err)
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

subroutine to_orbit_reading (orb, ele, axis, reading, err)

use random_mod

implicit none

type (coord_struct) orb
type (ele_struct) ele

real(rp) reading, noise_factor
real(rp) ran_num, x, y

integer axis

character(20) :: r_name = "to_orbit_reading"

logical err

!

reading = 0.0
call check_if_ele_is_monitor (ele, err)
if (err) return

err = .true.

call compute_bpm_transformation_numbers (ele)

if (ele%value(noise$) /= 0) then
  call ran_gauss (ran_num)
  noise_factor = ele%value(noise$) * ran_num
else
  noise_factor = 0
endif

x = orb%vec(1) - m_com%x_off
y = orb%vec(3) - m_com%y_off

if (axis == x_plane$) then
  reading = noise_factor + (x * m_com%M_m(1,1) + y * m_com%M_m(1,2))
elseif (axis == y_plane$) then
  reading = noise_factor + (x * m_com%M_m(2,1) + y * m_com%M_m(2,2))
else
  reading = 0.0
  call out_io (s_warn$, r_name, "This axis not supported for BPM reading!")
  call err_exit
endif

err = .false.

end subroutine to_orbit_reading

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_eta_reading (eta, ele, axis, reading, err)
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
!  err -- Logical: Set True if there is an error. False otherwise.
!-

subroutine to_eta_reading (eta_actual, ele, axis, reading, err)

use random_mod

implicit none

type (ele_struct) ele

real(rp) eta_actual(:)
real(rp) reading, noise_factor
real(rp) ran_num, x, y

integer axis

character(20) :: r_name = "to_eta_reading"

logical err

!

reading = 0.0

call check_if_ele_is_monitor (ele, err)
if (err) return
if (.not. ele%is_on) return

!

err = .true.

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

call compute_bpm_transformation_numbers (ele)

x = eta_actual(1)
y = eta_actual(2)

if (axis == x_plane$) then
  reading = noise_factor + (x * m_com%M_m(1,1) + y * m_com%M_m(1,2))
elseif (axis == y_plane$) then
  reading = noise_factor + (x * m_com%M_m(2,1) + y * m_com%M_m(2,2))
else
  reading = 0.0
  call out_io (s_warn$, r_name, "This axis not supported for BPM reading!")
  call err_exit
endif

err = .false.
  
end subroutine to_eta_reading

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!+
! Subroutine to_phase_and_coupling_reading (ele, mon, err)
!
! Find the measured coupling values given the actual ones
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
!  mon -- bpm_phase_coupling_struct: K and Cbar coupling parameters
!  err -- Logical: Set True if there is an error. False otherwise.
!-

subroutine to_phase_and_coupling_reading (ele, mon, err)

use random_mod

implicit none

type (ele_struct) ele
type (bpm_phase_coupling_struct) mon, lab

real(rp) ran_num(6), cbar_lab(2,2), denom, ratio
real(rp) q_a1_lab(2), q_a2_lab(2), q_b1_lab(2), q_b2_lab(2)
real(rp) q_a1_mon(2), q_a2_mon(2), q_b1_mon(2), q_b2_mon(2)

integer axis

character(40) :: r_name = "to_phase_and_coupling_reading"

logical err

! 

mon%K_22a = 0
mon%K_12a = 0
mon%K_11b = 0
mon%K_12b = 0

call check_if_ele_is_monitor (ele, err)
if (err) return
if (.not. ele%is_on) return

!

err = .true.

call compute_bpm_transformation_numbers (ele)

if (ele%value(noise$) /= 0) then
  if (ele%value(n_sample$) == 0 .or. ele%value(osc_amplitude$) == 0) then
    call out_io (s_error$, r_name, 'N_SAMPLE OR OSC_AMPLITUDE NOT SET!')
    return
  endif
  call ran_gauss (ran_num)
  ran_num = ran_num * ele%value(noise$) / (ele%value(n_sample$) * ele%value(osc_amplitude$))
endif

!

call c_to_cbar (ele, cbar_lab)
ratio = m_com%sqrt_beta_b / m_com%sqrt_beta_a

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

q_a1_mon = matmul(m_com%M_m, q_a1_lab)
q_a2_mon = matmul(m_com%M_m, q_a2_lab)
q_b1_mon = matmul(m_com%M_m, q_b1_lab)
q_b2_mon = matmul(m_com%M_m, q_b2_lab)

! q_mon to k_mon

denom = q_a1_mon(1)**2 + q_a2_mon(1)**2
mon%K_22a = (q_a1_mon(1) * q_a1_mon(2) + q_a2_mon(1) * q_a2_mon(2)) / denom
mon%K_12a = (q_a1_mon(1) * q_a2_mon(2) - q_a2_mon(1) * q_a1_mon(2)) / denom

denom = q_b1_mon(2)**2 + q_b2_mon(2)**2
mon%K_11b = (q_b1_mon(2) * q_b1_mon(1) + q_b2_mon(2) * q_b2_mon(1)) / denom
mon%K_12b = (q_b1_mon(2) * q_b2_mon(1) - q_b2_mon(2) * q_b1_mon(1)) / denom

! Add random terms to k_mon

if (ele%value(noise$) /= 0) then
  mon%K_22a = mon%K_22a + ran_num(1)
  mon%K_12a = mon%K_12a + ran_num(2)
  mon%K_11b = mon%K_11b + ran_num(3)
  mon%K_12b = mon%K_12b + ran_num(4)
endif

! k_mon to Cbar_mon

mon%cbar22_a = -mon%K_22a * ele%gamma_c / ratio 
mon%cbar12_a = -mon%K_12a * ele%gamma_c / ratio 
mon%cbar11_b =  mon%K_11b * ele%gamma_c * ratio
mon%cbar12_b = -mon%K_12b * ele%gamma_c * ratio

! Phase calc

mon%phi_a = ele%a%phi - atan2(q_a2_mon(1), q_a1_mon(1)) + ran_num(5)
mon%phi_b = ele%b%phi - atan2(q_b2_mon(2), q_b1_mon(2)) + ran_num(6)

err = .false.

end subroutine to_phase_and_coupling_reading

end module
