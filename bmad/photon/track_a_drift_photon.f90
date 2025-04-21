!+
! Subroutine track_a_drift_photon (orb, length, phase_relative_to_ref)
!
! Subroutine to track a particle as through a drift.
!
! Input:
!   orb      -- coord_struct: Orbit at start of the drift.
!   length   -- real(rp): Longitudinal length to drift through.
!   phase_relative_to_ref
!            -- logical: If true then E field phase shift is relative to ref particle.
!
! Output:
!   orb      -- coord_struct: Orbit at end of the drift
!-

subroutine track_a_drift_photon (orb, length, phase_relative_to_ref)

use bmad_routine_interface, except_dummy => track_a_drift_photon

implicit none

type (coord_struct) orb
real(rp) length, path_len, l, v2, dph
logical phase_relative_to_ref

! Check for lost

if (orb%vec(6) == 0) then
  orb%state = lost$
  return
endif

! Photon tracking uses a different coordinate system. 
! Notice that if orb%vec(6) is negative then the photon will be going back in time.

l = length  ! In case actual length argument is a component of orb.
path_len = l / orb%vec(6)
orb%vec(1) = orb%vec(1) + path_len * orb%vec(2)
orb%vec(3) = orb%vec(3) + path_len * orb%vec(4)
orb%vec(5) = orb%vec(5) + l
orb%s      = orb%s      + l
orb%t      = orb%t + path_len / c_light
orb%dt_ref = orb%dt_ref + path_len / c_light

if (phase_relative_to_ref) then
  v2 = orb%vec(2)**2 + orb%vec(4)**2
  if (v2 < 1d-4) then
    dph = abs(l) * v2 * (0.5_rp + 3 * v2 / 8 + 5 * v2**2 / 16)
  else
    dph = abs(path_len) - abs(l)
  endif
  orb%phase = orb%phase + sign(dph, path_len) * orb%p0c / (c_light * h_bar_planck)
else
  orb%phase = orb%phase + path_len * orb%p0c / (c_light * h_bar_planck)
endif
end subroutine track_a_drift_photon

