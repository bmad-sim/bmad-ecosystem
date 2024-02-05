!+
! Subroutine ion_kick (orbit, r_beam, n_beam_part, a_twiss, b_twiss, sig_ee, kick)
!
! Routine to return the kick felt by an ion due to the passage of a bunch.
!
! Note: If a_twiss%sigma or b_twiss%sigma are set then this routine uses these
! numbers. Otherwise, it computes the beam size based on the Twiss parameters and emittance.
!
! Input:
!   orbit       -- coord_struct: Ion position.
!   r_beam(2)   -- real(rp): Beam (x, y) position.
!   n_beam_part -- real(rp): Number of beam particles.
!   a_twiss     -- twiss_struct: Horizontal like beam twiss parameters.
!   b_twiss     -- twiss_struct: vertical like beam twiss parameters.
!   sig_ee      -- real(rp): Sigma_E/E beam energy spread.
!
! Output:
!   kick(3)     -- real(rp): (x, y, s) kick in m/sec.
!-

subroutine ion_kick (orbit, r_beam, n_beam_part, a_twiss, b_twiss, sig_ee, kick)

use bmad_interface, except_dummy => ion_kick

implicit none

type (coord_struct) orbit
type (twiss_struct) a_twiss, b_twiss

real(rp) r_beam(2), n_beam_part, sig_ee, kick(3)
real(rp) scale, k(2), dk(2,2), x, y, ion_weight, sig_x, sig_y

real(rp), parameter :: factor = r_p * c_light / twopi  ! 7.3e-11

character(*), parameter :: r_name = 'ion_kick'

!

if ((a_twiss%sigma == 0 .and. a_twiss%emit == 0) .or. (b_twiss%sigma == 0 .and. b_twiss%emit == 0)) then
  call out_io (s_fatal$, r_name, 'BOTH SIGMA AND EMITTANCE ARE ZERO FOR BEAM.')
  if (global_com%exit_on_error) call err_exit
  return
endif

ion_weight = mass_of(orbit%species)

sig_x = a_twiss%sigma
if (sig_x == 0) sig_x = sqrt(a_twiss%emit * a_twiss%beta + (a_twiss%eta*sig_ee)**2)

sig_y = a_twiss%sigma
if (sig_y == 0) sig_y = sqrt(b_twiss%emit * b_twiss%beta + (b_twiss%eta*sig_ee)**2)

scale = factor * n_beam_part  / ((sig_x + sig_y) * ion_weight)
x = orbit%vec(1) - r_beam(1)
y = orbit%vec(3) - r_beam(2)

call bbi_kick (x, y, [sig_x, sig_y], k, dk)

kick(1) = scale * k(1)
kick(2) = scale * k(2)
kick(3) = scale * (-a_twiss%alpha * a_twiss%emit + a_twiss%eta * a_twiss%etap * sig_ee**2) * dk(1,1) + &
                  (-b_twiss%alpha * b_twiss%emit + b_twiss%eta * b_twiss%etap * sig_ee**2) * dk(2,2)

end subroutine
