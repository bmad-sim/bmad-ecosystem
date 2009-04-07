!  ION_KICK    SUBROUTINE  BEAM-ION    C.DCS.LIB   DCS         96.7
!+
!     subroutine ION_KICK
!
!     subroutine to return the kick felt by an ion due to the
!     passage of a bunch.
!
!     This subroutine uses MKS units: the output kicks are in m/sec
!
!     n_elec     -- the number of electrons in a bunch
!     ion_weight  -- the weight of the ion in AMU
!     sig_x, sig_y -- rms widths of the electron bunch
!     alpha, beta  -- TWISS parameters
!-

!$Id$
!$Log$
!Revision 1.5  2003/07/09 01:29:30  dcs
!new bmad
!
!Revision 1.4  2002/02/23 20:34:45  dcs
!Modified for Single/Double Real Toggle
!
!Revision 1.3  2001/10/25 19:13:15  helms
!Added explicit variable declarations.
!mat_inv has been replaced by mat_inverse in BMAD (not DCSLIB)
!
!Revision 1.2  2001/09/27 17:47:06  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


subroutine ion_kick (x, y, x_kicker, y_kicker, s_kicker)

  use precision_def

  implicit none

  real(rp) r, factor, scale, kx, ky, kxx, kyy
  real(rp) x_position, y_position
  real(rp) x, y, x_kicker, y_kicker, s_kicker

  real(rp) n_elec, ion_weight, sig_x, sig_y
  real(rp) alpha_x, alpha_y, epsilon_x, epsilon_y
  real(rp) sig_ee, eta, etap

  common /ratio/ r
  common /lattice/ alpha_x, alpha_y, eta, etap, sig_x, sig_y

  common /kick_com/ n_elec, ion_weight, epsilon_x, epsilon_y, sig_ee

  parameter (factor = 7.319e-11)    ! = r_proton * v_light / (2 * pi)

!

  r = sig_y/sig_x
  scale = factor * n_elec  / ((sig_x + sig_y) * ion_weight)
  x_position = x / sig_x
  y_position = y / sig_y

  call bbi_kick (x_position + .01, y_position, kxx, ky)
  call bbi_kick (x_position, y_position + .01, kx, kyy)
  call bbi_kick (x_position, y_position, kx, ky)

  x_kicker = scale * kx
  y_kicker = scale * ky
  s_kicker = scale * ((-alpha_x*epsilon_x + eta*etap*sig_ee*sig_ee)* &
      (kxx-kx)/(.01*sig_x) - alpha_y*epsilon_y*(kyy-ky)/(.01*sig_y))

  return

  end
