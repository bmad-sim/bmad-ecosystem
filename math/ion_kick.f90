!+
! subroutine ion_kick (x, y, x_kicker, y_kicker, s_kicker)
!
! subroutine to return the kick felt by an ion due to the
! passage of a bunch.
!
! This subroutine uses MKS units: the output kicks are in m/sec
!
! n_elec     -- the number of electrons in a bunch
! ion_weight  -- the weight of the ion in AMU
! sig_x, sig_y -- rms widths of the electron bunch
! alpha, beta  -- TWISS parameters
!-

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

end subroutine
