!+
! subroutine ion_kick_2d (x, y, x_kicker, y_kicker)
!
! subroutine to return the transverse kick felt by a ion
! due to the passage of an electron bunch.
!
! This subroutine uses MKS units.
! The output kicks are in m/sec
!
! n_elec                -- the number of electrons in a bunch
! ion_mass_charge_ratio -- the weight / mass ratio of the ion in kg / Coul
! sig_x, sig_y          -- rms widths of the electron bunch
!-

subroutine ion_kick_2d (x, y, x_kicker, y_kicker)

  use precision_def

  implicit none

  real(rp) factor, scale, kx, ky, x_position, y_position
  real(rp) x, y, r, x_kicker, y_kicker

  real(rp) n_elec, ion_mass_charge_ratio, sig_x, sig_y

  common / ratio / r
  common / ion_kick_2d_com / n_elec, ion_mass_charge_ratio, sig_x, sig_y


! factor = r_proton * c_light / (2 * pi * N_avrogado * e_charge)

  parameter (factor = 1.534e-18 * 2.998e8 / &
                           (2 * 3.1415926 * 6.02e26 * 1.6e-19))

!

  r = sig_y/sig_x
  scale = factor * n_elec / ((sig_x + sig_y) * ion_mass_charge_ratio)
  x_position = x / sig_x
  y_position = y / sig_y

  call bbi_kick (x_position, y_position, kx, ky)

  x_kicker = scale * kx
  y_kicker = scale * ky

end subroutine
