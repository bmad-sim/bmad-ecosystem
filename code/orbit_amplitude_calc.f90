!+
! Subroutine orbit_amplitude_calc (ele, orb, amp_a, amp_b, &
!                                                  amp_na, amp_nb, particle)
!
! Subroutine to calculate the "invariant" amplitude of a particle at a 
! particular point in its orbit. This routine takes into account dispersion 
! and coupling. For uncoupled motion the horizontal motion is
!     x(s) = sqrt[2 * amp_a * beta_a(s)] * cos[phi(s)] + eta(s) * p_z(s)
! The amplitude is calculated at the exit end of ele.
! Notice that with the factor of 2 here the emittance is simply an average 
! of amp_a over all the particles in a bunch:
!     emit_a = <amp_a>
!
! The energy normalized amplitudes amp_na and amp_nb are used when the energy
! is changing. These are the true invariants. The relationship is:
!     amp_na = amp_a * gamma_E
!     amp_nb = amp_b * gamma_E
! where gamma_E is the usual energy gamma factor = Energy / mc^2.  
!
! Modules needed:
!   use bmad
!
! Input:
!   ele      -- Ele_struct: Element holding the Twiss parameters, 
!                  dispersion and coupling info.
!   orb      -- Coord_struct: Orbit coordinates at the exit end of ele.
!   particle -- Integer: Particle type (electron$, etc.).
!                 Needed for amp_na and amp_nb calc.
!
! Output:
!   amp_a  -- Real(rp), optional: a-mode amplitude
!   amp_b  -- Real(rp), optional: b-mode amplitude
!   amp_na -- Real(rp), optional: a-mode, energy normalized, amplitude.
!   amp_nb -- Real(rp), optional: b-mode, energy normalized, amplitude.
!-

subroutine orbit_amplitude_calc (ele, orb, amp_a, amp_b, &
                                              amp_na, amp_nb, particle)

  use bmad_struct
  use bmad_interface

  implicit none

  type (ele_struct) ele
  type (coord_struct) orb
  
  integer, optional :: particle
  real(rp), optional :: amp_a, amp_b, amp_na, amp_nb
  real(rp) v_mat(4,4), v_inv_mat(4,4), a_orb(4), amp

!

  call make_v_mats (ele, v_mat, v_inv_mat)
  a_orb = matmul (v_inv_mat, orb%vec(1:4)) - &
              orb%vec(6) * (/ ele%x%eta, ele%x%etap, ele%y%eta, ele%y%etap /)

  amp = (ele%x%gamma * a_orb(1)**2 + 2 * ele%x%alpha * a_orb(1)*a_orb(2) + &
                                                 ele%x%beta * a_orb(2)**2) / 2
  if (present(amp_a)) amp_a = amp
  if (present(amp_na)) amp_na = &
          amp * ele%value(beam_energy$) * (1 + orb%vec(6)) / mass_of(particle)

  amp = (ele%y%gamma * a_orb(3)**2 + 2 * ele%y%alpha * a_orb(3)*a_orb(4) + &
                                                 ele%y%beta * a_orb(4)**2) / 2
  if (present(amp_b)) amp_b = amp
  if (present(amp_nb)) amp_nb = &
          amp * ele%value(beam_energy$) * (1 + orb%vec(6)) / mass_of(particle)

end subroutine
