!+
! Function spin_omega (field, orbit, sign_z_vel, phase_space_coords), result (omega)
!
! Return the modified T-BMT spin omega vector:
!   dOmega/d|s|   With phase space coords.
!   dOmega/dt     With time_RK coords.
!
! With phase_space_coords, d|s| is positive in the longitudinal direction of propagation 
! independent of the sign of sign_z_vel.
!
! In a bend, the omega returned should be modified:
!   true_omega = (1 + g*x) * omega + [0, g, 0]
!
! Input:
!   field              -- em_field_struct: E and B fields.
!   orbit              -- coord_struct: Particle position in element body coords.
!   sign_z_vel         -- integer: Direction of the z-velocity wrt the element body coords. 
!                           With phase space coords: sign_z_vel = orbit%direciton * ele%orientation.
!                           With time RK coords:     sign_z_vel is not used (orbit%vec(6) has the correct sign.)
!   phase_space_coords -- logical, optional: Is coord in standard phase space coordinates or
!                           is it time Runge Kutta coords? Default is True.
!
! Output:
!   omega(3)   -- real(rp): If phase_space_coords: Omega_TBMT/v_z
!                           If not: Omega_TBMT
!-

function spin_omega (field, orbit, sign_z_vel, phase_space_coords) result (omega)

use equal_mod, dummy_except => spin_omega

implicit none

type (em_field_struct) :: field
type (coord_struct) :: orbit

real(rp) omega(3),  beta_vec(3), pc
real(rp) anomalous_moment, gamma, rel_p, e_particle, mc2, bz2

integer sign_z_vel

logical, optional :: phase_space_coords

! Want everything in units of eV

anomalous_moment = anomalous_moment_of(orbit%species)
mc2 = mass_of(orbit%species)

if (logic_option(.true., phase_space_coords)) then
  rel_p = 1 + orbit%vec(6)
  e_particle = orbit%p0c * rel_p / orbit%beta
  gamma = e_particle / mc2
  bz2 = rel_p**2 - orbit%vec(2)**2 - orbit%vec(4)**2
  if (bz2 < 0) then  ! Particle has unphysical velocity
    omega = 0
    return
  endif
  beta_vec = (orbit%beta / rel_p) * [orbit%vec(2), orbit%vec(4), sign_z_vel * sqrt(bz2)]

! Can happen in an e_gun with the particle starting from rest.
elseif (orbit%beta == 0 .and. orbit%vec(2) == 0 .and. orbit%vec(4) == 0 .and. orbit%vec(6) == 0) then
  e_particle = mc2
  beta_vec = 0
  gamma = 1

else
  e_particle = sqrt(orbit%vec(2)**2 + orbit%vec(4)**2 + orbit%vec(6)**2) / orbit%beta
  beta_vec = [orbit%vec(2), orbit%vec(4), orbit%vec(6)] / e_particle
  gamma = e_particle / mc2
endif

omega = c_light * (1/gamma + anomalous_moment) * field%B
omega = omega - c_light * (gamma * anomalous_moment * dot_product(beta_vec, field%B) / (gamma + 1)) * beta_vec
omega = omega - (anomalous_moment + 1/(1+gamma)) * cross_product(beta_vec, field%E)

if (bmad_com%electric_dipole_moment /= 0) then
  omega = omega + (bmad_com%electric_dipole_moment / 2) * &
            (field%E - (gamma * dot_product(beta_vec, field%E)/ (1 + gamma)) * beta_vec + &
             c_light * cross_product(beta_vec, field%B))
endif

if (logic_option(.true., phase_space_coords)) then
  omega = -(charge_of(orbit%species) / (mc2 * abs(beta_vec(3)))) * omega
else
  omega = -(charge_of(orbit%species) * c_light / mc2) * omega
endif

end function spin_omega


