!+ 
! Subroutine energy_to_kinetic (energy, particle, 
!                                           gamma, kinetic, beta, p0c, brho)
!
! Subroutine to calculate the kinetic energy, etc. from a particle's energy.
!
! Modules needed:
!   use accelerator
!
! Input:
!   energy   -- Real*8: Energy of the particle.
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   gamma   -- Real*8, optional: Gamma factor.
!   kinetic -- Real*8, optional: Kinetic energy
!   beta    -- Real*8, optional: velocity / c_light
!   p0c     -- Real*8, optional: Particle momentum
!   brho    -- Real*8, optional: Nominal B_field*rho_bend
!-

subroutine energy_to_kinetic (energy, particle, gamma, kinetic, beta, p0c, brho)

  use accelerator

  implicit none

  real*8, intent(in) :: energy
  real*8, intent(out), optional :: kinetic, beta, p0c, brho, gamma
  real*8 p0c_, mc2

  integer, intent(in) :: particle

!

  if (particle == positron$ .or. particle == electron$) then
    mc2 = e_mass
  else
    mc2 = p_mass
  endif

  p0c_    = sqrt(energy**2 - mc2**2)
  if (present(p0c))     p0c     = sqrt(energy**2 - mc2**2)
  if (present(beta))    beta    = p0c_ / energy  
  if (present(kinetic)) kinetic = energy - mc2
  if (present(brho))    brho    = p0c_ * 1d9 / c_light
  if (present(gamma))   gamma   = energy / mc2

end subroutine
