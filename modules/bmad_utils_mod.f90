#include "CESR_platform.inc"

module bmad_utils_mod

  use bmad_struct
  use bmad_interface

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine energy_to_kinetic (energy, particle, 
!                                           gamma, kinetic, beta, p0c, brho)
!
! Subroutine to calculate the kinetic energy, etc. from a particle's energy.
!
! Modules needed:
!   use bmad
!
! Input:
!   energy   -- Real(rp): Energy of the particle.
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   gamma   -- Real(rp), optional: Gamma factor.
!   kinetic -- Real(rp), optional: Kinetic energy
!   beta    -- Real(rp), optional: velocity / c_light
!   p0c     -- Real(rp), optional: Particle momentum
!   brho    -- Real(rp), optional: Nominal B_field*rho_bend
!-

subroutine energy_to_kinetic (energy, particle, &
                                            gamma, kinetic, beta, p0c, brho)

  implicit none

  real(rp), intent(in) :: energy
  real(rp), intent(out), optional :: kinetic, beta, p0c, brho, gamma
  real(rp) p0c_, mc2

  integer, intent(in) :: particle

!

  if (particle == positron$ .or. particle == electron$) then
    mc2 = m_electron
  else
    mc2 = m_proton
  endif

  p0c_ = sqrt(energy**2 - mc2**2)
  if (present(p0c))     p0c     = sqrt(energy**2 - mc2**2)
  if (present(beta))    beta    = p0c_ / energy  
  if (present(kinetic)) kinetic = energy - mc2
  if (present(brho))    brho    = p0c_ / c_light
  if (present(gamma))   gamma   = energy / mc2

end subroutine

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine wiggler_vec_potential (ele, energy, here, vec_pot)
!
! Subroutine to calculate the normalized vector potential at 
! a point for a wiggler. The normalized potental a_norm is defined by:
!      p_cononical = p_mv - a_norm
! The Gauge used here is the same one as used in PTC and has A_x = 0.
! 
! Modules needed:
!   use bmad
!
! Input:
!   ele     -- Ele_struct: wiggler element.
!   energy  -- Real(rdef): Particle energy.
!   here    -- Coord_struct: Coordinates for calculating the vector pot.
!
! Output:
!   vec_pot(3) -- Real(rdef): Normalized vector potential
!-

subroutine wiggler_vec_potential (ele, energy, here, vec_pot)

  implicit none

  type (ele_struct), target, intent(in) :: ele
  type (coord_struct), intent(in) :: here
  real(rdef), intent(in) :: energy
  real(rdef), intent(out) :: vec_pot(3)

  type (wig_term_struct), pointer :: t

  real(rdef) c_x, s_x, c_y, s_y, c_z, s_z
  real(rdef) x, y, s, coef

  integer i

!

  if (ele%key /= wiggler$) then
    print *, 'ERROR IN WIGGLER_VEC_POTENTIAL. ELEMENT NOT A WIGGLER: ', &
                                                                 ele%name
    call err_exit
  endif

!

  x = here%x%pos
  y = here%y%pos
  s = here%z%pos

  vec_pot = 0

  do i = 1, size(ele%wig_term)
    t => ele%wig_term(i)

      if (t%type == hyper_y$) then
        c_x = cos(t%kx * x)
        s_x = sin(t%kx * x)
      elseif (t%type == hyper_x$ .or. t%type == hyper_xy$) then
        c_x = cosh(t%kx * x)
        s_x = sinh(t%kx * x)
      else
        print *, 'ERROR IN WIGGLER_VEC_POTENTIAL: UNKNOWN TERM TYPE!'
        call err_exit
      endif

      if (t%type == hyper_y$ .or. t%type == hyper_xy$) then
        c_y = cosh (t%ky * y)
        s_y = sinh (t%ky * y)
      else
        c_y = cos (t%ky * y)
        s_y = sin (t%ky * y)
      endif

      c_z = cos (t%kz * s + t%phi_z)
      s_z = sin (t%kz * s + t%phi_z)

      coef = ele%value(polarity$) * t%coef

      vec_pot(2) = vec_pot(2) - coef  * (t%kz / (t%kx * t%ky)) * s_x * s_y * s_z
      vec_pot(3) = vec_pot(3) - coef  * (1 / t%kx)             * s_x * c_y * c_z
    enddo


end subroutine

end module
