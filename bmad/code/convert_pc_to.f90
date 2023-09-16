!+ 
! Subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, beta1, err_flag)
!
! Routine to calculate the energy, etc. from a particle's momentum.
!
! Input:
!   pc       -- Real(rp): Particle momentum
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   E_tot    -- Real(rp), optional: Total energy of the particle.
!   gamma    -- Real(rp), optional: Gamma factor.
!   kinetic  -- Real(rp), optional: Kinetic energy
!   beta     -- Real(rp), optional: velocity / c_light
!   brho     -- Real(rp), optional: Nominal B_field*rho_bend
!   beta1    -- Real(rp), optional: 1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, beta1, err_flag)

use sim_utils

implicit none

real(rp), intent(in) :: pc
real(rp), intent(out), optional :: E_tot, kinetic, beta, brho, gamma, beta1
real(rp) g2, mc2, E_tot_this 

integer, intent(in) :: particle
logical, optional :: err_flag

character(20) :: r_name = 'convert_pc_to'

!

if (present(err_flag)) err_flag = .false.

mc2 = mass_of(particle)
E_tot_this = sqrt(pc**2 + mc2**2)

if (present(E_tot))   E_tot   = E_tot_this
if (present(beta))    beta    = pc / E_tot_this
if (present(kinetic)) kinetic = E_tot_this - mc2
if (present(brho))    brho    = pc / c_light
if (present(gamma))   gamma   = E_tot_this / mc2
if (present(beta1))   beta1   = -sqrt_one(-(mc2/E_tot_this)**2)

if (present(err_flag)) err_flag = .false.

end subroutine convert_pc_to
