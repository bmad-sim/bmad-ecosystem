!+ 
! Subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, beta1, err_flag)
!
! Routine to calculate the momentum, etc. from a particle's total energy.
!
! Input:
!   E_tot    -- Real(rp): Total energy of the particle.
!   particle -- Integer: Type of particle. positron$, etc.
!
! Output:
!   gamma    -- Real(rp), optional: Gamma factor. Set to -1 for photons.
!   kinetic  -- Real(rp), optional: Kinetic energy
!   beta     -- Real(rp), optional: velocity / c_light
!   pc       -- Real(rp), optional: Particle momentum
!   brho     -- Real(rp), optional: Nominal B_field*rho_bend
!   beta1    -- Real(rp), optional: 1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, beta1, err_flag)

use sim_utils

implicit none

real(rp), intent(in) :: E_tot
real(rp), intent(out), optional :: kinetic, beta, pc, brho, gamma, beta1
real(rp) pc_new, mc2, g2

integer, intent(in) :: particle
logical, optional :: err_flag

character(24) :: r_name = 'convert_total_energy_to'

!

if (present(err_flag)) err_flag = .true.

mc2 = mass_of(particle)
if (E_tot < mc2 .or. E_tot == 0) then
  call out_io (s_abort$, r_name, 'ERROR: TOTAL ENERGY IS LESS THAN REST MASS: \es10.2\ ', E_tot)
  if (global_com%exit_on_error) call err_exit
  return
endif

pc_new = sqrt(E_tot**2 - mc2**2)
if (present(pc))     pc       = pc_new
if (present(beta))    beta    = pc_new / E_tot  
if (present(kinetic)) kinetic = E_tot - mc2
if (present(brho))    brho    = pc_new / c_light

if (present(gamma)) then
  if (mc2 == 0) then
    gamma = -1
  else
    gamma   = E_tot / mc2
  endif
endif

if (present(beta1)) beta1 = -sqrt_one(-(mc2/E_tot)**2)

if (present(err_flag)) err_flag = .false.

end subroutine convert_total_energy_to
