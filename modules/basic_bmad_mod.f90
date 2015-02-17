!+
! Module basic_bmad_mod
!
! Some basic routines independent of any bmad structures.
!-

module basic_bmad_mod

use sim_utils

integer, parameter :: n_pole_maxx = 21  ! maximum multipole order

! Species ID, mass, and charge.

integer, parameter :: not_set$ = -999

integer, parameter :: ref_particle$ = 6, anti_ref_particle$ = 7
integer, parameter :: pion_0$     = +5
integer, parameter :: pion_plus$  = +4
integer, parameter :: antimuon$   = +3
integer, parameter :: proton$     = +2
integer, parameter :: positron$   = +1
integer, parameter :: photon$     =  0
integer, parameter :: electron$   = -1
integer, parameter :: antiproton$ = -2
integer, parameter :: muon$       = -3
integer, parameter :: pion_minus$ = -4

character(20), parameter :: particle_name(-4:7) = [&
                'Pion_Minus       ', 'Muon             ', 'Antiproton       ', 'Electron         ', &
                'Photon           ', 'Positron         ', 'Proton           ', 'Antimuon         ', &
                'Pion_Plus        ', 'Pion_0           ', 'Ref_Particle     ', 'Anti_Ref_Particle']

integer, parameter :: charge_of(-4:5) = [-1, -1, -1, -1, 0, 1, 1, 1, 1, 0]
real(rp), parameter :: mass_of(-4:5) = [m_pion_charged, m_muon, m_proton, m_electron, 0.0_rp, &
                                m_electron, m_proton, m_muon, m_pion_charged, m_pion_0]

real(rp), parameter :: anomalous_moment_of(-4:5) = [0.0_rp, anomalous_mag_moment_muon, &
                        anomalous_mag_moment_proton, anomalous_mag_moment_electron, 0.0_rp, &
                        anomalous_mag_moment_electron, anomalous_mag_moment_proton, &
                        anomalous_mag_moment_muon, 0.0_rp, 0.0_rp]

integer, parameter :: antiparticle(-4:7) = [4, 3, 2, 1, 0, -1, -2, -3, -4, 5, 7, 6]

contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function is_true (param) result (this_true)
!
! Routine to translate from a real number to a boolian True or False.
! Translation: 0 = False, nonzero = True
!
! The typical use of this routine is for parameters in ele_struct%value(:) which
! is a real array. Some of the elements in the %value array are used to specify
! boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).
! 
! Input:
!   param    -- real(rp): Real number to be translated
!
! Output:
!   this_true -- logical: Set False if param is zero. True otherwise.
!-

function is_true (param) result (this_true)

implicit none

real(rp) param
logical this_true

!

this_true = (param /= 0)

end function is_true

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function is_false (param) result (this_false)
!
! Routine to translate from a real number to a boolian True or False.
! Translation: 0 = False, nonzero = True
!
! The typical use of this routine is for parameters in ele_struct%value(:) which
! is a real array. Some of the elements in the %value array are used to specify
! boolian attributes. For example, quadrupoles use ele%value(scale_multipoles$).
! 
! Input:
!   param    -- real(rp): Real number to be translated
!
! Output:
!   this_false -- logical: Set True if param is zero. False otherwise.
!-

function is_false (param) result (this_false)

implicit none

real(rp) param
logical this_false

!

this_false = (param == 0)

end function is_false

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Function species_name (species) result (species_str)
!
! Routine to return the string representation of the type of particle under consideration.
! This routine is similar in action to the particle_name array except that
! this routine does not have a problem with array index out of bounds.
!
! Module needed:
!   use bmad
!
! Input:
!   species     -- integer: Species. positron$, etc.
!
! Output:
!   species_str -- Character(12): String representation.
!-

Function species_name (species) result (species_str)

implicit none

integer species
character(12) species_str

!

if (species < lbound(particle_name, 1) .or. species > ubound(particle_name, 1)) then
  select case (species)
  case (not_set$); species_str = 'Not_Set!'
  case default;    species_str = 'UNKNOWN!'
  end select
  return
endif

species_str = particle_name(species)

end function species_name

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, dbeta, err_flag)
!
! Routine to calculate the momentum, etc. from a particle's total energy.
!
! Modules needed:
!   use bmad
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
!   dbeta    -- Real(rp), optional: 1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine convert_total_energy_to (E_tot, particle, gamma, kinetic, beta, pc, brho, dbeta, err_flag)

implicit none

real(rp), intent(in) :: E_tot
real(rp), intent(out), optional :: kinetic, beta, pc, brho, gamma, dbeta
real(rp) pc_new, mc2, g2

integer, intent(in) :: particle
logical, optional :: err_flag

character(24) :: r_name = 'convert_total_energy_to'

!

if (present(err_flag)) err_flag = .true.

mc2 = mass_of(particle)
if (E_tot < mc2 .or. E_tot == 0) then
  call out_io (s_abort$, r_name, 'ERROR: TOTAL ENERGY IS LESS THAN REST MASS:\f10.0\ ', E_tot)
  if (global_com%exit_on_error) call err_exit
  return
endif

pc_new = E_tot * sqrt(1.0 - (mc2/E_tot)**2)
if (present(pc))     pc     = pc_new
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

if (present(dbeta)) then
  if (E_tot/mc2 > 100) then
    g2 = (E_tot / mc2)**2
    dbeta = 1/(2*g2) + 1/(8*g2**2)
  else
    dbeta = 1 - pc_new / E_tot
  endif
endif

if (present(err_flag)) err_flag = .false.

end subroutine convert_total_energy_to

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, dbeta, err_flag)
!
! Routine to calculate the energy, etc. from a particle's momentum.
!
! Modules needed:
!   use bmad
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
!   dbeta    -- Real(rp), optional: 1 - beta. Equal to 1/(2*gamma^2) in ultra-rel limit.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine convert_pc_to (pc, particle, E_tot, gamma, kinetic, beta, brho, dbeta, err_flag)

implicit none

real(rp), intent(in) :: pc
real(rp), intent(out), optional :: E_tot, kinetic, beta, brho, gamma, dbeta
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

if (present(dbeta)) then
  if (E_tot/mc2 > 100) then
    g2 = (E_tot_this / mc2)**2
    dbeta = 1/(2*g2) + 1/(8*g2**2)
  else
    dbeta = 1 - pc / E_tot_this
  endif
endif

if (present(err_flag)) err_flag = .false.

end subroutine convert_pc_to

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Subroutine compute_even_steps (ds_in, length, ds_default, ds_out, n_step)
!
! Subroutine to compute a step size ds_out, close to ds_in, so that an 
! integer number of steps spans the length:
!   length = ds_out * n_step
!
! Modules needed:
!   use bmad
!
! Input:
!   ds_in      -- Real(rp): Input step size.
!   length     -- Real(rp): Total length.
!   ds_default -- Real(rp): Default to use if ds_in = 0.
!
! Output:
!   ds_out    -- Real(rp): Step size to use.
!   n_step    -- Integer: Number of steps needed.
!-

subroutine compute_even_steps (ds_in, length, ds_default, ds_out, n_step)

implicit none

real(rp) ds_in, length, ds_default, ds_out
integer n_step

!

ds_out = ds_in
if (ds_out == 0) ds_out = ds_default
n_step = nint(length / ds_out)
if (n_step == 0) n_step = 1
ds_out = length / n_step  

end subroutine compute_even_steps

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
! Function c_multi (n, m, no_n_fact, c_full) result (c_out)
!
! Subroutine to compute multipole factors:
!          c_multi(n, m) =  +/- ("n choose m")/n!
! This is used in calculating multipoles.
!
! Input:
!   n,m       -- Integer: For n choose m
!   no_n_fact -- Logical, optional: If present and true then
!                 c_out = +/- "n choose m".
!   c_full(:,:) --    real(rp), optional:  If present, will be populated with
!                     return entire c(n_pole_maxx,n_pole_maxx) matrix
!
! Output:
!   c_out  -- Real(rp): Multipole factor.
!-

function c_multi (n, m, no_n_fact, c_full) result (c_out)

implicit none

integer, intent(in) :: n, m
integer in, im

real(rp) c_out
real(rp), save :: n_factorial(0:n_pole_maxx)
real(rp), save :: c(0:n_pole_maxx, 0:n_pole_maxx)
real(rp), optional :: c_full(0:n_pole_maxx, 0:n_pole_maxx)

logical, save :: init_needed = .true.
logical, optional :: no_n_fact

! The magnitude of c(n, m) is number of combinations normalized by n!

if (init_needed) then

  c(0, 0) = 1

  do in = 1, n_pole_maxx
    c(in, 0) = 1
    c(in, in) = 1
    do im = 1, in-1
      c(in, im) = c(in-1, im-1) + c(in-1, im)
    enddo
  enddo

  n_factorial(0) = 1

  do in = 0, n_pole_maxx
    if (in > 0) n_factorial(in) = in * n_factorial(in-1)
    do im = 0, in
      c(in, im) = c(in, im) / n_factorial(in)
      if (mod(im, 4) == 0) c(in, im) = -c(in, im)
      if (mod(im, 4) == 3) c(in, im) = -c(in, im)
    enddo
  enddo

  init_needed = .false.

endif

!

if (logic_option (.false., no_n_fact)) then
  c_out = c(n, m) * n_factorial(n)
else
  c_out = c(n, m)
endif

if (present(c_full)) c_full = c

end function c_multi

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function mexp (x, m) result (this_exp)
!
! Returns x^m with 0^0 = 1.
!
! Modules needed:
!   use bmad
!
! Input:
!   x -- Real(rp): Number.
!   m -- Integer: Exponent.
!
! Output:
!   this_exp -- Real(rp): Result.
!-

function mexp (x, m) result (this_exp)

implicit none

real(rp) x, this_exp
integer m

!

if (m < 0) then
  this_exp = 0
elseif (m == 0) then
  this_exp = 1
else
  this_exp = x**m
endif

end function mexp

end module
