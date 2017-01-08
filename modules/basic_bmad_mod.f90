!+
! Module basic_bmad_mod
!
! Some basic routines independent of any bmad structures.
!-

module basic_bmad_mod

use sim_utils

implicit none

integer, parameter :: n_pole_maxx = 21  ! maximum multipole order

integer, parameter :: old_control_var_offset$ = 1000  ! For indexing into ele%control_var(:) array
integer, parameter :: var_offset$ = 2000              ! Important: var_offset$ > old_control_var_offset$
integer, parameter :: taylor_offset$ = 1000000000     ! Taylor term index offset.

type expression_atom_struct
  character(40) :: name = ''
  integer :: type = 0
  real(rp) :: value = 0
end type

type controller_var_struct
  character(40) :: name = ''
  real(rp) :: value = 0
  real(rp) :: old_value = 0
end type


contains

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+ 
! Subroutine reallocate_expression_stack (stack, n, exact)
!
! Routine to resize an expression stack array.
!
! Input:
!   stack(:)  -- expression_atom_stuct, allocatable: Existing stack array.
!   n         -- integer: Array size needed.
!   exact     -- logical, optional: If present and False then the size of the
!                  output array is permitted to be larger than n.
!                  Default is True.
! Output:
!   stack(:)  -- expression_atom_struct, allocatable: Resized stack.
!                  The stack info will be preserved.
!-

subroutine reallocate_expression_stack (stack, n, exact)

type (expression_atom_struct), allocatable :: stack(:)
type (expression_atom_struct), allocatable :: stack_temp(:)
integer n, nn
logical, optional :: exact

!

if (.not. allocated(stack)) allocate(stack(n))

if (size(stack) == n .or. size(stack) > n .and. .not. logic_option(.true., exact)) return

call move_alloc (stack, stack_temp)
allocate (stack(n))
nn = min(size(stack_temp), size(stack))
stack(1:nn) = stack_temp(1:nn)

end subroutine reallocate_expression_stack

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
  call out_io (s_abort$, r_name, 'ERROR: TOTAL ENERGY IS LESS THAN REST MASS: \es10.2\ ', E_tot)
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

if (present(dbeta)) dbeta = -sqrt_one(-(mc2/E_tot)**2)

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
if (present(dbeta))   dbeta   = -sqrt_one(-(mc2/E_tot_this)**2)

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
