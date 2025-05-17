!+
! Subroutine tao_de_optimizer (abort)
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches the measured data.
! 
! This subroutine is a wrapper for the opti_de optimizer in sim_utils.
! 'de' stands for 'differential evolution' see opti_de routine for
! more details.
!
! Output:
!   abort -- Logical: Set True if an user stop signal detected.
!-

subroutine tao_de_optimizer (abort)

use tao_interface, dummy => tao_de_optimizer
use tao_top10_mod, only: tao_var_write
use opti_de_mod
!!! use opti_de_openmp_mod  # And there is commented out code below.
!!! use omp_lib, only: omp_get_max_threads

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable :: var_vec(:), var_step(:)
real(rp) merit_start, merit_end, merit

integer i, n, gen, pop, n_var, population, status

character(20) :: r_name = 'tao_de_optimizer'
character(80) line

logical abort

! setup

abort = .false.

! put the variable values into an array for the optimizer

call tao_get_opt_vars (var_vec, var_step = var_step)
var_step = var_step * s%global%de_lm_step_ratio
n_var = size(var_vec)

population=s%global%de_var_to_population_factor*n_var
population = max(population, 20)
merit_start = tao_merit ()

! run the optimizer

write (line, '(a, i0)') 'Differential evolution optimizer, population: ', population
call out_io (s_blank$, r_name, line)

!!! if (omp_get_max_threads() == 1) then
  merit = opti_de (var_vec, s%global%n_opti_cycles, population, merit_wrapper, var_step, status)
!!! else
!!!  merit = opti_de_openmp (var_vec, s%global%n_opti_cycles, population, merit_wrapper, var_step, status)
!!! endif

print *, 'tao_de_optimizer merit for rank ', merit, s%mpi%rank

! cleanup after the optimizer

call tao_set_opt_vars (var_vec, s%global%optimizer_var_limit_warn)
merit_end = tao_merit ()
if (s%global%opti_write_var_file) call tao_var_write (s%global%var_out_file)

write (line, '(a, es14.6)') 'Merit start:', merit_start
call out_io (s_blank$, r_name, line)
write (line, '(a, es14.6)') 'Merit end:  ', merit_end
call out_io (s_blank$, r_name, line)

if (status /= 0) abort = .true.

contains

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function merit_wrapper (var_vec, status, iter_count) result (this_merit)
!
! Function called by opti_de to set the variables and return the merit value.
!
! Input:
!   var_vec(:) -- Array of variables.
!   iter_count -- Integer: Number of times this routine has been called since
!                 the start of the optimization modulo 1000.
!
! Output:
!   status     -- Integer: Set True if we want opti_de to halt.
!   this_merit -- Real(rp): Value of the merit function.
!   iter_count -- Integer: Increased by one modulo 1000 from the input number.
!-

function merit_wrapper (var_vec, status, iter_count) result (this_merit)

use tao_interface
use tao_top10_mod, only: tao_var_write
use input_mod

implicit none

real(rp) var_vec(:)
real(rp) this_merit, merit
real(rp), save :: merit_min, merit_min_out, merit_min_type

integer i, status, rank
integer iter_count
integer, save :: t0(8), t1(8), t_del(8), t_delta

character(80) line, line2, stars
character(*), parameter :: r_name = 'tao_de_optimizer'
character(1) char

logical calc_ok

! Init

if (iter_count == 0) then
  merit_min = 1e35
  merit_min_out = 1e35
  merit_min_type = 1e35
endif

! look for keyboard input to end optimization

status = 0  ! continue
if (tao_user_is_terminating_optimization()) status = 1

!

stars = '****************************************************'

call tao_set_opt_vars (var_vec, s%global%optimizer_var_limit_warn)

this_merit = tao_merit (calc_ok)
merit_min = min(merit_min, this_merit)

if (iter_count == 1000) then
  call date_and_time (values = t1)
  t_del = t1 - t0
  t_delta = t_del(7) + 60*(t_del(6) + 60*(t_del(5) + 24*(t_del(3) + 30*t_del(2)))) 
endif

if (this_merit <= 0.98*merit_min_type .or. t_delta > 10) then
  write (line, '(a, es14.6)') ' So far the minimum is ', merit_min
  if (s%global%opti_write_var_file) call tao_var_write (s%global%var_out_file)

  if (calc_ok) then
    call out_io (s_blank$, r_name, stars, line, stars)
  else
    write (line2, *) 'Computation had problems...'
    call out_io (s_blank$, r_name, stars, line, line2, stars)
  endif

  call date_and_time (values = t0)
  t_delta = 0
  merit_min_type = merit_min
endif

if (this_merit < 1e-10 .and. s%com%all_merit_weights_positive) then
  call out_io (s_blank$, r_name, stars, ' MERIT < 1E-10 ==> AT MINIMUM. QUITING HERE.', stars)
  status = 1
endif

if (this_merit <= 0.9*merit_min_out) then
  merit_min_out = this_merit
endif

iter_count = mod(iter_count, 1000) + 1

end function merit_wrapper

end subroutine tao_de_optimizer
