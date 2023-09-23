                                    
!+
! Subroutine tao_lmdif_optimizer (abort)
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches the measured data.
! 
! This subroutine is a wrapper for the mrqmin routine of Numerical Recipes.
! See the Numerical Recipes writeup for more details.
! 'lm' stands for Levenburg - Marquardt. Otherwise known as LMDIF. 
!
! Output:
!   abort -- Logical: Set True if an user stop signal detected or there is 
!            a problem with calculating the merit function.
!-

subroutine tao_lmdif_optimizer (abort)

use lmdif_mod
use tao_interface, dummy => tao_lmdif_optimizer
use tao_top10_mod, only: tao_var_write

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable :: merit_vec(:), weight(:)
real(rp), allocatable :: var_delta(:), var_value(:), var_at_min(:)
real(rp) merit, merit_at_min

integer i, j, k, n
integer n_data, n_var

logical :: abort, init_needed = .true.
logical at_end

character(20) :: r_name = 'tao_lmdif_optimizer'
character(80) :: line
character(1) char

! setup

abort = .false.
merit = tao_merit()

merit_at_min = merit

call tao_get_opt_vars (var_value, var_delta = var_delta, var_weight = weight)

n_var = size(var_delta)
n_data = n_var
do i = lbound(s%u, 1), ubound(s%u, 1)
  if (.not. s%u(i)%is_on) cycle
  n_data = n_data + count(s%u(i)%data(:)%useit_opt .and. s%u(i)%data(:)%weight /= 0)
enddo

call re_allocate (merit_vec, n_data)
call re_allocate (var_at_min, n_var)
var_at_min = var_value

if (s%global%n_opti_cycles <= 1.2*n_var+5) then
  call out_io (s_warn$, r_name, &
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
            '!!!!! FOR THE LMDIF OPTIMIZER TO WORK WELL, N_OPTI_CYCLES (\i0\) !!!!!', &
            '!!!!! MUST BE WELL ABOVE THE NUMBER OF VARIABLES (\i0\)          !!!!!', &
            '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', &
            i_array = (/ s%global%n_opti_cycles, n_var /))
endif

! run optimizer 

call initial_lmdif

call out_io (s_blank$, r_name, '  Cycle      Merit')

cycle_loop: do i = 1, s%global%n_opti_cycles

  merit_vec = 0
  merit_vec(1:n_var) = sqrt(weight) * var_delta
  k = n_var
  do n = lbound(s%u, 1), ubound(s%u, 1)
    u => s%u(n)
    if (.not. u%is_on) cycle
    do j = 1, size(u%data)
      if (.not. u%data(j)%useit_opt) cycle
      if (u%data(j)%weight == 0) cycle
      k = k + 1
      merit_vec(k) = sqrt(u%data(j)%weight) * u%data(j)%delta_merit
    enddo
  enddo

  call suggest_lmdif (var_value, merit_vec, s%global%lmdif_eps, s%global%n_opti_cycles, at_end)
  call tao_set_opt_vars (var_value, s%global%optimizer_var_limit_warn)
  merit = tao_merit()
  if (merit < merit_at_min) then
    merit_at_min = merit
    var_at_min = var_value
  endif
  write (line, '(i5, es14.4, es10.2)') i, merit
  call out_io (s_blank$, r_name, line)

  ! Should make 1d

  if (merit <= s%global%lmdif_negligible_merit) then
    call out_io (s_blank$, r_name, 'Merit value is negligible! (Smaller than global%lmdif_negligible_merit)', &
                                   'Stopping now.')
    abort = .true.
    exit cycle_loop
  endif

  ! look for keyboard input to end optimization

  abort = tao_user_is_terminating_optimization()
  if (abort) then
    call out_io (s_blank$, r_name, line)
    exit cycle_loop
  endif

  if (at_end) exit

enddo cycle_loop

! cleanup

if (.not. abort .and. i < s%global%n_opti_cycles) then
  call out_io (s_blank$, r_name, 'Optimizer at minimum. Stopping now.')
endif

if (merit > merit_at_min) then
  call out_io (s_blank$, r_name, 'Setting to minimum.')
  call tao_set_opt_vars (var_at_min, .true.)
  merit = tao_merit()
  write (line, '(i5, es14.4, es10.2)') i+1, merit
  call out_io (s_blank$, r_name, line)
endif

if (s%opti_write_var_file) call tao_var_write (s%global%var_out_file)
deallocate (var_at_min)

end subroutine


