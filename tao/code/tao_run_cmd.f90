!+
! Subroutine tao_run_cmd (which, abort)
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches a constraint or measured data.
!
! Input:
!   which -- Character(*): which optimizer to use. 
!             ' '        -- Same as last time
!             'de'       -- Differential Evolution.
!             'lm'       -- Levenberg - Marquardt (aka lmdif).
!             'custom'   -- Custom routine.
!
! Output:
!  abort -- Logical: Set True if the run was aborted by the user, an at minimum 
!             condition, a singular matrix condition, etc.. False otherwise.
!-

subroutine tao_run_cmd (which, abort)

use tao_interface, dummy => tao_run_cmd
use tao_lm_optimizer_mod, only: tao_lm_optimizer
use tao_svd_optimizer_mod, only: tao_svd_optimizer
use tao_geodesic_lm_optimizer_mod, only: tao_geodesic_lm_optimizer
!MPI use tao_mpi_mod

implicit none

type (tao_universe_struct), pointer :: u
type (tao_universe_calc_struct) u_calc(lbound(s%u,1):ubound(s%u,1))
real(rp), allocatable :: var_vec(:)
real(rp) merit0, merit
integer n_data, i, j, iu0, iu1

character(*)  which
character(40) :: r_name = 'tao_run_cmd', my_opti

logical abort

!

call tao_set_var_useit_opt()
call tao_set_data_useit_opt()

if (which /= '') then
  if (all (which /= tao_optimizer_name)) then
    call out_io (s_error$, r_name, 'OPTIMIZER NOT RECOGNIZED: ' // which)
    return
  endif
  s%global%optimizer = which
endif

call out_io (s_blank$, r_name, 'Optimizing with: ' // s%global%optimizer)
call out_io (s_blank$, r_name, "Type ``.'' to stop the optimizer before it's finished.")

call tao_get_opt_vars (var_vec)
if (size(var_vec) == 0) then
  call out_io (s_fatal$, r_name, 'No variables to vary!')
  abort = .true.
  return
endif

! Do not do radiation_integrals calc if not needed

s%com%have_datums_using_expressions = .false.
s%com%all_merit_weights_positive = (all(s%var%weight >= 0))

iu0 = lbound(s%u, 1); iu1 = ubound(s%u, 1)
do i = iu0, iu1
  u => s%u(i)
  u_calc(i) = u%calc

  u%calc%rad_int_for_data        = .false.
  u%calc%rad_int_for_plotting    = .false.
  u%calc%chrom_for_data          = .false.
  u%calc%chrom_for_plotting      = .false.
  u%calc%lat_sigma_for_data      = .false.
  u%calc%lat_sigma_for_plotting  = .false.
  u%calc%spin_matrices           = .false.
  u%calc%srdt_for_data           = 0

  do j = 1, size(u%data)
    if (.not. u%data(j)%useit_opt) cycle
    if (tao_rad_int_calc_needed(u%data(j)%data_type, u%data(j)%data_source)) u%calc%rad_int_for_data = .true.
    if (tao_chrom_calc_needed(u%data(j)%data_type, u%data(j)%data_source)) u%calc%chrom_for_data = .true.
    if (tao_lat_sigma_calc_needed(u%data(j)%data_type, u%data(j)%data_source)) u%calc%lat_sigma_for_data = .true.
    if (tao_spin_matrices_calc_needed(u%data(j)%data_type, u%data(j)%data_source)) u%calc%spin_matrices = .true.
    if (substr(u%data(j)%data_type,1,11) == 'expression:') s%com%have_datums_using_expressions = .true.
    u%calc%srdt_for_data = max(tao_srdt_calc_needed(u%data(j)%data_type, u%data(j)%data_source), u%calc%srdt_for_data)
  enddo

  s%com%all_merit_weights_positive = (s%com%all_merit_weights_positive .and. all(u%data%weight >= 0))
enddo

! See if there are any constraints

n_data = 0
do i = iu0, iu1
  n_data = n_data + count(s%u(i)%data(:)%useit_opt)
enddo

if (n_data == 0) then
  call out_io (s_error$, r_name, 'No data constraints defined for the merit function!')
  abort = .true.
  return
endif

! Optimize...

s%com%optimizer_running = .true.
merit0 = tao_merit()

do i = 1, s%global%n_opti_loops
  select case (s%global%optimizer)
  case ('de') 
    call tao_de_optimizer (abort)

  case ('lm') 
    call tao_lm_optimizer (abort)

  case ('lmdif')
    call tao_lmdif_optimizer (abort)

  case ('svd')
    call tao_svd_optimizer (abort)

  case ('geodesic_lm')
    call tao_geodesic_lm_optimizer (abort)

  case ('custom')
    if (.not. associated(tao_hook_optimizer_ptr)) then
      call out_io(s_error$, r_name, 'TAO_HOOK_OPTIMIZER_PTR NOT ASSOCIATED!')
      return
    endif
    call tao_hook_optimizer_ptr (abort)
  end select

  if (abort) exit

  merit = tao_merit()
  if (s%global%merit_stop_value /= 0 .and. s%global%merit_stop_value > merit) then
    call out_io (s_info$, r_name, 'Merit value below global%merit_stop_value. Stopping Optimization.')
    exit
  endif

  if (merit0 - merit <= abs(merit) * s%global%dmerit_stop_value) then
    call out_io (s_info$, r_name, 'Fractional change in Merit value in one loop below global%dmerit_stop_value. Stopping Optimization.')
    exit
  endif

  merit0 = merit
enddo

! We need a lattice recalc one last time since data not used in the 
! optimization has not been updated.

s%com%optimizer_running = .false.
s%u(:)%calc = u_calc

call tao_turn_on_special_calcs_if_needed_for_plotting ()

s%u(:)%calc%lattice = .true.

end subroutine

