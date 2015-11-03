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

use tao_mod, dummy => tao_run_cmd
use tao_lm_optimizer_mod, only: tao_lm_optimizer
use tao_svd_optimizer_mod, only: tao_svd_optimizer
use tao_var_mod, only: tao_get_opt_vars
use tao_geodesic_lm_optimizer_mod, only: tao_geodesic_lm_optimizer
!MPI use tao_mpi_mod

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable, save :: var_vec(:)
real(rp) merit
integer n_data, i, j, iu0, iu1

character(*)  which
character(40) :: r_name = 'tao_run_cmd', my_opti

logical abort
logical, allocatable :: do_rad_int_data(:), do_chrom_data(:)

!

call tao_set_var_useit_opt()
call tao_set_data_useit_opt()

if (all (which /= ['           ', 'de         ', 'lm         ', 'lmdif      ', &
                   'custom     ', 'svd        ', 'geodesic_lm'])) then
  call out_io (s_error$, r_name, 'OPTIMIZER NOT RECOGNIZED: ' // which)
  return
endif

if (which /= ' ') s%global%optimizer = which
call out_io (s_blank$, r_name, 'Optimizing with: ' // s%global%optimizer)
call out_io (s_blank$, r_name, &
              "Type ``.'' to stop the optimizer before it's finished.")

call tao_get_opt_vars (var_vec)
if (size(var_vec) == 0) then
  call out_io (s_fatal$, r_name, 'No variables to vary!')
  abort = .true.
  return
endif

! Do not do radiation_integrals calc if not needed

s%com%have_datums_using_expressions = .false.

iu0 = lbound(s%u, 1); iu1 = ubound(s%u, 1)
allocate (do_rad_int_data(iu0:iu1), do_chrom_data(iu0:iu1))

do i = iu0, iu1
  u => s%u(i)

  do_rad_int_data(i) = u%calc%rad_int_for_data
  do_chrom_data(i)   = u%calc%chrom_for_data

  u%calc%rad_int_for_data     = .false.
  u%calc%rad_int_for_plotting = .false.
  u%calc%chrom_for_data     = .false.
  u%calc%chrom_for_plotting = .false.

  do j = 1, size(u%data)
    if (.not. u%data(j)%useit_opt) cycle
    if (tao_rad_int_calc_needed(u%data(j)%data_type, u%data(j)%data_source)) &
                                                       u%calc%rad_int_for_data = .true.
    if (tao_chrom_calc_needed(u%data(j)%data_type, u%data(j)%data_source)) &
                                                       u%calc%chrom_for_data = .true.
    if (u%data(j)%data_type(1:11) == 'expression:') s%com%have_datums_using_expressions = .true.
  enddo
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

! Save time with orm analysis by turning off transfer matrix calc in
! everything but the common universe.

if (s%global%orm_analysis) then
  s%u(:)%calc%mat6 = .false.
  s%u(ix_common_uni$)%calc%mat6 = .true.
  s%u(:)%calc%track = .false.
  s%u(ix_common_uni$)%calc%track = .true.
endif

! Optimize...

s%com%optimizer_running = .true.

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
    call tao_hook_optimizer (abort)

  end select

  if (abort) exit
  if (s%global%merit_stop_value > 0 .and. s%global%merit_stop_value > tao_merit()) then
    call out_io (s_info$, r_name, 'Merit value below global%merit_stop_value. Stopping optimization.')
    exit
  endif

enddo

! We need a lattice recalc one last time since data not used in the 
! optimization has not been updated.

s%com%optimizer_running = .false.

s%u(:)%calc%rad_int_for_data = do_rad_int_data
s%u(:)%calc%chrom_for_data = do_chrom_data
call tao_turn_on_chrom_or_rad_int_calcs_if_needed_for_plotting ()
deallocate (do_rad_int_data, do_chrom_data)

if (s%global%orm_analysis) s%u(:)%calc%mat6 = .true.
s%u(:)%calc%lattice = .true.


end subroutine

