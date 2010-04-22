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

use tao_mod
use tao_var_mod
use tao_lm_optimizer_mod
use tao_svd_optimizer_mod

implicit none

real(rp), allocatable, save :: var_vec(:)
real(rp) merit
integer n_data, i

character(*)  which
character(40) :: r_name = 'tao_run_cmd', my_opti

logical abort

!

call tao_set_var_useit_opt()
call tao_set_data_useit_opt()

if (all (which /= ['      ', 'de    ', 'lm    ', 'lmdif ', 'custom', 'svd   '])) then
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

! See if there are any constraints

n_data = 0
do i = lbound(s%u, 1), ubound(s%u, 1)
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
  s%u(:)%mat6_recalc_on = .false.
  s%u(ix_common_uni$)%mat6_recalc_on = .true.
endif

! Optimize...

tao_com%optimizer_running = .true.

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

tao_com%optimizer_running = .false.
if (s%global%orm_analysis) s%u(:)%mat6_recalc_on = .true.
s%u(:)%lattice_recalc = .true.


end subroutine

