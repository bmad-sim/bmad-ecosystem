!+
! Subroutine tao_lmdif_optimizer ()
!
! Subrutine to minimize the merit function by varying variables until
! the "data" as calculated from the model matches the measured data.
! 
! This subroutine is a wrapper for the mrqmin routine of Numerical Recipes.
! See the Numerical Recipes writeup for more details.
! 'lm' stands for Levenburg - Marquardt. Otherwise known as LMDIF. 
!
! Note: LM assumes 
!
! Input:
!
! Output:
!-

subroutine tao_lmdif_optimizer ()

use tao_mod
use tao_dmerit_mod
use tao_top10_mod
use tao_var_mod
use single_char_input_mod
use lmdif_mod

implicit none

type (tao_universe_struct), pointer :: u

real(rp), allocatable, save :: merit_vec(:), weight(:)
real(rp), allocatable, save :: var_delta(:), var_value(:)
real(rp) merit

integer i, j, k, n
integer n_data, n_var

logical :: finished, init_needed = .true.
logical at_end

character(20) :: r_name = 'tao_lmdif_optimizer'
character(80) :: line
character(1) char

! setup

call tao_get_vars (var_value, var_delta = var_delta, var_weight = weight)
n_var = size(var_delta)

n_data = n_var
do i = 1, size(s%u)
  n_data = n_data + count(s%u(i)%data(:)%useit_opt .and. s%u(i)%data(:)%weight /= 0)
enddo

if (allocated(merit_vec)) deallocate(merit_vec)
allocate (merit_vec(n_data))

! run optimizer mrqmin from Numerical Recipes.

finished = .false.
call initial_lmdif

merit = tao_merit()

call out_io (s_blank$, r_name, '   Loop      Merit')

do i = 1, s%global%n_opti_cycles

  merit_vec(1:n_var) = sqrt(weight) * var_delta
  k = n_var
  do n = 1, size(s%u)
    u => s%u(n)
    do j = 1, size(u%data)
      if (.not. u%data(j)%useit_opt) cycle
      if (u%data(j)%weight == 0) cycle
      k = k + 1
      merit_vec(k) = sqrt(u%data(j)%weight) * u%data(j)%delta_merit
    enddo
  enddo

  call suggest_lmdif (var_value, merit_vec, s%global%lmdif_eps, s%global%n_opti_cycles, at_end)
  call tao_set_vars (var_value)
  write (line, '(i5, es14.4, es10.2)') i, tao_merit()
  call out_io (s_blank$, r_name, line)

  if (at_end) then
    s%global%optimizer_running = .false.
    call out_io (s_blank$, r_name, 'Optimizer at minimum. Stopping now.')
    exit
  endif

#ifndef CESR_WINCVF
  ! look for keyboard input to end optimization

  do
    call get_tty_char (char, .false., .false.) 
    if (char == '.') then
      s%global%optimizer_running = .false.
      finished = .true.
      call suggest_lmdif (var_value, merit_vec, s%global%lmdif_eps, i, at_end)
      call tao_set_vars (var_value)
      write (line, '(i5, es14.4, es10.2)') i, tao_merit()
      call out_io (s_blank$, r_name, line)
      call out_io (s_blank$, r_name, 'Optimizer stop signal detected.', &
                                                             'Stopping now.')
      call tao_var_write (s%global%var_out_file)
      return
    endif
    if (char == achar(0)) exit   ! only exit if there is no more input
  enddo
#endif

enddo

call tao_var_write (s%global%var_out_file)

end subroutine


