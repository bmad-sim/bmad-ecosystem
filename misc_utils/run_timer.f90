!+
! Subroutine run_timer (command, time)
!
! Subroutine to time a program.
!
! Modules needed:
!   use dcslib_interface
!
! Input:
!   command -- Character(*): 'START', 'READ', or 'STOP'
!   time    -- Real(rp), optional: Elapsed time in seconds.
!
! Example:
!     call run_timer ('START')       ! Start the timer
!     call run_timer ('READ', time)  ! Read the timer, Keep running.
!     call run_timer ('STOP', time)  ! Read the timer and do a reset.
!-

#include "CESR_platform.inc"

subroutine run_timer(command, time)

  use precision_def

  implicit none

  integer t0(8), t1(8), dt(8), days_in_year0, days_in_year1
  real(rp) time
  character(*) command

!

  select case (command)
  case ('START') 
    call date_and_time (values = t0)
    t0(3) = julian_day(t0(1), t0(2), t0(3), days_in_year0)

  case ('READ', 'STOP')
    call date_and_time (values = t1)
    t1(3) = julian_day(t1(1), t1(2), t1(3), days_in_year1)
    dt = t1 - t0
    dt(3) = dt(3) + dt(1) * days_in_year0
    time = 60*(60*(24*dt(3) + dt(5)) + dt(6)) + dt(7) + dt(8) / 1000.0

  case default
    print *, 'ERROR IN RUN_TIMER SUBROUTINE: INVALID COMMAND: ', command
    time = 1.0e10   ! some non-zero time to not bomb program.

  end select

!------------------------------------------------------------------------
contains

function julian_day (year, month, day, days_in_year) result (jday)

integer year, month, day, jday, days_in_year
integer :: month_len(12) = &
               (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

!

month_len(2) = 28
if (mod(year,   4) == 0) month_len(2) = 29 ! leap year
if (mod(year, 100) == 0) month_len(2) = 28 ! but not if century
if (mod(year, 400) == 0) month_len(2) = 29 ! but but yes it is a leap year

jday = day + sum(month_len(1:month-1))
days_in_year = 365 + month_len(2) - 28

end function

end subroutine

