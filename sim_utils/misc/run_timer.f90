!+
! Subroutine run_timer (command, time)
!
! Subroutine to time a program or return the absolute time in seconds
! from some reference date (default is January 1, 1970).
! Note: Deprecated command = 'STOP' is the same as command = 'READ'.
!
! Input:
!   command -- Character(*): 'START', 'READ', 'ABS', 'REF'
!   time    -- Real(rp), optional: Set reference year. Used with command = 'REF'.
!                 Month and day are always January 1. Default is 1970.0_rp. 
!
! Output:
!   time    -- Real(rp), optional: 
!                 If command = 'READ': Elapsed time from last 'START' in seconds.
!                 If command = 'ABS':  Time in seconds from reference year.
!                 If not set using 'REF', Reference is January 1, 1970.
!
! Example:
!     call run_timer ('START')          ! Start (reset) the timer.
!     call run_timer ('READ', time)     ! time = Elapsed time from last 'START'
!     call run_timer ('REF', 2000.0_rp) ! Reference date is January 1, 2000.
!     call run_timer ('ABS', time)      ! Time from reference date
!-

subroutine run_timer(command, time)

use precision_def

implicit none

integer, save :: t0(8), t1(8), dt(8), ref_year = 1970
logical, save :: timer_running = .false.
real(rp) time
character(*) command

!

select case (command)
case ('START') 
  call date_and_time (values = t0)
  t0(3) = julian_day(t0(1), t0(2), t0(3))
  timer_running = .true.

case ('READ', 'STOP')
  if (.not. timer_running) then
    print *, 'ERROR IN RUN_TIMER: "START" NEVER DONE.'
    time = 1.0d10   ! some non-zero time to not bomb program.
  endif

  call date_and_time (values = t1)
  t1(3) = julian_day(t1(1), t1(2), t1(3))
  dt = t1 - t0
  time = 60.0_rp*(60.0_rp*(24.0_rp*dt(3) + dt(5)) + dt(6)) + dt(7) + dt(8) / 1000.0_rp

case ('ABS')
  call date_and_time (values = t1)
  t1(3) = julian_day(t1(1), t1(2), t1(3)) - julian_day(ref_year, 0, 0)
  time = 60.0_rp*(60.0_rp*(24.0_rp*t1(3) + t1(5)) + t1(6)) + t1(7) + t1(8) / 1000.0_rp

case ('REF')
  ref_year = nint(time)

case default
  print *, 'ERROR IN RUN_TIMER SUBROUTINE: INVALID COMMAND: ', command
  if (global_com%exit_on_error) call err_exit

end select

!------------------------------------------------------------------------
contains

! Computes the Julian date given a Gregorian calendar date (year, month, day).
! Algorithm by Fliegel and van Flandern (1968)

function julian_day (year, month, day) result (jday)

implicit none

integer year, month, day, jday

! 

jday = Day - 32075 + 1461*(YEAR+4800+(MONTH-14)/12)/4 + &
             367*(MONTH-2-(MONTH-14)/12*12)/12 - 3*((YEAR+4900+(MONTH-14)/12)/100)/4

end function

end subroutine

