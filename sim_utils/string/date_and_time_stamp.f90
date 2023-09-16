!+
! Subroutine date_and_time_stamp (string, numeric_month, include_zone)
!
! Subroutine to return the current date and time in a character string.
!
! Input:
!   numeric_month -- logical, optional: If present and True, the month string
!                       will be a two digit number.
!   include_zone  -- logical, optional: If present and True, the time zone suffix will be added.
! Output:
!   string -- character(20) or character(26) with zone: Returns with current date and time
!             For example:
!                 call date_and_time_stamp(string)                 ! string = "1997-JUL-02 10:53:59"
!                 call date_and_time_stamp(string, .true., .true.) ! string = "1997-07-02 10:53:59 -0500"
!-

subroutine date_and_time_stamp (string, numeric_month, include_zone)

use sim_utils, except => date_and_time_stamp

implicit none

character(*) string
character(8) date
character(10) time
character(5) zone
integer values(8)
logical, optional :: numeric_month, include_zone

character(3) :: month(12) = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' ]

!

call date_and_time (date, time, zone, values)

if (logic_option(.false., numeric_month)) then
  string = date(1:4) // '-' // date(5:6) // '-' // date(7:8) // &
    ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6)
else
  string = date(1:4) // '-' // month(values(2)) // '-' // date(7:8) // &
    ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6)
endif

if (logic_option(.false., include_zone)) then
  string = trim(string) // ' ' // zone
endif

end subroutine
