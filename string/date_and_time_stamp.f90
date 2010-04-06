!+
! Subroutine date_and_time_stamp (string, numeric_month)
!
! Subroutine to return the current date and time in a character string.
!
! Input:
!   numeric_month -- Logical, optional: If present and True, the month string
!                       will be a two digit number.
! Output:
!   string -- Character(20): Returns with current date and time of the form:
!                 "yyyy-mmm-dd hh:mm:ss"
!             For example: 
!                 "1997-JUL-02 10:53:59"
!             If the numeric_month argument is present then string will be:
!                 "1997-07-02 10:53:59 "
!-

subroutine date_and_time_stamp (string, numeric_month)

  use sim_utils

  implicit none

  character(*) string
  character date*8, time*10, zone*5
  integer values(8)
  logical, optional :: numeric_month

  character(3) :: month(12) = (/ 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', &
                                 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC' /)

!

  call date_and_time (date, time, zone, values)

  if (logic_option(.false., numeric_month)) then
    string = date(1:4) // '-' // date(5:6) // '-' // date(7:8) // &
      ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6)
  else
    string = date(1:4) // '-' // month(values(2)) // '-' // date(7:8) // &
      ' ' // time(1:2) // ':' // time(3:4) // ':' // time(5:6)
  endif

end subroutine
