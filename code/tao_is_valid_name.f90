!+
! Function tao_is_valid_name (name, why_invalid) result (is_valid)
!
! Routine to check if a variable or data name is valid.
!
! Input:
!   name          -- character(*): Name to be checked.
!
! Output:
!   why_invalid   -- character(*): Why invalid description.
!   is_valid      -- logical: True if valid. False otherwise.
!-

function tao_is_valid_name (name, why_invalid) result (is_valid)

use sim_utils

implicit none

character(*) name, why_invalid
logical is_valid

!

is_valid = .false.
why_invalid = ''

if (index(name, '.') /= 0) then
  why_invalid = 'NAME CANNOT CONTAIN A PERIOD CHARACTER: ' // quote(name)
  return
endif

if (name == 'ALL' .or. name == 'all') then
  why_invalid = 'NAME CANNOT BE "ALL" OR "all"'
  return
endif

if (is_integer(name)) then
  why_invalid = 'NAME CANNOT BE AN INTEGER: ' // quote(name)
  return
endif

if (index(trim(name), ' ') /= 0 .or. index(trim(name), char(9)) /= 0) then
  why_invalid = 'NAME CANNOT CONTAIN A BLANK OR TAB CHARACTER: ' // quote(name)
  return
endif

is_valid = .true.

end function tao_is_valid_name
