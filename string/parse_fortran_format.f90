!+
! Subroutine parse_fortran_format (format_str, multiplyer, code, width, digits)
!
! Routine to parse a Fortran edit descriptor.
! This routine assumes that format_str will be a edit descriptor for a single entity like '3f10.6'.
! Strings multiple entities (any string with a comma) are not allowed.
!
! Also see: real_num_fortran_format
!
! Input:
!   format_str  -- character(*): Format string.
!
! Output:
!   multiplyer  -- integer: Multiplyer prefix before the code string. EG '5I2' -> multiplyer = 5, 
!                   'es12.6' -> multiplyer = 1 (default).
!   code        -- character(*): Format code. EG 'T37' -> code = 'T'. Set to '' if there is a decode error.
!   width       -- integer: Field width. EG 'A23' -> width = 23. Set to -1 if not present (EG: '34X')
!   digits      -- integer: number of digits after the decimal place for real numbers. EG '4f10.3' -> digits = 3.
!                   Set to -1 if not present.
!-

subroutine parse_fortran_format (format_str, multiplyer, code, width, digits)

use sim_utils_interface, dummy => parse_fortran_format

implicit none

integer multiplyer, width, digits
integer i, i0, nn
character(*) format_str, code

! Init

multiplyer = 1
code = ''
width = -1
digits = -1
nn = len(format_str)

! ignore leading spaces

do i = 1, nn
  if (format_str(i:i) /= ' ') exit
enddo

if (i == nn + 1) return

! Look for multiplyer

i0 = i
do i = i0, nn
  if (index('1234567890', format_str(i:i)) == 0) exit
enddo

if (i /= i0) then
  if (.not. is_integer(format_str(i0:i-1))) return
  read (format_str(i0:i-1), '(i)') multiplyer
endif

! Look for code

i0 = i
do i = i0, nn
  if (index('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ', format_str(i:i)) == 0) exit
enddo

if (i == i0) return
code = format_str(i0:i-1)
if (no_more(i)) return

! Look for width

i0 = i
do i = i0, nn
  if (index('1234567890', format_str(i:i)) == 0) exit
enddo

if (i == i0 .or. .not. is_integer(format_str(i0:i-1))) then
  code = ''   ! Decode error
  return
endif

read (format_str(i0:i-1), '(i)') width
if (no_more(i)) return

! Look for digits

if (format_str(i:i) /= '.') then
  code = ''   ! Decode error
  return
endif

i0 = i + 1
do i = i0, nn
  if (index('1234567890', format_str(i:i)) == 0) exit
enddo

if (i == i0 .or. .not. is_integer(format_str(i0:i-1))) then
  code = ''   ! Decode error
  return
endif

read (format_str(i0:i-1), '(i)') digits
if (no_more(i)) return


!----------------------------------------
contains

function no_more(i1) result (nothing_else)

logical nothing_else
integer i1

!

if (i1 == nn + 1) then
  nothing_else = .true.
  return
endif

nothing_else = (format_str(i1:i1) == '')
if (nothing_else .and. format_str(i1+1:) /= '') code = ''   ! Decode error

end function

end subroutine
