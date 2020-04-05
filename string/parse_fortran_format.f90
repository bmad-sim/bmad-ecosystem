!+
! Subroutine parse_fortran_format (format_str, multiplyer, power, code, width, digits)
!
! Routine to parse a Fortran edit descriptor.
! This routine assumes that format_str will be a edit descriptor for a single entity like '3f10.6'.
! format_str must only have a single edit descriptor. That is, format_str must not have any commas.
!
! Also see: real_num_fortran_format
!
! Input:
!   format_str  -- character(*): Format string.
!
! Output:
!   multiplyer  -- integer: Multiplyer prefix before the code string. EG '5I2' -> multiplyer = 5, 
!                   'es12.6' -> multiplyer = 1 (default).
!   power       -- integer: Power scale factor. EG: '3pf12.2' -> power = 3. Set to 0 if no power factor present.
!   code        -- character(*): Format code. EG 't37' -> code = 'T'. Will be upper case. Set to '' if there is a decode error.
!                   Note: If a power factor is present, code will *not* contain "p". EG: '3pf12.5' -> code = 'f'.
!   width       -- integer: Field width. EG 'A23' -> width = 23. Set to -1 if not present (EG: '34X')
!   digits      -- integer: number of digits after the decimal place for real numbers. EG '4f10.3' -> digits = 3.
!                   Set to -1 if not present.
!-

subroutine parse_fortran_format (format_str, multiplyer, power, code, width, digits)

use sim_utils_interface, dummy => parse_fortran_format

implicit none

integer multiplyer, width, digits
integer i1, i0, nn, power, n_parens
character(*) format_str, code
logical found_parens

! Init

multiplyer = 1
code = ''
width = -1
digits = -1
power = 0
n_parens = 0
nn = len(format_str)
if (nn == 0) return

! ignore leading spaces

do i1 = 1, nn
  if (format_str(i1:i1) /= ' ') exit
enddo

call check_for_open_parens (i1)
if (i1 > nn) return

! Look for multiplyer

if (index('1234567890', format_str(i1:i1)) /= 0) then
  i0 = i1
  do i1 = i0, nn
    if (index('1234567890', format_str(i1:i1)) == 0) exit
  enddo
  read (format_str(i0:i1-1), *) multiplyer
endif

! Look for power factor

call check_for_open_parens (i1, found_parens)
if (i1 > nn .or. (found_parens .and. i1 >= nn)) return

if (found_parens .and. index('1234567890', format_str(i1:i1)) /= 0) then ! If something like '3(4pf12.4)' then 
  i0 = i1
  do i1 = i0, nn
    if (index('1234567890', format_str(i1:i1)) == 0) exit
  enddo
  read (format_str(i0:i1-1), *) power
  if (i1 > nn-1) return
  if (format_str(i1:i1) /= 'p' .and. format_str(i1:i1) /= 'P') return
  i1 = i1 + 1

elseif (format_str(i1:i1) == 'p' .or. format_str(i1:i1) == 'P') then  ! If something like '4pf12.4'
  power = multiplyer
  multiplyer = 1
  i1 = i1 + 1
endif

! Look for code

if (i1 > nn .or. (found_parens .and. i1 >= nn)) return

i0 = i1
do i1 = i0, nn
  if (.not. is_alphabetic(format_str(i1:i1))) exit
enddo

if (i1 == i0) return
code = upcase(format_str(i0:i1-1))
if (no_more(i1)) return

! Look for width

i0 = i1
do i1 = i0, nn
  if (index('1234567890', format_str(i1:i1)) == 0) exit
enddo

if (i1 == i0 .or. .not. is_integer(format_str(i0:i1-1))) then
  code = ''   ! Decode error
  return
endif

read (format_str(i0:i1-1), *) width
if (no_more(i1)) return

! Look for digits

if (format_str(i1:i1) /= '.') then
  code = ''   ! Decode error
  return
endif

i0 = i1 + 1
do i1 = i0, nn
  if (index('1234567890', format_str(i1:i1)) == 0) exit
enddo

if (i1 == i0 .or. .not. is_integer(format_str(i0:i1-1))) then
  code = ''   ! Decode error
  return
endif

read (format_str(i0:i1-1), *) digits
if (no_more(i1)) return

!----------------------------------------
contains

function no_more(i1) result (nothing_else)

logical nothing_else
integer i1

!

do
  if (i1 == nn + 1) then
    nothing_else = .true.
    if (n_parens /= 0) code = ''
    return
  endif

  if (format_str(i1:i1) == ')') then
    n_parens = n_parens - 1
    i1 = i1 + 1
  else
    exit
  endif
enddo

nothing_else = (format_str(i1:i1) == '')
if (nothing_else .and. format_str(i1+1:) /= '') code = ''   ! Decode error

end function no_more

!----------------------------------------
! contains

subroutine check_for_open_parens(i1, found_parens)

integer i1
logical, optional :: found_parens

! 

if (present(found_parens)) found_parens = .false.
if (i1 > nn) return
if (format_str(i1:i1) /= '(') return

if (present(found_parens)) found_parens = .true.
n_parens = n_parens + 1
i1 = i1 + 1

end subroutine check_for_open_parens

end subroutine
