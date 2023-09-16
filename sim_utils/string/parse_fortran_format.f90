!+
! Subroutine parse_fortran_format (format_str, n_repeat, power, descrip, width, digits)
!
! Routine to parse a Fortran edit descriptor.
!
! This routine assumes that format_str will be a edit descriptor for a single entity like '3f10.6'.
! format_str must only have a single edit descriptor.
!
! Also see: real_num_fortran_format
!
! Input:
!   format_str  -- character(*): Format string.
!
! Output:
!   n_repeat    -- integer: Repeat Count prefix before the descrip string. EG '5I2' -> n_repeat = 5, 
!                   'es12.6' -> n_repeat = 1 (default). '3x' -> n_repeat = 3, width = 1.
!   power       -- integer: Power scale factor. EG: '3pf12.2' -> power = 3. Set to 0 if no power factor present.
!   descrip     -- character(*): Format descrip. EG 't37' -> descrip = 'T'. Will be upper case. Set to '' if there is a decode error.
!                   Note: If a power factor is present, descrip will *not* contain "p". EG: '3pf12.5' -> descrip = 'F'.
!                   descrip will be blank if there is an error.
!   width       -- integer: Field width. EG 'A23' -> width = 23. Set to -1 if not present (EG: '34X')
!   digits      -- integer: number of digits after the decimal place for real numbers. EG '4f10.3' -> digits = 3.
!                   Set to -1 if not present.
!-

subroutine parse_fortran_format (format_str, n_repeat, power, descrip, width, digits)

use sim_utils_interface, dummy => parse_fortran_format

implicit none

integer n_repeat, width, digits
integer n, i1, nn, power, n_parens, ios
character(*) format_str, descrip
character(40) fmt
logical found_parens

! Init

n_repeat = 1
descrip = ''
width = -1
digits = -1
power = 0
n_parens = 0

call string_trim(format_str, fmt, nn)
if (nn == 0) return

! Example formats:
!   "7(-3p, f12.3)"  "7pes12.1" "2pes12.1" "2(-3pf12.1)"

! Look for "p"

i1 = max(index(fmt, 'p'), index(fmt, 'P'))
if (i1 /= 0) then
  do n = i1-1, 1, -1
    if (index('-1234567890', fmt(n:n)) == 0) exit
  enddo
  if (n == i1-1) return
  read(fmt(n+1:i1-1), *, iostat = ios) power
  if (ios /= 0) return
  fmt(n+1:) = fmt(i1+1:)
endif

! Remove commas and parens.

n = 1
nn = len_trim(fmt)

do i1 = 1, nn
  select case (fmt(n:n))
  case (' ', ',', '(', ')')
    fmt(n:) = fmt(n+1:)
  case default
    n = n + 1
    if (n > nn) exit
  end select
enddo

! Look for n_repeat

do n = 1, nn
  if (index('1234567890', fmt(n:n)) == 0) exit
enddo

if (n > 1) then
  read (fmt(:n-1), *, iostat = ios) n_repeat
  if (ios /= 0) return
  fmt = fmt(n:)
endif

! Look for descrip

do n = 1, nn
  if (.not. is_alphabetic(fmt(n:n))) exit
enddo

if (n == 1) return
descrip = upcase(fmt(:n-1))
fmt = fmt(n:)
if (fmt == '') return

! Look for width

do n = 1, nn
  if (index('1234567890', fmt(n:n)) == 0) exit
enddo

if (n > 1) then
  read (fmt(1:n-1), *, iostat = ios) width
  if (ios /= 0) then
    descrip = ''
    return
  endif
  fmt = fmt(n:)
endif

! Look for digits

if (fmt(1:1) /= '.') return
read (fmt(2:), *, iostat = ios) digits
if (ios /= 0) descrip = ''


end subroutine
