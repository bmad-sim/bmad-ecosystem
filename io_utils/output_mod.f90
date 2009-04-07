#include "CESR_platform.inc"

module output_mod

use cesr_interface
use precision_def
use parallel_mod
use utilities_mod

! Message levels: Status level flags for messages.

integer, parameter :: s_blank$   = -1 ! Information message. The routine name is not printed.
integer, parameter :: s_info$    = 0  ! Informational message.
integer, parameter :: s_dinfo$   = 1  ! Info message (w/timestamp).
integer, parameter :: s_success$ = 2  ! Successful completion.
integer, parameter :: s_warn$    = 3  ! Warning of a possible problem.
integer, parameter :: s_dwarn$   = 5  ! Warning of a possible problem (w/timestamp).
integer, parameter :: s_error$   = 7  ! An error as occurred [EG: bad user input] (w/ timestamp).
integer, parameter :: s_fatal$   = 8  ! A fatal error has occurred so that computations
                                      ! cannot be continued but the program will try to
                                      ! reset itself and keep running (w/timestamp).
integer, parameter :: s_abort$   = 9  ! A severe error has occurred and
                                      ! the program is being aborted (w/timestamp).


type output_mod_com_struct
  logical :: do_print(-1:9) = .true.
  integer :: file_unit(-1:9) = 0
#if defined (CESR_VMS)
  integer :: indent_num(-1:9) = (/ 1, 1, 5, 5, 5, 5, 5, 5, 5, 5, 5 /)
#else
  integer :: indent_num(-1:9) = (/ 0, 0, 4, 4, 4, 4, 4, 4, 4, 4, 4 /)
#endif
end type

type (output_mod_com_struct), save, private :: output_com

private output_line4, output_int, output_real, output_logical
private header_io, find_format, output_lines, insert_numbers

!+
! Subroutine out_io
!
! Subroutine to print to the terminal for command line type programs.
! The idea is that for programs with a gui this routine can be easily
! replaced with another routine.
!
! This routine is an overloaded name for:
!   output_real (level, routine_name, line, r_num)
!   output_int (level, routine_name, line, i_num)
!   output_logical (level, routine_name, line, l_num)
!   output_lines (level, routine_name, lines, r_array, i_array, l_array)
!   output_line4 (level, routine_name, &
!                     line1, line2, line3, line4, r_array, i_array, l_array)
!
! Numbers are encoded in lines using the syntax "\<fmt>\" 
! where <fmt> is the desired format. For example:
!   r_array(1:4) = (/ 1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp /)
!   call out_io (s_info$, routine_name, "4 Numbers: \4f8.3\", r_array = r_array)
! A valid format
!   a) Must not contain any spaces.
!   b) Must start with an optional number followed by a letter:
!        'f', 'e', 'i', 'l', 'g', 'o', or. 'z'
! Anything that does not look like a valid format construct is ignored.
! Thus a literal "\" is allowed if it doesn't look like a valid format.
!
! To direct output to a file use the routine:
!   output_direct (do_print, file_unit)
!
! Modules needed:
!   use output_mod
!
! Input:
!   level -- Integer: Status level flags for messages.
!       s_blank$   -   Information message. The routine name is not printed.
!       s_info$    -   Informational message.
!       s_dinfo$   -   Info message (w/timestamp).
!       s_success$ -   Successful completion.
!       s_warn$    -   Warning of a possible problem.
!       s_dwarn$   -   Warning of a possible problem (w/timestamp).
!       s_error$   -   An error as occurred [EG: bad user input] (w/ timestamp).
!       s_fatal$   -   A fatal error has occurred so that computations
!                        cannot be continued but the program will try to
!                        reset itself and keep running (w/timestamp).
!       s_abort$   -   A severe error has occurred and
!                        the program is being aborted (w/timestamp).
!   routine_name -- Character(*): Name of the calling routine.
!   line         -- Character(*), Line to print.
!   lines(:)     -- Character(*), Lines to print.
!   line1        -- Character(*): First line to print.
!   line2        -- Character(*), optional: Second line to print.
!   line3        -- Character(*), optional: Third line to print.
!   line4        -- Character(*), optional: Firth line to print.
!   r_num        -- Real(rp): Real Number to print.
!   i_num        -- Integer: Integer to print.
!   l_num        -- Logical: Logical to print.
!   r_array(:)   -- Real(rp), optional: Real numbers to print.
!   i_array(:)   -- Integer, optional: Integer numbers to print.
!   r_array(:)   -- Logical, optional: Logicals to print.
!-

interface out_io
  module procedure output_line4
  module procedure output_lines
  module procedure output_real
  module procedure output_int
  module procedure output_logical
end interface

contains

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine output_direct (file_unit, do_print, min_level, max_level)
!
! Subroutine to set where the output goes when out_io is called.
! Output may be sent to the terminal screen, written to a file, or both.
! Once set, the settings remain until the next call to output_direct.
! 
! Settings can be made on a message status level by level basis.
! "getf out_io" will show a list of the message status levels.
!
! Modules needed:
!   use output_mod
!
! Input:
!   file_unit -- Integer, optional: Unit number for writing to a file. 
!                  0 => No writing (initial setting).
!   do_print  -- Logical, optional: If True (initial setting) then 
!                  print output at the TTY.
!   min_level -- Integer, optional: Minimum message status level to apply to. 
!                  Default is s_blank$
!   max_level -- Integer, optional: Maximum message status level to apply to. 
!                  Default is s_abort$
!-

subroutine output_direct (file_unit, do_print, min_level, max_level)

implicit none

logical, optional :: do_print
integer, optional :: file_unit, min_level, max_level
integer i

!

do i = integer_option(s_blank$, min_level), integer_option(s_abort$, max_level)
  if (present(do_print)) output_com%do_print(i) = do_print
  if (present(file_unit)) output_com%file_unit(i) = file_unit
enddo

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine out_io_line_out (line, level, indent)
! 
! Subroutine to print and/or write a line to a file.
!
! Input:
!   line       -- Character(*): Line to be outputted.
!   level      -- Integer: Message level.
!   indent     -- Logical, optional: If True then indent the line. Default is True.
!                   The number of spaces is set by the output_direct routine.
!-

subroutine out_io_line_out (line, level, indent)

implicit none

character(*) line
integer level, ix
logical, optional :: indent
character(80) :: blank = ''

!

if (logic_option(.true., indent)) then
  ix = output_com%indent_num(level)
  if (output_com%do_print(level)) &
            write (*, '(2a)') blank(1:ix), trim(line)
  if (output_com%file_unit(level) /= 0) &
            write (output_com%file_unit(level), '(2a)') blank(1:ix), trim(line)
else
#if defined (CESR_VMS)
  if (output_com%do_print(level)) write (*, '(1x, a)') trim(line)
#else
  if (output_com%do_print(level)) write (*, '(a)') trim(line)
#endif
  if (output_com%file_unit(level) /= 0) write (output_com%file_unit(level), '(a)') trim(line)
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine output_real (level, routine_name, line, r_num)
!
! Subroutine to print to the terminal for command line type programs.
! This routine is overloaded by the routine: out_io. See out_io for more details.
!-

subroutine output_real (level, routine_name, line, r_num)

implicit none

character(*) routine_name, line
character(20) fmt
character(len(line)+100) this_line

real(rp) r_num
integer level
integer ix1, ix2, n_prefix
logical found

!

if (global_rank /= 0) return
call header_io (level, routine_name)

call find_format (line, n_prefix, fmt, ix1, ix2, found)
if (found) then
  fmt = '(a, ' // trim(fmt) // ', a)'
  write (this_line, fmt) line(:ix1), r_num, trim(line(ix2:))
  call out_io_line_out (this_line, level)
else
  call out_io_line_out (line, level)
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine output_int (level, routine_name, line, i_num)
!
! Subroutine to print to the terminal for command line type programs.
! This routine is overloaded by the routine: out_io. See out_io for more details.
!-

subroutine output_int (level, routine_name, line, i_num)

implicit none

character(*) routine_name, line
character(20) fmt
character(len(line)+100) this_line

integer i_num
integer level
integer ix1, ix2, n_prefix
logical found

!

if (global_rank /= 0) return
call header_io (level, routine_name)

call find_format (line, n_prefix, fmt, ix1, ix2, found)
if (found) then
  fmt = '(a, ' // trim(fmt) // ', a)'
  write (this_line, fmt) line(:ix1), i_num, trim(line(ix2:))
  call out_io_line_out(this_line, level)
else
  call out_io_line_out(line, level)
endif


end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine output_logical (level, routine_name, line, l_num)
!
! Subroutine to print to the terminal for command line type programs.
! This routine is overloaded by the routine: out_io. See out_io for more details.
!-

subroutine output_logical (level, routine_name, line, l_num)

implicit none

character(*) routine_name, line
character(20) fmt
character(len(line)+100) this_line

logical l_num
integer level
integer ix1, ix2, n_prefix
logical found

!

if (global_rank /= 0) return
call header_io (level, routine_name)

call find_format (line, n_prefix, fmt, ix1, ix2, found)
if (found) then
  fmt = '(a, ' // trim(fmt) // ', a)'
  write (this_line, fmt) line(:ix1), l_num, trim(line(ix2:))
  call out_io_line_out(this_line, level)
else
  call out_io_line_out(line, level)
endif

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine output_line4 (level, routine_name, &
!                     line1, line2, line3, line4, r_array, i_array, l_array)
!
! Subroutine to print to the terminal for command line type programs.
! This routine is overloaded by the routine: out_io. See out_io for more details.
!-

subroutine output_line4 (level, routine_name, line1, line2, line3, line4, &
                                                     r_array, i_array, l_array)

implicit none

character(*) routine_name, line1
character(*), optional :: line2, line3, line4
character(40) fmt

real(rp), optional :: r_array(:)
integer, optional :: i_array(:)
logical, optional :: l_array(:)

integer level, nr, ni, nl

!

if (global_rank /= 0) return
call header_io (level, routine_name)

nr = 0; ni = 0; nl = 0  ! number of numbers used.

call insert_numbers (level, fmt, nr, ni, nl, line1, r_array, i_array, l_array)
if (.not. present(line2)) return
call insert_numbers (level, fmt, nr, ni, nl, line2, r_array, i_array, l_array)
if (.not. present(line3)) return
call insert_numbers (level, fmt, nr, ni, nl, line3, r_array, i_array, l_array)
if (.not. present(line4)) return
call insert_numbers (level, fmt, nr, ni, nl, line4, r_array, i_array, l_array)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine output_lines (level, routine_name, lines, r_array, i_array, l_array)
!
! Subroutine to print to the terminal for command line type programs.
! This routine is overloaded by the routine: out_io. See out_io for more details.
!-

subroutine output_lines (level, routine_name, lines, r_array, i_array, l_array)

implicit none

character(*) routine_name, lines(:)
character(40) fmt

real(rp), optional :: r_array(:)
integer, optional :: i_array(:)
logical, optional :: l_array(:)

integer level, i, nr, ni, nl

!

if (global_rank /= 0) return
call header_io (level, routine_name)

nr = 0; ni = 0; nl = 0  ! number of numbers used.

do i = 1, size(lines)
  call insert_numbers (level, fmt, nr, ni, nl, lines(i), r_array, i_array, l_array)
enddo

end subroutine


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine insert_numbers (level, fmt, nr, ni, nl, line_in, r_array, i_array, l_array)
!
! Subroutine to insert numbers into a string and print it.
! This is an internal subroutine not meant for general use.
!
! r_array, i_array, and l_array are encoded in lines using the syntax "\<fmt>\" 
! where <fmt> is the desired format. For example:
!   r_array(1:4) = (/ 1.0_rp, 2.0_rp, 3.0_rp, 4.0_rp /)
!   lines(1) = "4 Numbers: \4f8.3\"
! To use a literal "\" in line use the syntax "\\"
!
! Input:
!   level        -- Integer: Message status level.
!   fmt          -- Character(*): Format for line.
!   nr           -- Integer: Index to next available element in r_array.
!   ni           -- Integer: Index to next available element in i_array.
!   nl           -- Integer: Index to next available element in l_array.
!   line_in      -- Character(*): Input line.
!   r_array(:)   -- Real(rp), optional: Real numbers to print.
!   i_array(:)   -- Integer, optional: Integer numbers to print.
!   l_array(:)   -- Logical, optional: Logicals to print.
!
! Output:
!   nr           -- Integer: Index of last used element in r_array.
!   ni           -- Integer: Index of last used element in i_array.
!   nl           -- Integer: Index of last used element in l_array.
!-

subroutine insert_numbers (level, fmt, nr, ni, nl, line_in, r_array, i_array, l_array)

implicit none

character(*) line_in
character(*) fmt
character(40) fmt2
character(1) descrip
character(len(line_in)+100) this_line

real(rp), optional :: r_array(:)
integer, optional :: i_array(:)
logical, optional :: l_array(:)

integer level
integer nr, ni, nl
integer ix1, ix2, nn
logical found

! 

this_line = line_in

if (any ( (/ present(r_array), present(i_array), present(l_array) /) )) then

  do

    call find_format (this_line, nn, fmt2, ix1, ix2, found, descrip)
    if (.not. found) exit

    select case (descrip)
    case ('L')
      fmt2 = '(a, ' // trim(fmt2) // ', a)'
      write (this_line, fmt2) this_line(:ix1), l_array(nl+1:nl+nn), trim(this_line(ix2:))
      nl = nl + nn

    case ('I', 'O', 'Z')
      fmt2 = '(a, ' // trim(fmt2) // ', a)'
      write (this_line, fmt2) this_line(:ix1), i_array(ni+1:ni+nn), trim(this_line(ix2:))
      ni = ni + nn

    case ('E', 'G', 'F')
      fmt2 = '(a, ' // trim(fmt2) // ', a)'
      write (this_line, fmt2) this_line(:ix1), r_array(nr+1:nr+nn), trim(this_line(ix2:))
      nr = nr + nn
      cycle

    case default
      write (this_line, '(3a)') this_line(:ix1), '######', trim(this_line(ix2:))

    end select

  enddo

endif

call out_io_line_out (this_line, level)

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine find_format (line, n_prefix, fmt, ix1, ix2, found, descrip)
!
! Subroutine to find the format sub-string.
! This is an internal subroutine not meant for general use.
!
! Input:
!   line -- Character(*): Input string containing the format sub-string.
!   ix0  -- Integer: Index of line to start from. 
!             That is, only line(ix0+1: is searched for a format specifer.
!
! Output:
!   n_prefix -- Integer: Format repeat count
!   fmt      -- Character(*): Extracted format sub-string.
!   ix1      -- Integer: Index in line just before the format sub-string.
!   ix2      -- Integer: Index in line just after the format sub-string.
!   found    -- Logical: Set true if a format is found. False otherwise.
!   descrip  -- Character(1): Optional: Edit descriptor character.
!-

subroutine find_format (line, n_prefix, fmt, ix1, ix2, found, descrip)

implicit none

character(*) line, fmt
character(1), optional :: descrip
character(1) des
integer i, j, ix1, ix2, i1, i2, ix0, n_prefix, n_fmt
logical found

! Init

ix0 = 0
found = .false.
n_prefix = 1
des = ' '

! loop until we find a valid format.

main_loop: do

! Find first \

  i1 = index(line(ix0+1:), '\') + ix0    ! 'find backslash
  if (i1 == ix0 .or. i1 > len(line)-2) return  ! no format found
  ix0 = i1

! Find next \. There must be no spaces inbetween

  i2 = index(line(i1+1:), '\') + i1     ! 'look for second backslash 
  if (i2 == 0) return                   ! no second "\" found
  if (i2 == i1+1) cycle main_loop       ! \\ is not valid

  fmt = line(i1+1:i2-1)
  call str_upcase (fmt, fmt)

  n_fmt = i2 - i1 - 1
  do i = 1, n_fmt
    if (index('0123456789EFGIOLZS.', fmt(i:i)) == 0) cycle main_loop
  enddo

  do i = 1, n_fmt
    if (index('0123456789', fmt(i:i)) == 0) then
      des = fmt(i:i)
      if (i /= 1) read (fmt(1:i-1), *) n_prefix
      exit
    endif
  enddo

  if (index('EFGIOZ', des) == 0) cycle main_loop
  if (i < n_fmt .and. index('0123456789.S', fmt(i+1:i+1)) == 0) cycle main_loop

  do j = i+2, n_fmt
    if (index('0123456789.', fmt(j:j)) == 0) cycle main_loop
  enddo

  ix1 = i1 - 1
  ix2 = i2 + 1
  found = .true.
  if (present(descrip)) descrip = des
  return

enddo main_loop

end subroutine

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! Subroutine header_io (level, routine_name)
!
! Internal routine for output_line4, etc.
!-


Subroutine header_io (level, routine_name)

implicit none

character(*) routine_name
character(20) date_time
integer level

!

call date_and_time_stamp (date_time)

select case (level)

case (s_blank$)

case (s_info$)
  call out_io_line_out('[INFO] ' // trim(routine_name) // ':', level, .false.)

case (s_dinfo$)
  call out_io_line_out('[INFO ' // trim(date_time) // '] ' // trim(routine_name) // ':', level, .false.)

case (s_success$)
  call out_io_line_out('[SUCCESS] ' // trim(routine_name) // ':', level, .false.)

case (s_warn$)
  call out_io_line_out('[WARNING] ' // trim(routine_name) // ':', level, .false.)

case (s_dwarn$)
  call out_io_line_out('[WARNING | ' // trim(date_time) // '] ' // trim(routine_name) // ':', level, .false.)

case (s_error$)
  call out_io_line_out('[ERROR | ' // trim(date_time) // '] ' // trim(routine_name) // ':', level, .false.)

case (s_fatal$)
  call out_io_line_out('[FATAL | ' // trim(date_time) // '] ' // trim(routine_name) // ':', level, .false.)

case (s_abort$)
  call out_io_line_out('[ABORT | ' // trim(date_time) // '] ' // trim(routine_name) // ':', level, .false.)

end select

end subroutine

end module
