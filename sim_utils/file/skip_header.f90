!+
! Subroutine SKIP_HEADER (IX_UNIT, ERROR_FLAG)
!
! Subroutine to find the first line of data in a file. The first line of
! data is defined to be the first line whose first non-blank character
! is one of the following set: {+, -, ., 0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
!
! Input:
!     IX_UNIT = unit number of file.
!
! Output:
!     ERROR_FLAG = .false. for normal return
!                = .true. if no data is found
!-

subroutine skip_header (ix_unit, error_flag)

  use precision_def

  implicit none

  integer ix_unit, ix

  character(78) :: line
  character(13) :: set = '+-.0123456789' 

  logical error_flag

!

  do

    read (ix_unit, '(a)', end = 9000) line
    call string_trim(line, line, ix)
    if (ix /= 0 .and. index(set, line(1:1)) /= 0) then
      backspace (unit = ix_unit)
      error_flag = .false.
      return
    endif

  enddo

!

9000  continue
  error_flag = .true.

end subroutine
