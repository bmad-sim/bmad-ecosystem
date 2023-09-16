!+
! subroutine file_get_open (string, dflt_file_name, file_name, file_unit, readonly)
!
! Routine to query the User for the name of an input file (a la FILE_GET).
! The subroutine will first look at the command line to see if anything was
! entered and will use this as the FILE_NAME. If nothing was entered the
! subroutine will query for a file name by using STRING.
! The subroutine will then open an old file with unit FILE_UNIT.
! The input file will be opened and the unit number returned.
!
! Input:
!   string         -- Character*(*): Query string.
!   dflt_file_name -- Character*(*): Default file name if no name is given.
!   readonly       -- Logical: If true then file will be: readonly, shared.
!
! Output:
!   file_name -- Character: File name supplied by the user.
!   file_unit -- Integer: Unit number for the file.
!
! Example:
!   call file_get_open ('input file:', '[cesr.default].in', file_name, funit, .true.)
!-           

subroutine file_get_open (string, dflt_file_name, file_name, file_unit, readonly)

  use sim_utils_interface, except => file_get_open

  implicit none

  character*(*) file_name, dflt_file_name, string

  integer file_unit, ios

  logical readonly

! get input file name if typed on the same line as the run command

  do 
    call file_get (string, dflt_file_name, file_name)

    file_unit = lunget()

    if (readonly) then
      open (unit = file_unit, file = file_name, status = 'old', &
                                         action = 'read', iostat = ios)
    else
      open (unit = file_unit, file = file_name, status = 'old', iostat = ios)
    endif

    if (ios /= 0) then
      print *, 'ERROR FROM FILE_GET: ERROR OPENING: ', trim(file_name)
      print *, '      TRY AGAIN ...'
      cycle
    endif

    return

  enddo

end subroutine
