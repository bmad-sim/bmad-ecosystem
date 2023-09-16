!+
! Subroutine FILE_GET (STRING, DFLT_FILE_NAME, FILE_NAME)
!
! Routine to query the User for the name of a file.
! The subroutine will first look at the command line to see if anything was
! entered and will use this as the FILE_NAME. If nothing was entered the
! subroutine will query for a file name by using STRING.
!
! Input:
!   STRING         -- Character*(*): Query string.
!   DFLT_FILE_NAME -- Character*(*): Default file name if no name is given.
!                  -- There are 3 parts to the default which are treated
!                     separately: 
!                            The directory, The file name, The file suffix
!
! Output:
!     FILE_NAME -- Character: File name supplied by the user.
!
! Example: (This uses a default directory and suffix)
!     call ('Input file:', '[CESR.DEFAULTS].IN', file_name)
!-           

subroutine file_get (string, dflt_file_name, file_name)

use precision_def

implicit none

character(*) file_name, dflt_file_name, string
character dflt_dir*60, dflt_name*60, dflt_suffix*60

integer nret, ix, iargc

logical :: try_foreign = .true.

! parse default

ix = index (dflt_file_name, ']')
if (ix /= 0) then
  dflt_dir = dflt_file_name(:ix)
  dflt_name = dflt_file_name(ix+1:)
else
  dflt_dir = ' '
  dflt_name = dflt_file_name
endif

ix = index (dflt_name, '.')
if (ix /= 0) then
  dflt_suffix = dflt_name(ix:)
  dflt_name = dflt_name(:ix-1)
else
  dflt_suffix = ' '
endif

! get input file name if typed on the same line as the run command.
! logical TRY_FOREIGN is used so if FILE_GET is called a second time
! we don't try to use the command line twice

if (try_foreign) then
  nret = iargc()
  if (nret > 0) call getarg(1, file_name)
else
  nret = 0
endif

! if no name typed then issue a query at the terminal

if (nret == 0) then

  if (len_trim(dflt_dir) /= 0) print *, 'Note: Default directory is: ', &
                                                              trim(dflt_dir)

  if (len_trim(dflt_name) /= 0) then
    write (*, '(1x, 4a, 1x)', advance = 'NO') &
                      trim(string), ' <CR = ', trim(dflt_name), '> '
    read '(a)', file_name
    if (len_trim(file_name) == 0) file_name = dflt_file_name
  else
    do
      write (*, '(1x, 2a)', advance = 'NO') trim(string), ' '
      read '(a)', file_name
      if (len_trim(file_name) == 0) then
        print *, 'Please enter a file name...'
      else
        exit
      endif
    enddo
  endif
endif

call file_suffixer (file_name, file_name, dflt_suffix, .false.)

if (len_trim(dflt_dir) /= 0) &
          call file_directorizer (file_name, file_name, dflt_dir, .false.)

end subroutine
