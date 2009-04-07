#include "CESR_platform.inc"

module filename_mod

contains

!+
! Subroutine FullFileName(filename, outfile, valid) 
!
! Description:  
!   Returns the full filename with a leading environment variable expanded
!   into the proper form (full path on Unix, logical on VMS).  NOTE: It is
!   the responsibility of the calling routine to provide a large enough
!   string to hold a fully expanded file name on return.
!
! Acknowledgments:
!   This routine draws on the routine trnfnm.F in the CLEO libraries.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   filename -- Character(*): Name requiring expansion of an environment
!                  variable.  The name can be of the form:
!                      LOGICAL:foo.bar        (VMS form)
!                      $ENV_VARIABLE/foo.bar  (UNIX form) 
!                  Hardwired paths are also accepted but strongly discouraged 
!                  because of the resulting platform dependence (see above).
! Output:
!   outfile  -- Character(*): Expanded name.
!   valid    -- Logical, optional: Set False if, under UNIX, an environmental 
!                 variable in the filename does not exist. Set True otherwise.
!
! Examples:  
!    Assume we have a variable DUMMY which is defined as an environment
!    variable on Unix and a logical on VMS.
!      On Unix:  DUMMY = /home/cesr/dummy
!      On VMS:   DUMMY = [cesr.dummy]
!    Then the following behaviors will result when using fullfilename:
!    (Note: Logicals are automatically interpreted by the system under VMS) 
!
!      System    filename                    outfile
!      ------   ------------------------    ---------------------------
!      UNIX     'DUMMY:foo.bar'             '/home/cesr/dummy/foo.bar'
!      UNIX     'DUMMY:[a.b]foo.bar'        '/home/cesr/dummy/a/b/foo.bar'
!      UNIX     '$DUMMY/foo.bar'            '/home/cesr/dummy/foo.bar'
!      UNIX     '/home/cesr/dummy/foo.bar'  '/home/cesr/dummy/foo.bar'
!      UNIX     '[cesr.dummy]foo.bar'        NOT a valid Unix file name.
!
!      VMS      'DUMMY:foo.bar'             'DUMMY:foo.bar'
!      VMS      '$DUMMY/foo.bar'            'DUMMY:foo.bar'
!      VMS       $DUMMY/../a/foo.bar'       'DUMMY:[-.a]foo.bar'
!      VMS      '/home/cesr/dummy/foo.bar'  '[home.cesr.dummy]foo.bar
!      VMS      'sub/foo.bar'               '[.sub]foo.bar'
!      VMS      './sub/foo.bar'             '[.sub]foo.bar'
!      VMS      '../sub/foo.bar'            '[-.sub]foo.bar'
!      VMS      '[cesr.dummy]foo.bar'       '[cesr.dummy]foo.bar' 
!      VMS      './'                        '[]' 
! 
! Author     :  M. Palmer   9/20/01
!-

subroutine FullFileName (filename, outfile, valid)

  implicit none

! Argument Declarations         

  character(*) filename
  character(*) outfile
  logical, optional :: valid

  character(16) :: r_name = 'FullFileName'
  character(len(outfile)) ExpandName        ! Expanded Name

  Integer InLen, iDollar, iColon, iSlash, iLeftB, iRightB

#if defined(CESR_VMS)
  Integer i, ix

#else
  Integer i
  External GetEnv
  integer explen
  logical, save :: rcsini = .true.
#endif

! Assign input to output by default

  if (present(valid)) valid = .false.
  outfile = FileName

!     Get length of input file name

  InLen = Len_Trim(FileName)

!     

  If (InLen .gt. len(outfile)) then
    Write(6, '(a12, a, i4, a)') r_name, &
        '[Warning]:  Filename too large (', len(outfile), ' char. limit):'
    Write(6, '(a)') FileName
    Return
  Endif

! Locate special characters (first dollar-sign, first colon, first slash)

  iDollar = Index(FileName(:InLen), '$' )
  iColon  = Index(FileName(:InLen), ':' )
  iSlash  = Index(FileName(:InLen), '/' )

  ExpandName = FileName

!-----------------------------------------------------------------------
! Translation on UNIX

#if defined(CESR_UNIX)

! A UNIX-style environment variable will have a leading '$' 

  If (iDollar == 1) then
        
! Expand Unix-style environment variable names
! Environment variable specifies the full name

    If (iSlash == 0) then
          
      Call GetEnv(Filename(2:InLen), ExpandName)
      if (len_trim(ExpandName) == 0) return

! Environment variable plus short name specifies the full name

    Elseif (iSlash > 2) then
      Call GetEnv(Filename(2:iSlash-1), ExpandName)
      ExpLen = Len_Trim(ExpandName)
      if (ExpLen == 0) return
      ExpandName(ExpLen+1:) = FileName(iSlash:InLen)
    Endif

! VMS-style logicals will have a trailing colon

  else

    if (iColon > 1) then
      Call GetEnv(Filename(1:iColon-1), ExpandName)
      ExpLen = Len_Trim(ExpandName)
      if (ExpLen == 0) return
      ExpandName(ExpLen+1:) = '/' // FileName(iColon+1:InLen)
    endif

    
    iLeftB  = index(ExpandName, '[')
    iRightB = index(ExpandName, ']')    

    if (iLeftB == 0) then
      if (iRightB /= 0) return
    else
      if (iRightB < iLeftB) return
      do i = iLeftB+1, iRightB
        if (ExpandName(i:i) == '.') ExpandName(i:i) = '/'
        if (ExpandName(i:i) == ']') ExpandName(i:i) = '/'
      enddo
      ExpandName = ExpandName(1:iLeftB-1) // ExpandName(iLeftB+1:)
    endif

  Endif

  outfile = ExpandName
  if (present(valid)) valid = .true.

!-----------------------------------------------------------------------
! Translation on VMS

#elif defined(CESR_VMS)
      
!   '$DUMMY/foo.bar'     --> 'DUMMY:foo.bar'
!   '$DUMMY/a/b/foo.bar' --> 'DUMMY:[a.b]foo.bar'

  if (iDollar == 1) then  ! has environment variable
    ExpandName = ExpandName(2:InLen)
    if (iColon /= 0) return  ! Error
    ix = Index (ExpandName, '/' )  
    if (ix /= 0) then
      ! If there are multiple "/" then need "[...]" construct.
      if (index (ExpandName(ix+1:), '/') == 0) then
        ExpandName = ExpandName(:ix-1) // ':' // ExpandName(ix+1:)
      else
        call dir_unix_to_vms (ExpandName(ix:))
        ExpandName = ExpandName(:ix-1) // ':' // ExpandName(ix:)
      endif
    endif

! Translation from UXIX to VMS for names that have a slash but do not 
! look like they have a logical directory. 

  elseif (iColon == 0 .and. iSlash /= 0 .and. index(ExpandName, ']') == 0) then
    call dir_unix_to_vms (ExpandName)
  endif

  outfile = ExpandName
  if (present(valid)) valid = .true.

!-----------------------------------------------------------------------
! Translation on WINDOWS (from unix to windows)

#elif defined(CESR_WINCVF)

! A UNIX-style environment variable will have a leading '$' 

  If ( iDollar == 1 ) then
        
! Expand Unix-style environment variable names
! Environment variable specifies the full name

    If (iSlash == 0) then
          
      Call GetEnv(Filename(2:InLen), ExpandName)
      if (Len_Trim(ExpandName) == 0) return

! Environment variable plus short name specifies the full name

    Elseif (iSlash > 2) then
      Call GetEnv(Filename(2:iSlash-1), ExpandName)
      ExpLen = Len_Trim(ExpandName)
      if (ExpLen == 0) return
      ExpandName(ExpLen+1:) = FileName(iSlash:InLen)
    Endif

! VMS-style logicals will have a trailing colon

  Elseif (iColon > 1) then
    Call GetEnv(Filename(1:iColon-1), ExpandName)
    ExpLen      = Len_Trim(ExpandName)
    if (ExpLen == 0) return
    ExpandName(ExpLen+1:) = '/' // FileName(iColon+1:InLen)
    outfile = ExpandName
  Endif
  
  
  do i=1, Len_Trim(ExpandName)
     If (ExpandName(i:i)=='/') ExpandName(i:i)='\' ! '
  end do

  outfile = ExpandName
  if (present(valid)) valid = .true.

#else

!-----------------------------------------------------------------------
! Translation on other platforms.  NONE PRESENTLY DEFINED

#endif

!-----------------------------------------------------------------------
contains

subroutine dir_unix_to_vms (outfile)

character(*) outfile

integer i, ix, isl, isl2

! Special case

if (outfile == './') then
  outfile = '[]'
  return
endif

!

do
  if (outfile(1:2) /= './') exit
  outfile = outfile(3:)
enddo

isl = index (outfile, '/')

if (isl == 0) return
    
do 
  i = index(outfile, '..')
  if (i == 0) exit
  outfile = outfile(:i-1) // '-' // outfile(i+2:)
enddo

if (outfile(1:1) == '/') then
  outfile = '[' // outfile(2:)
elseif (outfile(1:1) == '-') then
  outfile = '[' // outfile
else
  outfile = '[.' // outfile
endif

i = index(outfile, '/')
if (i /= 0) then
  do
    i = index(outfile, '/')
    if (i == 0) exit
    outfile = outfile(:i-1) // '.' // outfile(i+1:)
    ix = i
  enddo
  outfile = outfile(:ix-1) // ']' // outfile(ix+1:)
endif

end subroutine

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function SplitFileName(filename, path, basename, is_relative) result (ix_char)
!
! Routine to take filename and splits it into its constituent parts, 
! the directory path and the base file name.  The return 
! value, ix_char, is the length of the path variable. If the 
! specified file is local, the return value is 0.
!
! Module needed:
!   use sim_utils
!
! Input:  
!   FileName    -- Character(*): Filename to be split.
!
! Output:
!   Path        -- Character(*): Path to file with terminating /, ], or :
!   BaseName    -- Character(*): Base filename, no path.
!   is_relative -- Logical, optional: True if pathname is relative. False otherwise.
!   ix_char     -- Integer: Number of characters in path string.
!
! Authors:  
!   M. Palmer   2001/9/23
!   D. Sagan    2006/8/23     Added is_relative      
!-

function SplitFileName(FileName, Path, BaseName, is_relative) result (ix_char)

  use output_mod

  implicit none

  character(*) FileName, Path, BaseName
  logical, optional :: is_relative

  character(16), parameter :: r_name = 'SplitFileName'
  character(62), parameter :: AtoZ = &
        'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890'

  Integer InLen, ix_char, ix
  Integer iBracket, iColon, iSlash

! Get length of input file name

  InLen = Len_Trim(FileName)

! Locate special characters (last dollar-sign, last colon, last slash)

  ix   = 0
  iBracket = Scan( FileName(:InLen), ']', .True. )
  If (iBracket .gt. ix) ix = iBracket
  iColon   = Scan( FileName(:InLen), ':', .True. )
  If (iColon .gt. ix) ix = iColon
  iSlash   = Scan( FileName(:InLen), '/', .True. )
  If (iSlash .gt. ix) ix = iSlash

      
  Path     = FileName(:ix)
  BaseName = FileName(ix+1:InLen)
  call string_trim(path, path, ix)
  ix_char = ix

! find if path name is relative or absolute.

  if (.not. present(is_relative)) return

  if (path(1:1) == ' ') then
    is_relative = .true.
  elseif (path(1:1) == '/') then
    is_relative = .false.
  elseif (path(1:1) == '~') then
    is_relative = .false.
  elseif (index(path, ':') /= 0) then
    is_relative = .false.
  elseif (path(1:1) == '.') then
    is_relative = .true.
  elseif (index(AtoZ, path(1:1)) /= 0) then
    is_relative = .true.
  elseif (path(1:2) == '[.') then
    is_relative = .true.
  elseif (path(1:2) == '[-') then
    is_relative = .true.
  elseif (path(1:1) == '[' .and. index(AtoZ, path(2:2)) /= 0) then
    is_relative = .false.
  elseif (path(1:1) == '$') then
    is_relative = .false.
  else
    is_relative = .false.
    call out_io (s_warn$, r_name, &
                    'CANNOT TELL IF PATH IS ABSOLUTE OR RELATIVE: ' // path, &
                    'ASSUMING ABSOLUTE. PLEASE CONTACT DCS! ')
  endif
      
End function

!--------------------------------------------------------------------
!--------------------------------------------------------------------

!-
! Subroutine simplify_path (name_in, name_out)
!
! Routine to remove './' and '/dir/..' constructs.
!
! Input:
!   name_in -- Character(*): Input path.
!
! Output:
!   name_out -- Character(*): Simplified path
!-

subroutine simplify_path (name_in, name_out)

implicit none

character(*) name_in, name_out
integer i, ix, ix_last

! First get rid of any './' constructs.

name_out = name_in

ix_last = 1
do
  ix = index(name_out(ix_last:), './')
  if (ix == 0) exit
  ix = ix + ix_last - 1
  if (ix == 1) then  ! beginning of string
    name_out = name_out(3:)
    cycle
  endif
  if (name_out(ix-1:ix-1) == '/') then
    name_out = name_out(:ix-1) // name_out(ix+2:)
  else
    ix_last = ix + 1
  endif
enddo


! Second get rid of any '/dir/..' constructs.

ix_last = 1
loop: do 
  ix = index(name_out(ix_last:), '/..')
  if (ix == 0) return
  ix = ix + ix_last - 1
  if (name_out(ix+3:ix+3) /= '/') then
    ix_last = ix + 3
    cycle
  endif
  do i = ix-1, 1, -1
    if (name_out(i:i) == '/') then
      name_out = name_out(:i-1) // name_out(ix+3:)
      cycle loop
    endif
  enddo
  name_out = name_out(ix+3:)
enddo loop

end subroutine

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Subroutine append_subdirectory (dir, sub_dir, dir_out)
!
! Subroutine to combine a directory specification with a 
! subdirectory specification to form a complete directory
! string.
!
! Examples:
!     System      dir        sub_dir     dir_out
!     UNIX       "abc"      "def"       "abc/def"
!     UNIX       "abc"      "/def"      "abc/def"
!     UNIX       "abc/"     "def"       "abc/def"
!     UNIX       "abc/"     "/def"      "abc/def"
!     VMS        "[abc]"    "[def]"     ""  (Error)
!     VMS        "[abc]"    "[.def]"    "[abc.def]"
!     VMS        "[abc]"    "def"       "[abc]def"
!     VMS        ""         "[.def]"    "[.def]"
!
! Modules needed:
!   use filename_mod
!
! Input:
!   dir     -- Character(*): Head directory.
!   sub_dir -- Character(*): Subdirectory.
!
! Output:
!   dir_out -- Character(*): Complete directory spec.
!   err     -- Logical: Set True if there is an error. False other wise
!-

subroutine append_subdirectory (dir, sub_dir, dir_out, err)

use output_mod

implicit none

character(*) dir, sub_dir, dir_out
character(len(dir_out)) temp
character(24) :: r_name = 'AppendSubDirectory'

integer n_dir

logical err

! Easy cases

err = .false.

if (dir == "") then
  dir_out = sub_dir
  return
elseif (sub_dir == "") then
  dir_out = dir
  return
endif

err = .true.

! VMS version

#if defined(CESR_VMS)

n_dir = len_trim(dir)

if (dir(n_dir:n_dir) /= ']') then
  call out_io (s_fatal$, r_name, 'BAD DIRECTORY STRUCTURE: ' // dir)
  dir_out = 'XXX'
  return
endif

if (sub_dir(1:1) == '[') then
  if (sub_dir(2:2) /= '.') then
    call out_io (s_fatal$, r_name, 'BAD SUB-DIRECTORY STRUCTURE: ' // sub_dir)
    dir_out = 'XXX'
    return
  endif
  temp = dir(:n_dir-1) // sub_dir(2:)
else
  temp = dir(:n_dir) // sub_dir(:)
endif

dir_out = temp

! non-VMS version

#else

n_dir = len_trim(dir)

if (dir(n_dir:n_dir) == '/' .and. sub_dir(1:1) == '/') then
  temp = dir(1:n_dir) // sub_dir(2:)
elseif ((dir(n_dir:n_dir) == '/' .and. sub_dir(1:1) /= '/') .or. &
        (dir(n_dir:n_dir) /= '/' .and. sub_dir(1:1) == '/')) then
  temp = dir(1:n_dir) // sub_dir
else
  temp = dir(1:n_dir) // '/' // sub_dir
endif

dir_out = temp

#endif

err = .false.

end subroutine

end module
