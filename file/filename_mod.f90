module filename_mod

use output_mod
use sim_utils_interface

contains

!+
! Subroutine FullFileName(filename, outfile, valid) 
!
! Description:  
!   Returns the full filename with a leading environment variable expanded
!   into the proper form (full path on Unix).  NOTE: It is
!   the responsibility of the calling routine to provide a large enough
!   string to hold a fully expanded file name on return.
!
! Also see: SplitFileName
!
! Input:
!   filename -- Character(*): Name requiring expansion of an environment
!                  variable.  The name can be of the form:
!                      $ENV_VARIABLE/foo.bar  (UNIX form)
!
! Output:
!   outfile  -- Character(*): Expanded name.
!   valid    -- Logical, optional: Set False if, under UNIX, an environmental 
!                 variable in the filename does not exist. Set True otherwise.
!
! Examples:  
!    Assume we have a variable DUMMY which is defined as an environment
!    variable
!      DUMMY = /home/cesr/dummy
!    Then the following behaviors will result when using fullfilename:
!
!      Filename                     outfile
!      -----------------------      ---------------------------
!      '$DUMMY/foo.bar'             '/home/cesr/dummy/foo.bar'
!      '/home/cesr/dummy/foo.bar'   '/home/cesr/dummy/foo.bar'
!-

subroutine FullFileName (filename, outfile, valid)

implicit none

! Argument Declarations         

character(*) filename
character(*) outfile
logical, optional :: valid

character(*), parameter :: r_name = 'FullFileName'
character(len(outfile)) expandname, name0        ! Expanded Name

integer InLen, iDollar, iColon, iSlash
integer i, ix
integer explen

! Get length of input file name

if (present(valid)) valid = .false.

InLen = len_trim(filename)

if (InLen > len(outfile)) then
  call out_io (s_error$, r_name, 'Outfile argument string length too small: \i0\ ', len(outfile))
  return
endif

! Assign input to output by default

outfile = unquote(filename)
if (outfile == '') return  ! Blank is invalid
name0 = outfile
expandname = outfile

! Locate special characters (first dollar-sign, first colon, first slash)

iDollar = Index(name0(:InLen), '$' )
iColon  = Index(name0(:InLen), ':' )
iSlash  = Index(name0(:InLen), '/' )


!-----------------------------------------------------------------------
! Translation on WINDOWS (from unix to windows)

#if defined(CESR_WINCVF)

! A UNIX-style environment variable will have a leading '$' 

if (iDollar == 1) then
      
! Expand Unix-style environment variable names
! Environment variable specifies the full name

  if (iSlash == 0) then
        
    call GetEnv(filename(2:InLen), expandname)
    if (Len_Trim(expandname) == 0) return

! Environment variable plus short name specifies the full name

  elseif (iSlash > 2) then
    call GetEnv(filename(2:iSlash-1), expandname)
    ExpLen = Len_Trim(expandname)
    if (ExpLen == 0) return
    expandname(ExpLen+1:) = filename(iSlash:InLen)
  endif

! VMS-style logicals will have a trailing colon

elseif (iColon > 1) then
  call GetEnv(filename(1:iColon-1), expandname)
  ExpLen      = Len_Trim(expandname)
  if (ExpLen == 0) return
  expandname(ExpLen+1:) = '/' // filename(iColon+1:InLen)
endif


do i=1, Len_Trim(expandname)
   if (expandname(i:i)=='/') expandname(i:i)='\' ! '
end do

outfile = expandname
if (present(valid)) valid = .true.

!-----------------------------------------------------------------------
! Translation on UNIX

#else

! Tilde

if (name0(1:1) == '~') then
  call GetEnv('HOME', expandname)
  if (expandname == '') return

  if (name0(2:2) == '/') then
    expandname = trim(expandname) // name0(2:)
  else
    ix = index(expandname, '/', back = .true.)
    if (ix == 0) then
      expandname = trim(expandname) // '/' // name0
    else
      expandname = expandname(1:ix) // name0
    endif
  endif

! A UNIX-style environment variable will have a leading '$' 

elseif (iDollar == 1) then
  if (iSlash == 0) then
        
    call GetEnv(name0(2:InLen), expandname)
    if (len_trim(expandname) == 0) return

  ! Environment variable plus short name specifies the full name

  elseif (iSlash > 2) then
    call GetEnv(name0(2:iSlash-1), expandname)
    ExpLen = Len_Trim(expandname)
    if (ExpLen == 0) return
    expandname(ExpLen+1:) = name0(iSlash:InLen)
  endif
endif

outfile = expandname
if (present(valid)) valid = .true.

#endif

end subroutine FullFileName

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
!+
! Function Splitfilename(filename, path, basename, is_relative) result (ix_char)
!
! Routine to take filename and splits it into its constituent parts, 
! the directory path and the base file name.  The return 
! value, ix_char, is the length of the path variable. if the 
! specified file is local, the return value is 0.
!
! Also see: Fullfilename
!
! Input:  
!   filename    -- Character(*): Filename to be split.
!
! Output:
!   path        -- Character(*): Path to file with terminating /, ], or :
!   basename    -- Character(*): Base filename, no path.
!   is_relative -- Logical, optional: True if pathname is relative. False otherwise.
!   ix_char     -- Integer: Number of characters in path string.
!-

function Splitfilename(filename, path, basename, is_relative) result (ix_char)

implicit none

character(*) filename, path, basename
logical, optional :: is_relative

character(*), parameter :: r_name = 'Splitfilename'

Integer InLen, ix_char, ix
Integer iBracket, iColon, iSlash

! Get length of input file name

InLen = Len_Trim(filename)

! Locate special characters (last dollar-sign, last colon, last slash)

ix = 0

#if defined(CESR_WINCVF)
iBracket = Scan(filename(:InLen), ']', .True.)
if (iBracket .gt. ix) ix = iBracket
iColon   = Scan(filename(:InLen), ':', .True.)
if (iColon .gt. ix) ix = iColon
#endif

iSlash   = Scan(filename(:InLen), '/', .True.)
if (iSlash .gt. ix) ix = iSlash

path     = filename(:ix)
basename = filename(ix+1:InLen)
call string_trim(path, path, ix)
ix_char = ix

! find if path name is relative or absolute.

if (present(is_relative)) is_relative = file_name_is_relative(path)

end function Splitfilename

!-----------------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!+
! Function file_name_is_relative (file_name) result (is_rel)
! 
! Routine to determine if a file name is absolute or relative.
! Examples:
!     abc       ! Relative
!     $ENV/abc  ! Absolute
!     /nfs/opt  ! Absolute
!
! Note: This routine works for Unix and has not been extended to Windows.
!
! Input:
!   file_name -- Character(*): Name of file or path.
!
! Output:
!   is_rel    -- Logical: True if relative, False otherwise
!-

function file_name_is_relative (file_name) result (is_rel)

implicit none

character(*) file_name
character(len(file_name)) name
character(*), parameter :: r_name = 'file_name_is_relative'
integer ix
logical is_rel

!

call string_trim(file_name, name, ix)

if (name(1:1) == ' ') then
  is_rel = .true.
elseif (name(1:1) == '/') then
  is_rel = .false.
elseif (name(1:1) == '.') then
  is_rel = .true.
#if defined(CESR_WINCVF)
elseif (name(1:1) == '\') then         !'
  is_rel = .false.
elseif (index(name, ':\') /= 0) then   !'
  is_rel = .false.
elseif (index(name, ':/') /= 0) then   !'
  is_rel = .false.
#else
elseif (name(1:1) == '$') then
  is_rel = .false.
elseif (name(1:1) == '~') then
  is_rel = .false.
#endif
! can have file names like "#abc#" so take everything else to be relative.
else
  is_rel = .true.
endif

end function file_name_is_relative 
      
!-----------------------------------------------------------------------------
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

  if (name_out(ix+3:ix+3) /= '/' .and. name_out(ix+3:ix+3) /= ' ') then
    ix_last = ix + 3
    cycle
  endif

  if (ix > 2) then
    if (name_out(ix-2:ix-1) == '..') then
      ix_last = ix + 3
      cycle
    endif
  endif

  do i = ix-1, 1, -1
    if (name_out(i:i) == '/') then
      name_out = name_out(:i) // name_out(ix+4:)
      cycle loop
    endif
  enddo
  name_out = name_out(ix+4:)
enddo loop

end subroutine simplify_path

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

implicit none

character(*) dir, sub_dir, dir_out
character(len(dir_out)) temp
character(*), parameter :: r_name = 'Append_SubDirectory'

integer n_dir

logical err

! Easy cases

err = .false.

if (dir == "" .or. dir == './' .or. dir == '.') then
  dir_out = sub_dir
  return
elseif (sub_dir == "") then
  dir_out = dir
  return
endif

err = .true.

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

err = .false.

end subroutine append_subdirectory

end module
