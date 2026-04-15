!+
! module tao_preprocess_mod
!
! Transparently runs the Python preprocessor (util/tao_preprocess.py) on
! tao init files before Fortran opens them. The preprocessor expands
! &symbolic_number references into numeric literals, allowing symbolic
! names to appear in real-valued fields (e.g. datum() meas) where the
! Fortran namelist reader cannot accept them.
!
! Lifecycle:
!   - session_init() is called lazily on the first tao_preprocess_file() call.
!     It locates the Python script and creates a PID-keyed temp directory
!     under /tmp (e.g. /tmp/tao_preprocess_12345).
!   - tao_preprocess_cleanup() removes that temp directory. It should be
!     called at tao exit (e.g. from deallocate_everything).
!
! If anything goes wrong (Python not found, env vars not set, write
! errors, etc.) the module silently falls back to returning the original
! file path so behavior is unchanged. A one-time informational message
! is printed when init fails so the user knows preprocessing is unavailable.
!
! Note: All path arguments passed to the shell are single-quoted to avoid
! breakage from spaces or special characters. Paths containing literal
! single quotes are not supported (extremely rare in practice).
!-

module tao_preprocess_mod

use, intrinsic :: iso_c_binding, only: c_int
use, intrinsic :: iso_fortran_env, only: int64
use tao_interface

implicit none

private
public :: tao_preprocess_file
public :: tao_preprocess_cleanup

! C binding for getpid -- used to derive a deterministic temp-dir name.
interface
  function c_getpid() bind(c, name='getpid') result(pid)
    import :: c_int
    integer(c_int) :: pid
  end function c_getpid
end interface

character(400), save :: tmp_dir      = ''
character(400), save :: symbols_file = ''
character(400), save :: script_path  = ''
logical, save        :: init_attempted = .false.
logical, save        :: init_ok        = .false.

contains

!+
! Subroutine tao_preprocess_file (path_in, path_out)
!
! Run the Python preprocessor on path_in and write the result to a
! session-wide temporary directory. Returns the preprocessed path in
! path_out. On any failure, path_out is set to path_in so the caller
! transparently uses the unmodified file.
!
! Before spawning Python, the file is scanned for '&symbolic_number'.
! If the file does not contain that string AND the symbols file is
! empty or absent, the Python call is skipped for performance.
!
! Input:
!   path_in  -- character(*): Path to the original init file.
!
! Output:
!   path_out -- character(*): Path to the preprocessed file, or path_in
!               on any error / when preprocessing is not needed.
!-

subroutine tao_preprocess_file (path_in, path_out)

character(*), intent(in)  :: path_in
character(*), intent(out) :: path_out

character(3000) :: cmd
character(400)  :: base_name
character(2000) :: line
integer :: exit_stat, cmd_stat, ix, iu, ios
integer(int64) :: sym_size
logical :: exists, need_python, has_symbolic

path_out = path_in

! Nothing to do for blank names.
if (len_trim(path_in) == 0) return

! Only preprocess files that actually exist (tao_open_file will emit
! its own "not found" error for anything else).
inquire(file=path_in, exist=exists)
if (.not. exists) return

! Lazy init of the session temp dir / script path.
if (.not. init_attempted) call session_init()
if (.not. init_ok) return

! --- Quick-skip check (issue #10) ---
! Scan the file for '&symbolic_number'. If this file doesn't contain it
! AND the symbols file is empty or absent, there is nothing for Python
! to do, so skip the (expensive) process spawn.
has_symbolic = .false.
iu = lunget()
open(iu, file=path_in, status='old', action='read', iostat=ios)
if (ios == 0) then
  do
    read(iu, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (index(upcase(line), '&SYMBOLIC_NUMBER') /= 0) then
      has_symbolic = .true.
      exit
    endif
  end do
  close(iu)
endif

need_python = has_symbolic
if (.not. need_python) then
  ! Even without &symbolic_number in THIS file, a prior file may have
  ! defined symbols that need substitution here. Check whether the
  ! symbols file exists and has content.
  inquire(file=symbols_file, exist=exists, size=sym_size)
  if (exists .and. sym_size > 2) need_python = .true.  ! '{}' = 2 bytes = empty dict
endif

if (.not. need_python) return

! Build output path: tmp_dir/<basename>.
ix = index(path_in, '/', back=.true.)
if (ix == 0) then
  base_name = path_in
else
  base_name = path_in(ix+1:)
endif
path_out = trim(tmp_dir) // '/' // trim(base_name)

! Invoke the preprocessor with single-quoted paths to prevent shell
! injection and to handle paths with spaces. Redirect stderr to
! /dev/null to keep tao output clean; errors cause fallback below.
cmd = "python3 '" // trim(script_path) // "' '" // trim(path_in) // "' '" // &
      trim(path_out) // "' --symbols-file '" // trim(symbols_file) // &
      "' 2>/dev/null"

call execute_command_line(cmd, wait=.true., exitstat=exit_stat, cmdstat=cmd_stat)

if (cmd_stat /= 0 .or. exit_stat /= 0) then
  path_out = path_in
  return
endif

inquire(file=path_out, exist=exists)
if (.not. exists) then
  path_out = path_in
  return
endif

end subroutine tao_preprocess_file

!+
! Subroutine session_init ()
!
! One-time setup: locate util/tao_preprocess.py and create a PID-keyed
! temp directory (/tmp/tao_preprocess_<pid>). Prints a one-time
! informational message via out_io if setup fails.
!
! Input:  None.
! Output: None. Side effects: sets module save variables (script_path,
!         tmp_dir, symbols_file, init_attempted, init_ok).
!-

subroutine session_init ()

character(400) :: dist_base, mkdir_cmd
character(20)  :: pid_str
character(*), parameter :: r_name = 'tao_preprocess'
integer :: exit_stat, cmd_stat
logical :: exists

init_attempted = .true.

call get_environment_variable('DIST_BASE_DIR', dist_base)
if (len_trim(dist_base) == 0) then
  call get_environment_variable('ACC_ROOT_DIR', dist_base)
endif
if (len_trim(dist_base) == 0) then
  call out_io(s_info$, r_name, 'tao preprocessor not available: DIST_BASE_DIR / ACC_ROOT_DIR not set')
  return
endif

script_path = trim(dist_base) // '/util/tao_preprocess.py'
inquire(file=script_path, exist=exists)
if (.not. exists) then
  call out_io(s_info$, r_name, 'tao preprocessor not available: script not found at ' // trim(script_path))
  return
endif

! Derive temp dir from PID -- deterministic and unique per process.
write(pid_str, '(i0)') int(c_getpid())
tmp_dir = '/tmp/tao_preprocess_' // trim(pid_str)

mkdir_cmd = "mkdir -p '" // trim(tmp_dir) // "' 2>/dev/null"
call execute_command_line(mkdir_cmd, wait=.true., exitstat=exit_stat, cmdstat=cmd_stat)
if (cmd_stat /= 0 .or. exit_stat /= 0) then
  call out_io(s_info$, r_name, 'tao preprocessor not available: could not create temp dir ' // trim(tmp_dir))
  return
endif

symbols_file = trim(tmp_dir) // '/symbols.json'
init_ok = .true.

end subroutine session_init

!+
! Subroutine tao_preprocess_cleanup ()
!
! Remove the session temp directory created by session_init.
! Safe to call even if init was never attempted or failed.
! Should be called at tao exit (e.g. from deallocate_everything).
!
! Input:  None.
! Output: None. Side effects: removes the session temp directory and
!         resets module save variables so a later call can re-init.
!-

subroutine tao_preprocess_cleanup ()

character(500) :: rm_cmd
integer :: exit_stat, cmd_stat

if (len_trim(tmp_dir) == 0) return

rm_cmd = "rm -rf '" // trim(tmp_dir) // "' 2>/dev/null"
call execute_command_line(rm_cmd, wait=.true., exitstat=exit_stat, cmdstat=cmd_stat)

tmp_dir = ''
symbols_file = ''
init_ok = .false.
init_attempted = .false.

end subroutine tao_preprocess_cleanup

end module tao_preprocess_mod
