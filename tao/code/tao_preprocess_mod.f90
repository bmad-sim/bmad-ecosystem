!+
! module tao_preprocess_mod
!
! Transparently runs the Python preprocessor (util/tao_preprocess.py) on
! tao init files before Fortran opens them. The preprocessor expands
! &symbolic_number references into numeric literals, allowing symbolic
! names to appear in real-valued fields (e.g. datum() meas) where the
! Fortran namelist reader cannot accept them.
!
! If anything goes wrong (Python not found, env vars not set, write
! errors, etc.) the module silently falls back to returning the original
! file path so behavior is unchanged.
!-

module tao_preprocess_mod

use sim_utils

implicit none

private
public :: tao_preprocess_file

character(400), save :: tmp_dir = ''
character(400), save :: symbols_file = ''
character(400), save :: script_path = ''
logical, save :: init_attempted = .false.
logical, save :: init_ok = .false.

contains

!+
! Subroutine tao_preprocess_file (path_in, path_out)
!
! Run the Python preprocessor on path_in and write the result to a
! session-wide temporary directory. Returns the preprocessed path in
! path_out. On any failure, path_out is set to path_in so the caller
! transparently uses the unmodified file.
!-

subroutine tao_preprocess_file (path_in, path_out)

character(*), intent(in) :: path_in
character(*), intent(out) :: path_out

character(1000) :: cmd
character(400) :: base_name, dir_name
integer :: exit_stat, cmd_stat, ix
logical :: exists

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

! Build output path: tmp_dir/<basename>.
ix = index(path_in, '/', back=.true.)
if (ix == 0) then
  base_name = path_in
else
  base_name = path_in(ix+1:)
endif
path_out = trim(tmp_dir) // '/' // trim(base_name)

! Invoke the preprocessor. Redirect stderr to /dev/null to keep tao
! output clean; errors cause fallback to the original file below.
cmd = 'python3 ' // trim(script_path) // ' ' // trim(path_in) // ' ' // &
      trim(path_out) // ' --symbols-file ' // trim(symbols_file) // &
      ' 2>/dev/null'

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
! One-time setup: locate the util/tao_preprocess.py script and create a
! session-wide temp directory. Sets init_ok = .true. on success.
!-

subroutine session_init ()

character(400) :: dist_base
integer :: exit_stat, cmd_stat, iu, ios
logical :: exists

init_attempted = .true.

call get_environment_variable('DIST_BASE_DIR', dist_base)
if (len_trim(dist_base) == 0) then
  call get_environment_variable('ACC_ROOT_DIR', dist_base)
endif
if (len_trim(dist_base) == 0) return

script_path = trim(dist_base) // '/util/tao_preprocess.py'
inquire(file=script_path, exist=exists)
if (.not. exists) return

! Use mktemp -d to create a unique per-session directory. Write the path
! to a small file we can read back from Fortran.
iu = lunget()
open(iu, status='scratch')  ! just to reserve it; won't actually use
close(iu)

block
  character(400) :: path_file
  character(500) :: make_cmd
  call get_unique_name(path_file)
  make_cmd = 'mktemp -d /tmp/tao_preprocess.XXXXXX > ' // trim(path_file) // ' 2>/dev/null'
  call execute_command_line(make_cmd, wait=.true., exitstat=exit_stat, cmdstat=cmd_stat)
  if (cmd_stat /= 0 .or. exit_stat /= 0) return

  iu = lunget()
  open(iu, file=path_file, status='old', action='read', iostat=ios)
  if (ios /= 0) return
  read(iu, '(a)', iostat=ios) tmp_dir
  close(iu, status='delete')
  if (ios /= 0 .or. len_trim(tmp_dir) == 0) return
end block

symbols_file = trim(tmp_dir) // '/symbols.json'
init_ok = .true.

end subroutine session_init

!+
! Produce a unique filesystem path in /tmp to use as a short-lived
! scratch file for capturing shell output.
!-

subroutine get_unique_name (path)
character(*), intent(out) :: path
integer :: pid
integer :: vals(8)
real :: rnd
call date_and_time(values=vals)
call random_number(rnd)
write(path, '(a, i0, a, i0)') '/tmp/tao_pp_init_', vals(7)*1000+vals(8), '_', int(rnd*1000000)
end subroutine get_unique_name

end module tao_preprocess_mod
