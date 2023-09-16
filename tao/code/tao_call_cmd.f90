!+
! Subroutine tao_call_cmd (file_name, cmd_arg)
!
! Routine to open a tao command file.
! 
! Input:
!   file_name  -- Character(*): Name of the tao command file.
!   cmd_arg(9) -- Character(*), optional: Command file arguments.
!
! Output:
!   s%global%cmd_file -- Integer: Logical unit number of the 
!                                   command file.
!-

subroutine tao_call_cmd (file_name, cmd_arg)

use tao_interface, dummy => tao_call_cmd

implicit none

type (tao_command_file_struct), allocatable :: tmp_cmd_file(:)

character(*) file_name
character(*), optional :: cmd_arg(:)
character(200) full_name, basename
character(*), parameter :: r_name = 'tao_call_cmd'

integer iu, nl, nl0, ix
logical err, is_relative, valid, reset_at_end

! Call options

reset_at_end = .false.

if (file_name(1:2) == '-p') then
  if (index('-ptc', trim(file_name)) /= 1 .or. cmd_arg(1) == '' .or. cmd_arg(2) /= '') then
    call out_io (s_error$, r_name, 'MALFORMED PTC COMMAND FILE CALL')
    return
  endif

  call read_ptc_command77 (cmd_arg(1))

elseif (file_name(1:2) == '-n') then
  if (index('-no_calc', trim(file_name)) /= 1 .or. .not. present(cmd_arg)) then
    call out_io (s_error$, r_name, 'MALFORMED "-no_calc" OPTION.')
    return
  endif

  reset_at_end = .true.
  file_name = cmd_arg(1)
  cmd_arg = [cmd_arg(2:), cmd_arg(1)]
  cmd_arg(size(cmd_arg)) = ''
endif

! 

! Open the command file and store the unit number

nl0 = s%com%cmd_file_level
nl = nl0 + 1
s%com%cmd_file_level = nl

! reallocate cmd_file array

if (ubound(s%com%cmd_file, 1) < nl) then
  call move_alloc (s%com%cmd_file, tmp_cmd_file)
  allocate (s%com%cmd_file(0:nl))
  s%com%cmd_file(0:nl0) = tmp_cmd_file
endif

!

call fullfilename (file_name, full_name, valid)
ix = splitfilename (full_name, s%com%cmd_file(nl)%dir, basename, is_relative)
if (is_relative) then
  call append_subdirectory (s%com%cmd_file(nl0)%dir, s%com%cmd_file(nl)%dir, s%com%cmd_file(nl)%dir, err)
  if (err) then
    s%com%cmd_file_level = nl0
    call tao_abort_command_file()
    return
  endif
  call append_subdirectory (s%com%cmd_file(nl0)%dir, full_name, full_name, err)
endif

iu = lunget()
call tao_open_file (full_name, iu, full_name, s_error$)
if (iu == 0) then ! open failed
  s%com%cmd_file_level = nl0
  call tao_abort_command_file()
  return
endif

s%com%cmd_file(nl)%ix_unit = iu
s%com%cmd_file(nl)%full_name = full_name

! Save command arguments.

if(present(cmd_arg)) then
  s%com%cmd_file(nl)%cmd_arg = cmd_arg
else
  s%com%cmd_file(nl)%cmd_arg = ' '
endif

s%com%cmd_file(nl)%reset_at_end      = reset_at_end
s%com%cmd_file(nl)%lattice_calc_save = s%global%lattice_calc_on
s%com%cmd_file(nl)%plot_save         = s%global%plot_on

end subroutine 
