!+
! Subroutine tao_show_cmd (what)
!
! Show information on variable, parameters, elements, etc.
!
! Input:
!   what  -- Character(*): What to show.
!-

subroutine tao_show_cmd (what)

use tao_command_mod, dummy => tao_show_cmd

implicit none

type (out_io_output_direct_struct) out_dir_state

integer iu, ix, n, nl, ios
integer :: n_write_file = 0            ! used for indexing 'show write' files

character(*) what
character(200) file_name, fname
character(100) result_id
character(20) switch
character(len(what)) what2
character(n_char_show), allocatable :: lines(:)
character(*), parameter :: r_name = 'tao_show_cmd'

logical opened, err, doprint, err_out, valid

! Init

what2 = what
opened = .false.
doprint = s%com%print_to_terminal
err_out = .true.

! See if the results need to be written to a file.

do
  call tao_next_switch (what2, [character(16):: '-append', '-write', '-noprint', '-no_err_out'], .false., switch, err, ix)
  if (err) return
  if (switch == '') exit

  select case (switch)
  case ('-noprint')
    doprint = .false.

  case ('-no_err_out')
    err_out = .false.

  case ('-append', '-write')
    file_name = what2(:ix)
    call string_trim(what2(ix+1:), what2, ix)

    ix = index(file_name, '*')
    if (ix /= 0) then
      n_write_file = n_write_file + 1
      write (file_name, '(a, i3.3, a)') file_name(1:ix-1), n_write_file, trim(file_name(ix+1:))
    endif

    iu = lunget()
    call fullfilename(file_name, fname, valid)
    if (.not. valid) then
      call out_io (s_error$, r_name, 'NOT A VALID FILE NAME: ' // file_name)
      return
    endif

    if (switch == '-append') then
      open (iu, file = fname, position = 'APPEND', status = 'UNKNOWN', recl = n_char_show, iostat = ios)
    else
      open (iu, file = fname, status = 'REPLACE', recl = n_char_show, iostat = ios)
    endif

    if (ios /= 0) then
      call out_io (s_error$, r_name, 'CANNOT OPEN FILE FOR WRITING: ' // fname)
      return
    endif

    opened = .true.
  end select
end do

! Get results.
! Result_id is for tao_show_this to show exactly what it did.
! This info can be helpful in tao_hook_show_cmd.

call output_direct (get = out_dir_state)
if (opened .and. err_out) call output_direct (iu)  ! tell out_io to write to a file

call tao_show_this (what2, result_id, lines, nl)  
call tao_hook_show_cmd (what2, result_id, lines, nl)

if (opened) call output_direct (iu)  ! tell out_io to write to a file

if (nl > 0) then
  if (result_id == 'ERROR') then
    call out_io (s_error$, r_name, lines(1:nl))
  else
    call output_direct (print_and_capture = doprint)
    call out_io (s_blank$, r_name, lines(1:nl))
  endif
endif

! Finish

call output_direct (set = out_dir_state)

if (opened) then
  close (iu)
  call out_io (s_blank$, r_name, '', 'Written to file: ' // file_name)
endif

end subroutine
