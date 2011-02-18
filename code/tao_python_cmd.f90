!+
! Subroutine tao_python_cmd (input_str)
!
! Print information in a form easily parsed by a scripting program like python.
!
! Input:
!   input_str  -- Character(*): What to show.
!-

subroutine tao_python_cmd (input_str)

use tao_mod

implicit none

character(n_char_show), allocatable, save :: lines(:)
type (tao_plot_struct), pointer :: p
type (tao_plot_array_struct), allocatable, save :: plot(:)

character(*) input_str
character(200) line
character(20) cmd, command
character(20) :: r_name = 'tao_python_cmd'
character(20) :: cmd_names(5) = &
          ['visible_plots ', 'template_plots', 'graph         ', 'curve      ', &
           'plot          ']

integer i, j, ix, ix_cmd, nl

logical err, print_flag

!

call string_trim(input_str, line, ix)
cmd = line(1:ix)
call string_trim(line(ix+1:), line, ix)

call match_word (cmd, cmd_names, ix, matched_name = command)
if (ix == 0) then
  print *, '***PYTHON WHAT? WORD NOT RECOGNIZED: ' // command
  return
endif

if (ix < 0) then
  print *, '***PYTHON COMMAND? AMBIGUOUS: ' // command
  return
endif

nl = 0
allocate (lines(200))

select case (command)

!----------------------------------------------------------------------
! Print info on a given plot

case ('plot')

  call tao_find_plots (err, line, 'BOTH', plot, print_flag = .false.)
  if (err) return

  if (allocated(plot)) then
  endif

!----------------------------------------------------------------------
! print list of visible plot names

case ('visible_plots')

!----------------------------------------------------------------------
! print list of plot templates

case ('template_plots')
  do i = 1, size(s%template_plot)
    p => s%template_plot(i)
    if (p%name == '') cycle
    if (p%name == 'scratch') cycle
    if (allocated(p%graph)) then
        nl=nl+1; write (lines(nl), '(100a)') &
                          trim(p%name(1:20)), (':.', trim(p%graph(j)%name), j = 1, size(p%graph))
    else
      nl=nl+1; write (lines(nl), '(3x, a)') p%name 
    endif
  enddo

!----------------------------------------------------------------------

case default

  print *, "***INTERNAL ERROR, SHOULDN'T BE HERE!"

end select

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine
