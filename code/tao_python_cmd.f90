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

character(*) input_str
character(200) line_in
character(20) cmd

character(20) :: cmd_names(3) = ['visible_plots', 'graph', 'curve']

integer ix, ix_cmd, nl

!

call string_trim(input_str, line_in, ix)
cmd = line_in(1:ix)
call string_trim(line_in(ix+1:), line_in, ix)

call match_word (cmd, cmd_names, ix, matched_name = cmd)
if (ix == 0) then
  print *, 'PYTHON WHAT? WORD NOT RECOGNIZED: ' // cmd
  return
endif

if (ix < 0) then
  print *, 'PYTHON CMD? AMBIGUOUS: ' // cmd
  return
endif

select case (cmd)

!----------------------------------------------------------------------
! print list of visible plot names

case ('visible_plots')

!----------------------------------------------------------------------
! print list of plot templates

case ('template_plots')


case default

  print *, "INTERNAL ERROR, SHOULDN'T BE HERE!"

end select


end subroutine
