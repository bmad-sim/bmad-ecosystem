!+
! Subroutine ps2gif (ps_file, gif_file, kill_ps_file)
!
! Routine to convert a PS file to a GIF file.
! This routine uses the 'gs' and 'ppmtogif' system commands.
!
! Modules needed:
!   use sim_utils
!
! Input:
!   ps_file      -- Character(*): Name of existing PS file.
!   gif_file     -- Character(*): Name of GIF file to be created.
!   kill_ps_file -- Logical, optional: If True then delete the PS file at the end.
!                     Default is False.
!-

subroutine ps2gif (ps_file, gif_file, kill_ps_file)

use sim_utils, only: system_command, logic_option

implicit none

character(*) ps_file, gif_file
logical, optional :: kill_ps_file

!

call system_command ( &
             'gs -q -sDEVICE=pbm -sOutputFile=temp999.pbm -dNOPAUSE - < ' // ps_file)
call system_command ('ppmtogif temp999.pbm > ' // gif_file)
call system_command ('rm -f temp999.pbm')
if (logic_option(.false., kill_ps_file)) call system_command ('rm -f ' // ps_file)

end subroutine
