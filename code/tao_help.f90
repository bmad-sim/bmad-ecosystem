!+
! Subroutine tao_help (cmd_line)
!
! Online help for TAO commmands. 
! Interfaces with the documentation.
!
! Input:
!   help_what   -- Character(*): command to query
!
! Output:
!   none
!
!-

subroutine tao_help (help_what)

use tao_mod

implicit none


character(*) :: help_what

character(16) :: r_name = "TAO_HELP"
character(200) lines(5)
integer nl

lines(1) = "TAO: The Tool for Accelerator Optics."
lines(2) = "This is an alpha version: don't expect things to work"
lines(3) = "Please submit bugs to dcs16@cornell.edu or js344@cornell.edu"
lines(4) = "See http://www.lepp.cornell.edu/~dcs/bmad for documentation."
nl = 4

call out_io (s_blank$, r_name, lines(1:nl))

end subroutine tao_help
