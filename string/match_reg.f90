!+
! function match_reg (str, pat)
!
! Function for matching with regular expressions
!
! Modules needed:
!   use sim_utils
!
! Input:
!   str   -- Character(*): string to test
!   pat   -- Character(*): pattern (regular expression)
!
! Returns .true. if string matches pattern; .false. otherwise.
! *Note: strings are trimmed before comparison
!-

#include "CESR_platform.inc"

function match_reg(str, pat)
  implicit none
  logical :: match_reg
  character(*), intent(in) :: str, pat
  integer stat

  interface
     function match_reg_c(str, pat)
       implicit none
       integer match_reg_c
       character(*) str, pat
     end function match_reg_c
  end interface

  match_reg = .false.

#if defined(CESR_VMS)
  write(*,*) "function match_reg not supported on VMS"
  call err_exit

#else
  stat = match_reg_c(trim(str)//char(0),trim(pat)//char(0))
  if (stat==1) match_reg = .true.
#endif

end function match_reg
