!+
! function match_reg (str, pat)
!
! Function for matching with regular expressions
!
! Input:
!   str   -- Character(*): string to test
!   pat   -- Character(*): pattern (regular expression)
!
! Returns .true. if string matches pattern; .false. otherwise.
! *Note: strings are trimmed before comparison
!-

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

stat = match_reg_c(trim(str)//char(0),trim(pat)//char(0))
if (stat==1) match_reg = .true.

end function match_reg
