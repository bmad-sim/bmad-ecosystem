!+
! Subroutine tao_fixer (word1, word2)
!
! Fixer commands which is one of:
!   <ele_name>
!   -active {<who>}
!   -store {<who>}
!
! Input:
!   word1     -- character(*): First word of command.
!   word2     -- character(*): Secton word of command.
!-

subroutine tao_fixer (word1, word2)

use tao_interface

implicit none

type (ele_struct), pointer :: ele
character(*) word1, word2
character(8) action
character(*), parameter :: r_name = 'tao_fixer'


!

if (index('-active', trim(word1)) > 1) then
  action = 'active'
elseif (index('-store', trim(word1)) > 1) then
endif

end subroutine
