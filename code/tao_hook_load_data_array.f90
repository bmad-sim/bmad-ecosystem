!+
! Subroutine tao_hook_load_data_array (s, found)
!-

subroutine tao_hook_load_data_array (s, found)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
logical found

!

found = .false.

end subroutine