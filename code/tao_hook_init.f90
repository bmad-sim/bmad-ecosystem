!+
! Subroutine tao_hook_init (s)
!
! Custom subroutine to initialize the tao structures.
! This is a dummy subroutine that can be over written.
!
! Output:
!   s -- Tao_super_universe_struct:
!-

subroutine tao_hook_init (s)

use tao_mod

implicit none

type (tao_super_universe_struct) s

!

return

end subroutine
