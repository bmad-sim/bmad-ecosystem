!+
! subroutine tao_set_var_useit_opt (s)
!
! Calculate which variables will be used by optimizer.
! This is indicated by s%u(:)%var(:)%useit_opt
!
! Input/Output:
!   s     -- type_super_universe_struct
!-

subroutine tao_set_var_useit_opt (s)

use tao_mod

implicit none

type (tao_super_universe_struct), target :: s
type (tao_var_struct), pointer :: v(:)

!

v => s%var
v(:)%useit_opt = v(:)%good_opt .and. v(:)%exists .and. &
                            v(:)%good_user .and. v(:)%good_var

end subroutine tao_set_var_useit_opt
