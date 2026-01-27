!+
! subroutine tao_set_var_useit_opt ()
!
! Calculate which variables will be used by optimizer.
! This is indicated by s%u(:)%var(:)%useit_opt
!
! Input:
!   s     -- type_super_universe_struct
!
! Output:
!   s     -- type_super_universe_struct
!-

subroutine tao_set_var_useit_opt ()

use tao_interface, dummy => tao_set_var_useit_opt

implicit none

type (tao_var_struct), pointer :: v(:)

!

v => s%var
v(:)%useit_opt = v(:)%good_opt .and. v(:)%exists .and. v(:)%good_user .and. v(:)%good_var

end subroutine tao_set_var_useit_opt
