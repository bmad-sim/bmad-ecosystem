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

integer i, j
logical good_opt

!

do i = 1, size(s%u)

!return if no variables
  if (.not. associated (s%u(i)%v1_var)) return 

  do j = 1, size(s%u(i)%v1_var)
    good_opt = s%u(i)%v1_var(j)%good_opt
    v => s%u(i)%v1_var(j)%v
    if (.not. associated(v)) cycle
    v(:)%useit_opt = good_opt .and. v(:)%exists .and. &
                        v(:)%good_user .and. v(:)%good_var
  enddo
enddo

end subroutine tao_set_var_useit_opt
