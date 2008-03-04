!+
! Subroutine tao_hook_init_global (is_set)
!
! This routine is part of a collection of hook routines 
! used to bypass the use of an initialization file.
!-

subroutine tao_hook_init_global (is_set)

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine



!+
! Subroutine tao_hook_init_connected_uni (is_set)
!-

subroutine tao_hook_init_connected_uni (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine




!+
! Subroutine tao_hook_init_data (is_set) 
!-

subroutine tao_hook_init_data (is_set) 

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine




!+
! Subroutine tao_hook_init_var(is_set) 
!-

subroutine tao_hook_init_var(is_set) 

use tao_mod

implicit none

logical is_set

!

is_set = .false.

end subroutine

!+
!
!-

subroutine tao_hook_init_beam (is_set)

implicit none

logical is_set

!

is_set = .false.

end subroutine

