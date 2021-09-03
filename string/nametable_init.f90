!+
! Subroutine nametable_init(nametable, n_min, n_max)
! 
! Routine to initialize a nametable_struct instance.
!
! Input:
!   n_min    -- integer, optional: Lower bound of array to be indexed. Default is 1.
!   n_max    -- integer, optional: Upper bound of array. Default is n_min - 1.
!
! Output:
!   nametable -- nametable_struct: Variable to be initialized.
!-

subroutine nametable_init (nametable, n_min, n_max)

use re_allocate_mod

implicit none

type (nametable_struct) nametable
integer, optional :: n_min, n_max
integer :: n0, n1

!

n0 = integer_option(1, n_min)
n1 = max(n0-1, integer_option(n0-1, n_max))

nametable%n_min = n0
nametable%n_max = n1

n1 = max(n1, n0+10)
call re_allocate2(nametable%index, n0, n1, .false.)
call re_allocate2(nametable%name, n0, n1, .false.)

end subroutine nametable_init
