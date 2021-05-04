!+
! Subroutine tao_re_allocate_expression_info (info, n, exact, init_val)
!
! Routine to reallocate an array of tao_expression_info_struct structs.
! Note: The data of the array is preserved but data at the end of the
! array will be lost if n is less than the original size of the array
!
! Input:
!   info(:) -- tao_expression_info_struct, allocatable:
!   n       -- Integer: Size wanted.
!   exact   -- Logical, optional: If present and False then the size of 
!                 the output array is permitted to be larger than n. 
!                 Default is True.
!
! Output:
!   info(:) -- Real(rp), Allocatable: Allocated array with size(re) >= n.
!-

subroutine tao_re_allocate_expression_info (info, n, exact)

use tao_struct

implicit none

type (tao_expression_info_struct), allocatable :: info(:), temp_info(:)

integer, intent(in) :: n
integer n_save, n_old

logical, optional :: exact

!

if (allocated(info)) then
  n_old = size(info)
  if (n == n_old) return
  if (.not. logic_option(.true., exact) .and. n < n_old) return
  call move_alloc (info, temp_info)
  allocate (info(n))
  n_save = min(n, n_old)
  info(1:n_save) = temp_info(1:n_save)
  deallocate (temp_info)  

else
  allocate (info(n))
endif

end subroutine tao_re_allocate_expression_info
