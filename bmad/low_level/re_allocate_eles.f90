!+
! Subroutine re_allocate_eles (eles, n, save_old, exact)
!
! Routine to allocate an array of ele_pointer_structs.
!
! Input:
!   eles(:)  -- ele_pointer_struct, allocatable: Array of element pointers with possible old data.
!   n        -- Integer: Array size to set.
!   save_old -- Logical, optional: If present and True then save the old data.
!   exact    -- Logical, optional: If present and True then eles will have size = n
!                 If False (default), reallcation will not be done if eles is already large enough
!
! Output:
!   eles(:) -- ele_pointer_struct, allocatable: Array of element pointers.
!-

subroutine re_allocate_eles (eles, n, save_old, exact)

use bmad_struct

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (ele_pointer_struct), allocatable :: l_temp(:)
integer n, n_old, nn
logical, optional :: save_old, exact

!

if (.not. allocated(eles)) then
  allocate (eles(n))
  return
endif

if (size(eles) == n .or. (.not. logic_option(.false., exact) .and. size(eles) >= n)) then
  if (.not. logic_option (.false., save_old)) eles = ele_pointer_struct()
  return
endif

!

if (logic_option (.false., save_old)) then
  call move_alloc (eles, l_temp)
  allocate (eles(n))
  nn = min(size(l_temp), n)
  eles(1:nn) = l_temp(1:nn)
  deallocate (l_temp)
else
  deallocate (eles)
  allocate (eles(n))
endif

end subroutine re_allocate_eles

