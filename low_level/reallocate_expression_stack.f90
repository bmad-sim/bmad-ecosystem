!+ 
! Subroutine reallocate_expression_stack (stack, n, exact)
!
! Routine to resize an expression stack array.
!
! Input:
!   stack(:)  -- expression_atom_stuct, allocatable: Existing stack array.
!   n         -- integer: Array size needed.
!   exact     -- logical, optional: If present and False then the size of the
!                  output array is permitted to be larger than n.
!                  Default is True.
! Output:
!   stack(:)  -- expression_atom_struct, allocatable: Resized stack.
!                  The stack info will be preserved.
!-

subroutine reallocate_expression_stack (stack, n, exact)

use bmad_struct

implicit none

type (expression_atom_struct), allocatable :: stack(:)
type (expression_atom_struct), allocatable :: stack_temp(:)
integer n, nn
logical, optional :: exact

!

if (.not. allocated(stack)) allocate(stack(n))

if (size(stack) == n .or. size(stack) > n .and. .not. logic_option(.true., exact)) return

call move_alloc (stack, stack_temp)
allocate (stack(n))
nn = min(size(stack_temp), size(stack))
stack(1:nn) = stack_temp(1:nn)

end subroutine reallocate_expression_stack
