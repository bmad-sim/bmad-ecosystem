!+
! Subroutine unlink_fieldmap (cartesian_map, cylindrical_map, gen_grad_map, grid_field)
!
! Subroutine to unlink the field components of an element.
!
! Input:
!   cartesian_map(:)   -- cartesian_map_struct, pointer, optional: cartesian_map array.
!   cylindrical_map(:) -- cylindrical_map_struct, pointer, optional: cylindrical_map array.
!   gen_grad_map(:)    -- gen_grad_map_struct, pointer: Gen grad array.
!   grid_field(:)      -- grid_field_struct, pointer, optional: grid_field array.
!
! Output:
!   cartesian_map(:)   -- cartesian_map_struct, pointer, optional: cartesian_map array.
!   cylindrical_map(:) -- cylindrical_map_struct, pointer, optional: cylindrical_map array.
!   gen_grad_map(:)    -- gen_grad_map_struct, pointer: Gen grad array.
!   grid_field(:)      -- grid_field_struct, pointer, optional: grid_field array.
!-

subroutine unlink_fieldmap (cartesian_map, cylindrical_map, gen_grad_map, grid_field)

use bmad_struct

type (cartesian_map_struct), pointer, optional :: cartesian_map(:)
type (cylindrical_map_struct), pointer, optional :: cylindrical_map(:)
type (gen_grad_map_struct), pointer, optional :: gen_grad_map(:)
type (grid_field_struct), pointer, optional :: grid_field(:)

integer i

! In theory, %prt should always be associated but a bad digested file can leave %ptr unassociated.

if (present(cartesian_map)) then
  do i = 1, size(cartesian_map)
    if (.not. associated(cartesian_map(i)%ptr)) cycle
    cartesian_map(i)%ptr%n_link = cartesian_map(i)%ptr%n_link - 1
    if (cartesian_map(i)%ptr%n_link == 0) deallocate (cartesian_map(i)%ptr)
  enddo
  deallocate (cartesian_map)
endif

!

if (present(cylindrical_map)) then
  do i = 1, size(cylindrical_map)
    if (.not. associated(cylindrical_map(i)%ptr)) cycle
    cylindrical_map(i)%ptr%n_link = cylindrical_map(i)%ptr%n_link - 1
    if (cylindrical_map(i)%ptr%n_link == 0) deallocate (cylindrical_map(i)%ptr)
  enddo
  deallocate (cylindrical_map)
endif

! Note: gen_grad_map does not use links.

if (present(gen_grad_map)) then
  do i = 1, size(gen_grad_map)
    deallocate (gen_grad_map(i)%gg)
  enddo
  deallocate (gen_grad_map)
endif

!

if (present(grid_field)) then
  do i = 1, size(grid_field)
    if (.not. associated(grid_field(i)%ptr)) cycle
    grid_field(i)%ptr%n_link = grid_field(i)%ptr%n_link - 1
    if (grid_field(i)%ptr%n_link == 0) deallocate (grid_field(i)%ptr)
  enddo
  deallocate (grid_field)
endif

end subroutine unlink_fieldmap

