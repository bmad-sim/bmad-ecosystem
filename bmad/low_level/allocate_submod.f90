submodule (bmad_routine_interface) allocate_submod

contains

!----------------------------------------------------------------------------------------------------
!+
! Subroutine allocate_grid_field (g_field, n_gf)
!
! Routine to allocate a new grid_field in the g_field(:) array
!
! Input:
!   g_field(:)      -- grid_field_struct, pointer: Array of grid fields.
!
! Output:
!   g_field(:)      -- grid_field_struct, pointer: Array of grid fields.
!   n_gf            -- integer: size of g_field after allocation
!-



module procedure allocate_grid_field ! (g_field, n_gf)

!  type (grid_field_struct), pointer :: g_field(:)
!  integer n_gf

type (grid_field_struct), pointer :: g_temp(:)

!

if (.not. associated(g_field)) then
  allocate (g_field(1))
  n_gf = 1
  return
endif

n_gf = size(g_field) + 1
g_temp => g_field
allocate (g_field(n_gf))
g_field(1:n_gf-1) = g_temp

end procedure

end submodule
