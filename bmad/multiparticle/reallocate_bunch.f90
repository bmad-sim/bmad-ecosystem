!+
! Subroutine reallocate_bunch (bunch, n_particle, save)
! 
! Subroutine to reallocate particles within a bunch_struct.
!
! Input:
!   n_particle -- Integer: Number of particles. Must be non-negative.
!   save       -- logical, optional: If present and True then save the old bunch info.
!
! Output:
!   bunch      -- bunch_struct: Allocated bunch_struct structure.
!-

subroutine reallocate_bunch (bunch, n_particle, save)

use equal_mod, except_dummy => reallocate_bunch

implicit none

type (bunch_struct) bunch
type (coord_struct), allocatable :: particle(:)

integer i, n_particle, n
logical, optional :: save

! Deallocate if needed

if (allocated(bunch%particle)) then
  if (size(bunch%particle) == n_particle) return
  call move_alloc(bunch%particle, particle)
endif

if (allocated(bunch%ix_z)) deallocate(bunch%ix_z)
allocate (bunch%particle(n_particle), bunch%ix_z(n_particle))
bunch%ix_z = 0

if (logic_option(.false., save) .and. allocated(particle)) then
  n = min(n_particle, size(particle))
  bunch%particle(1:n) = particle(1:n)
endif

if (allocated(particle)) then
  deallocate(particle)
endif

end subroutine reallocate_bunch
