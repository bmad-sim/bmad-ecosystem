!+
! Subroutine reallocate_bunch (bunch, n_particle)
! 
! Subroutine to reallocate particles within a bunch_struct.
!
! Input:
!   n_particle -- Integer: Number of particles. Must be non-negative.
!
! Output:
!   bunch -- bunch_struct: Allocated bunch_struct structure.
!-

subroutine reallocate_bunch (bunch, n_particle)

use equal_mod, except_dummy => reallocate_bunch

implicit none

type (bunch_struct) bunch

integer i, n_particle

! Deallocate if needed

if (allocated(bunch%particle)) then
  if (size(bunch%particle) /= n_particle) deallocate (bunch%particle, bunch%ix_z)
endif

if (.not. allocated(bunch%particle)) then
  allocate (bunch%particle(n_particle), bunch%ix_z(n_particle))
  bunch%ix_z = 0
endif

end subroutine reallocate_bunch
