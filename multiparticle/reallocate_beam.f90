!+
! Subroutine reallocate_beam (beam, n_bunch, n_particle)
! 
! Subroutine to reallocate memory within a beam_struct.
!
! If n_bunch = 0 then beam%bunch will be deallocated.
!
! Rule: Bmad routines (except for reallocate_bunch) are allowed to assume that if 
! beam%bunch(:) is allocated, beam%bunch(i)%particle(:) is also allocated.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   n_bunch    -- Integer: Number of bunches.
!   n_particle -- Integer, optional: Number of particles. Must be non-negative.
!                   If not present, beam%bunch(i)%particle(:) will be in an undefined state.
!
! Output:
!   beam -- beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_beam (beam, n_bunch, n_particle)

use equal_mod, except_dummy => reallocate_beam

implicit none

type (beam_struct) beam

integer i, n_bunch
integer, optional :: n_particle

! Deallocate if needed

if (allocated(beam%bunch)) then
  if (n_bunch == 0 .or. size(beam%bunch) /= n_bunch) deallocate (beam%bunch)
endif

if (n_bunch == 0) return
  
! Allocate

if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))

if (present(n_particle)) then
  do i = 1, n_bunch
    call reallocate_bunch (beam%bunch(i), n_particle)
  enddo
endif

end subroutine reallocate_beam
