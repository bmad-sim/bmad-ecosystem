!+
! Subroutine reallocate_beam (beam, n_bunch, n_particle)
! 
! Subroutine to reallocate memory within a beam_struct.
!
! If n_bunch = 0 then all macro beam pointers will be deallocated.
! Rule: If beam%bunch(:) is allocated, beam%bunch(i)%particle(:) will be allocated.
!
! Modules needed:
!   use beam_mod
!
! Input:
!   n_bunch    -- Integer: Number of bunches.
!   n_particle -- Integer: Number of particles. Must be non-negative.
!
! Output:
!   beam -- beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_beam (beam, n_bunch, n_particle)

use basic_bmad_interface, except_dummy => reallocate_beam

implicit none

type (beam_struct) beam

integer i, n_bunch, n_particle

! Deallocate if needed

if (allocated(beam%bunch)) then
  if (n_bunch == 0 .or. size(beam%bunch) /= n_bunch) deallocate (beam%bunch)
endif

if (n_bunch == 0) return
  
! Allocate

if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))

do i = 1, n_bunch
  call reallocate_bunch (beam%bunch(i), n_particle)
enddo

end subroutine reallocate_beam
