!+
! Subroutine reallocate_beam (beam, n_bunch, n_particle, save)
! 
! Subroutine to reallocate memory within a beam_struct.
!
! If n_bunch = 0 then beam%bunch will be deallocated.
!
! Rule: Bmad routines (except for reallocate_bunch) are allowed to assume that if 
! beam%bunch(:) is allocated, beam%bunch(i)%particle(:) is also allocated.
!
! Input:
!   beam       -- beam_struct: Beam bunches are saved if save = True.
!   n_bunch    -- Integer: Number of bunches.
!   n_particle -- Integer, optional: Number of particles. Must be non-negative. 
!                   If save = True then the number of particles in existing bunches will not be touched.
!                   If not present, beam%bunch(i)%particle(:) will be in an undefined state.
!   save       -- logical, optional: If present and True then save the old bunch info
!
! Output:
!   beam       -- beam_struct: Allocated beam_struct structure.
!-

subroutine reallocate_beam (beam, n_bunch, n_particle, save)

use equal_mod, except_dummy => reallocate_beam

implicit none

type (beam_struct) beam, temp_beam
integer i, n, n_save, n_bunch
integer, optional :: n_particle
logical, optional :: save

! Save present bunches

n_save = 0
if (logic_option(.false., save) .and. allocated(beam%bunch)) then
  n_save = size(beam%bunch)
  call move_alloc(beam%bunch, temp_beam%bunch)
endif  

! Deallocate if needed

if (allocated(beam%bunch)) then
  if (n_bunch == 0 .or. size(beam%bunch) /= n_bunch) deallocate (beam%bunch)
endif

if (n_bunch == 0) return
  
! Allocate

if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))

n = min(n_bunch, n_save)
if (n > 0) beam%bunch(1:n) = temp_beam%bunch(1:n)

if (present(n_particle)) then
  do i = n_save+1, n_bunch
    call reallocate_bunch (beam%bunch(i), n_particle)
  enddo
endif

end subroutine reallocate_beam
