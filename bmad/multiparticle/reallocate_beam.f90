!+
! Subroutine reallocate_beam (beam, n_bunch, n_particle, extend)
! 
! Subroutine to reallocate memory within a beam_struct.
!
! If extend = False (the default), the number of particles in all bunches will be n_particle and
! all particle structures will be initialized (all old data lost).
! If extend = True, the bunches 1 to n_bunch-1 will not be tounched and the last bunch with
! index n_bunch will be initialized to have n_particles.
!
! Note:To extend the number of particles in a given bunch use reallocate_bunch.
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

subroutine reallocate_beam (beam, n_bunch, n_particle, extend)

use equal_mod, except_dummy => reallocate_beam

implicit none

type (beam_struct) beam, temp_beam
integer i, n, n_old, n_bunch
integer, optional :: n_particle
logical, optional :: extend

!

if (n_bunch == 0) then
  if (allocated(beam%bunch)) deallocate(beam%bunch)
  return
endif

! Save present bunches

n_old = 0
if (logic_option(.false., extend) .and. allocated(beam%bunch)) then
  n_old = size(beam%bunch)
  call move_alloc(beam%bunch, temp_beam%bunch)
endif  

! Deallocate and allocate if needed

if (allocated(beam%bunch)) then
  if (size(beam%bunch) /= n_bunch) deallocate (beam%bunch)
endif

if (.not. allocated (beam%bunch)) allocate (beam%bunch(n_bunch))

!

n = min(n_bunch, n_old)
if (n > 0) beam%bunch(1:n) = temp_beam%bunch(1:n)

if (present(n_particle)) then
  do i = n_old+1, n_bunch
    call reallocate_bunch (beam%bunch(i), n_particle)
  enddo
endif

end subroutine reallocate_beam
