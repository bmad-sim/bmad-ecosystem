!+
! Subroutine synrad3d_photon_init_custom (orbit, ix_branch, lat, photons, n_photon_generated, n_photon_passed, finished)
!
! Routine that can be customized for generating initial photon positions.
! RULE: Only modify arguments that appear in the Output argument list. 
!
! Input:
!   lat                 -- lat_struct: Lattice photon will be tracked through (includes wall information).
!   photons(:)          -- sr3d_photon_track_struct: List of photons passing filter tests
!   n_photon_generated  -- integer: Number of photons generated excluding this one.
!   n_photon_passed     -- ingeger: Number of photons that have passed filter the tests.
!
! Output:
!   orbit         -- coord_struct: Photon initial orbit.
!   ix_branch     -- integer: Index of lattice branch to start photon at.
!   finished      -- logical: Synrad3D will run until this argument is set to True.
!                       When set to True, the current photon will not be tracked.
!-

subroutine synrad3d_photon_init_custom (orbit, ix_branch, lat, photons, n_photon_generated, n_photon_passed, finished)

use synrad3d_struct, dummy => synrad3d_photon_init_custom

implicit none

type (coord_struct) orbit
type (lat_struct) lat
type (sr3d_photon_track_struct) photons(:)

integer ix_branch, n_photon_generated, n_photon_passed
logical finished

! orbit%vec =  [x, vx, y, vy, z, vz]
! The z coordinate is ignored.

orbit%vec = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp]  
orbit%s  = 0   ! Longitudinal position
ix_branch = 0  ! Lattice branch index. 

! Synrad3D will run until the finished argument is set to True.
! photons, n_photon_generated and n_photon_passed can be used to decide when to terminate the program.

finished = .false.

end subroutine synrad3d_photon_init_custom 
