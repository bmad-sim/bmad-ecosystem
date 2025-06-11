!+
! Subroutine orbit_reference_energy_correction (orbit, dp0c, mat6, make_matrix)
!
! Change orbit%p0c and dependent parameters due to change in reference energy.
!
! Input:
!   orbit         -- coord_struct: Coordinates to correct.
!   dp0c          -- real(rp): Change in reference p0c.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix before correction.
!   make_matrix   -- logical, optional: Propagate the transfer matrix? Default is false.
!
! Output:
!   orbit         -- coord_struct: Coordinates to correct.
!   mat6(6,6)     -- real(rp), optional: Transfer matrix transfer matrix including correction.
!-

subroutine orbit_reference_energy_correction (orbit, dp0c, mat6, make_matrix)

use bmad_struct

implicit none

type (coord_struct) :: orbit

real(rp), optional :: mat6(6,6)
real(rp) dp0c, p0c_new, p_rel

integer particle_at
logical, optional :: make_matrix

character(*), parameter :: r_name = 'orbit_reference_energy_correction'

!

p0c_new = orbit%p0c + dp0c
p_rel = orbit%p0c / p0c_new

if (logic_option(.false., make_matrix)) then
  mat6(2,:) = p_rel * mat6(2,:)
  mat6(4,:) = p_rel * mat6(4,:)
  mat6(6,:) = p_rel * mat6(6,:)
endif

!

orbit%vec(2) = orbit%vec(2) * p_rel
orbit%vec(4) = orbit%vec(4) * p_rel
orbit%vec(6) = (orbit%vec(6) * orbit%p0c - dp0c) / p0c_new
orbit%p0c = p0c_new

end subroutine orbit_reference_energy_correction

