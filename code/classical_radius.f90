!+
! Function classical_radius (species) result (radius)
!
! Routine to return the classical radius of a particle.
! In particular, if the particle is an electron, the classical electron radius is returned.
!
! Input:
!   species     -- Integer: Species of particle.
!
! Output:
!   radius      -- real(rp): Classical radius.
!-

function classical_radius (species) result (radius)

use sim_utils

integer species
real(rp) radius

!

radius = classical_radius_factor * charge_of(species)**2 / mass_of(species)

end function
