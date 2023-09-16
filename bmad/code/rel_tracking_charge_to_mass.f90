!+
! Function rel_tracking_charge_to_mass (orbit, ref_species) result (rel_charge)
!
! Routine to determine the relative charge/mass of the particle being 
! tracked relative to the charge of the reference particle.
!
! Input:
!   orbit         -- coord_struct: Particle position structure.
!   ref_species   -- integer: Reference species
!
! Output:
!   rel_charge -- real(rp): Relative charge/mass
!-

function rel_tracking_charge_to_mass (orbit, ref_species) result (rel_charge)

use bmad_struct

implicit none

type (coord_struct) orbit
integer ref_species
real(rp) rel_charge

!

if (orbit%species == photon$ .or. ref_species == photon$) then
  rel_charge = 0
elseif (orbit%species == ref_species) then
  rel_charge = 1
else
  rel_charge = (charge_of(orbit%species) / mass_of(orbit%species)) / (charge_of(ref_species) / mass_of(ref_species))
endif

end function rel_tracking_charge_to_mass

