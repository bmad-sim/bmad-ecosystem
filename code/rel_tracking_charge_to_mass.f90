!+
! Function rel_tracking_charge_to_mass (orbit, param) result (rel_charge)
!
! Routine to determine the relative charge/mass of the particle being 
! tracked relative to the charge of the reference particle.
!
! Input:
!   orbit -- coord_struct: Particle position structure.
!   param -- lat_param_struct: Structure holding the reference particle id.
!
! Output:
!   rel_charge -- real(rp): Relative charge/mass
!-

function rel_tracking_charge_to_mass (orbit, param) result (rel_charge)

use bmad_struct

implicit none

type (coord_struct) orbit
type (lat_param_struct) param
real(rp) rel_charge

!

if (orbit%species == photon$ .or. param%particle == photon$) then
  rel_charge = 0
elseif (orbit%species == param%particle) then
  rel_charge = 1
else
  rel_charge = (charge_of(orbit%species) / mass_of(orbit%species)) / &
                 (charge_of(param%particle) / mass_of(param%particle))
endif

end function rel_tracking_charge_to_mass

