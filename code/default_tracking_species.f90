!+
! Function default_tracking_species (param) result (species)
!
! Routine to return the default species for tracking through a lattice branch.
!
! Input:
!   param  -- lat_param_struct: Parameters for a lattice branch.
!
! Output:
!   species -- Integer: Default species to be used for tracking.
!-

function default_tracking_species (param) result (species)

use bmad_struct

implicit none

type (lat_param_struct) param
integer species

!

if (param%default_tracking_species == ref_particle$) then
  species = param%particle
elseif (param%default_tracking_species == anti_ref_particle$) then
  species = antiparticle(param%particle)
else
  species = param%default_tracking_species
endif

end function default_tracking_species

