!+
! subroutine calculate_synrad_power (lat, orb, direction, power, walls, sr_param, use_ele_ix)
!
! subroutine to calculate the synch radiation power
!   hitting wall segments from all elements in the lat
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   lat        -- lat_struct with twiss propagated and mat6s made
!   orb(0:*)   -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction  -- integer: +1 for in direction of s
!                         -1 for against s
!   walls      -- walls_struct: both walls with outlines ready
!   sr_param   -- synrad_param_struct: Contains lat name,
!                     vert emittance, and beam current
!   use_ele_ix -- calc power from only this element number, or all if 0
!
! Output:
!   power(:) -- ele_power_struct: power radiated from a lat ele
!   walls    -- wall_struct: both wall with power information
!                         
!-

subroutine calculate_synrad_power (lat, orb, direction, power, walls, sr_param, use_ele_ix)

use synrad_struct
use synrad_interface, except => calculate_synrad_power

implicit none

type (lat_struct) lat
type (coord_struct) orb(0:)
type (walls_struct), target :: walls
type (wall_struct), pointer :: negative_x_wall, positive_x_wall
type (synrad_param_struct) sr_param
type (ele_power_struct) power(:)

integer direction, use_ele_ix, ie

! set pointers

positive_x_wall => walls%positive_x_wall
negative_x_wall => walls%negative_x_wall

! initialize all accumulated power to 0

power(1:lat%n_ele_max)%at_wall = 0
power(1:lat%n_ele_max)%radiated = 0

! loop over all elements

do ie = 1, lat%n_ele_track
  if (use_ele_ix > 0 .and. ie /= use_ele_ix) cycle
  if (sr_param%debug) print '(i6, 2x, a)', ie, lat%ele(ie)%name
  call ele_synrad_power (lat, ie, orb, direction, power, walls, sr_param)
enddo

end subroutine
