!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine calculate_sr_power (lat, orb, direction, power, &
!                                  walls, gen)
!
! subroutine to calculate the synch radiation power
!   hitting wall segments from all elements in the lat
!
! Modules needed:
!   use sr_mod
!
! Input:
!   lat   -- lat_struct with twiss propagated and mat6s made
!   orb(0:*) -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction -- integer: +1 for in direction of s
!                         -1 for against s
!   walls  -- walls_struct: both walls with outlines ready
!   gen    -- synrad_param_struct: Contains lat name,
!                     vert emittance, and beam current
!
! Output:
!   power(:) -- ele_power_struct: power radiated from a lat ele
!   walls   -- wall_struct: both wall with power information
!                         
!-

subroutine calculate_sr_power (lat, orb, direction, power, &
     walls, gen)

  use sr_struct
  use sr_interface

  implicit none

  type (lat_struct) lat
  type (coord_struct) orb(0:)
  type (walls_struct), target :: walls
  type (wall_struct), pointer :: negative_x_wall, positive_x_wall
  type (synrad_param_struct) gen
  type (ele_power_struct) power(:)

  integer direction, ie

! set pointers
  positive_x_wall => walls%positive_x_wall
  negative_x_wall => walls%negative_x_wall

! initialize all accumulated power to 0

  power(1:lat%n_ele_max)%at_wall = 0
  power(1:lat%n_ele_max)%radiated = 0

! loop over all elements

  do ie = 1, lat%n_ele_track
!    print *, lat%ele(ie)%name,',',lat%ele(ie)%type, &
!         ' ele ',ie,' of ',lat%n_ele_track
    call ele_sr_power (lat, ie, orb, direction, power, walls, gen)
  enddo

end subroutine
