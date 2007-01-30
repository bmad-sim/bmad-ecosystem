!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! subroutine calculate_sr_power (lat, orb, direction, power, &
!                                  inside, outside, gen)
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
!   inside  -- wall_struct: inside wall with outline ready
!   outside -- wall_struct: outside wall with outline ready
!   gen    -- general_lat_param_struct: Contains lat name,
!                     vert emittance, and beam current
!
! Output:
!   power(*)  -- ele_power_struct: power radiated from a lat ele
!   inside  -- wall_struct: inside wall with power information
!   outside -- wall_struct: outside wall with power information
!                         
!-

subroutine calculate_sr_power (lat, orb, direction, power, &
                                                       inside, outside, gen)

  use sr_struct
  use sr_interface

  implicit none

  type (lat_struct) lat
  type (coord_struct) orb(0:*)
  type (wall_struct) inside, outside
  type (general_lat_param_struct) gen
  type (ele_power_struct) power(*)

  integer direction, ie

! initialize all accumulated power to 0

  power(1:lat%n_ele_max)%at_wall = 0
  power(1:lat%n_ele_max)%radiated = 0

! loop over all elements

  do ie = 1, lat%n_ele_track
!    print *, lat%ele(ie)%name,',',lat%ele(ie)%type, &
!         ' ele ',ie,' of ',lat%n_ele_track
    call ele_sr_power (lat, ie, orb, direction, power, inside, outside, gen)
  enddo

end subroutine
