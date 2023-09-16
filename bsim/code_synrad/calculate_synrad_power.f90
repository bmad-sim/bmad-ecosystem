!+
! subroutine calculate_synrad_power (branch, orb, direction, power, walls, sr_param, ix_ele1, ix_ele2)
!
! subroutine to calculate the synch radiation power
!   hitting wall segments from all elements in the lattice
!
! Modules needed:
!   use synrad_mod
!
! Input:
!   branch        -- branch_struct with twiss propagated and mat6s made
!   orb(0:*)   -- coord_struct: orbit of particles to use as 
!                             source of ray
!   direction  -- integer: +1 for in direction of s
!                         -1 for against s
!   walls      -- walls_struct: both walls with outlines ready
!   sr_param   -- synrad_param_struct: Contains lattice name,
!                     vert emittance, and beam current
!   ix_ele1    -- Integer: Element index for start of region over which to calculate the SR power.
!                   The region includes this element.
!   ix_ele2    -- Integer: Element index for End of region over which to calculate the SR poser.
!
! Output:
!   power(:) -- ele_power_struct: power radiated from a lattice ele
!   walls    -- wall_struct: both wall with power information
!                         
!-

subroutine calculate_synrad_power (branch, orb, direction, power, walls, sr_param, ix_ele1, ix_ele2)

use synrad_struct
use synrad_interface, except => calculate_synrad_power

implicit none

type (branch_struct), target :: branch
type (coord_struct) orb(0:)
type (walls_struct), target :: walls
type (wall_struct), pointer :: negative_x_wall, positive_x_wall
type (synrad_param_struct) sr_param
type (ele_power_struct) power(:)

real(rp) time

integer direction, ix_ele1, ix_ele2, ie

! set pointers

positive_x_wall => walls%positive_x_wall
negative_x_wall => walls%negative_x_wall

! initialize all accumulated power to 0

power(1:branch%n_ele_max)%at_wall = 0
power(1:branch%n_ele_max)%radiated = 0

! loop over all elements

do ie = ix_ele1, ix_ele2
  if (sr_param%debug) then
    print '(i6, 2x, a)', ie, branch%ele(ie)%name
    call run_timer ('START')
  endif

  call ele_synrad_power (branch, ie, orb, direction, power, walls, sr_param)

  if (sr_param%debug) then
    call run_timer ('READ', time)
    print '(a, f10.2)', 'dTime: (min)', time/60
  endif
enddo

end subroutine
