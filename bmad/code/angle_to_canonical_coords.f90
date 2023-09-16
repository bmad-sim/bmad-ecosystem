!+
! Subroutine angle_to_canonical_coords (orbit, coord_type)
!
! Routine to convert from angle coordinates to canonical (x, px, y, py, z, pz) coordinates.
!
! Input:
!   orbit       -- coord_struct: Orbit in angular coordinates.
!   coord_type  -- character(*), optional: Angular coordinates type
!                     '' (default): (x, x' = dx/ds, y, y' = dy/ds, z, pz) 
!                     'ZGOUBI':     (x, x' = dx/ds, y, y' = dy/ds, dt = -z / (beta * c), pz) 
!
! Output:
!   orbit       -- coord_struct: Orbit in canonical coordinates.
!-

subroutine angle_to_canonical_coords (orbit, coord_type)

use bmad_struct
implicit none

type (coord_struct) orbit
real(rp) f
character(*), optional :: coord_type

!

f = (1 + orbit%vec(6)) / sqrt(1 + orbit%vec(2)**2 + orbit%vec(4)**2)
orbit%vec(2) = orbit%vec(2) * f
orbit%vec(4) = orbit%vec(4) * f

if (present(coord_type)) then
  select case (coord_type)
  case ('ZGOUBI')
    orbit%vec(5) = -orbit%vec(5) * (orbit%beta * c_light)
  case ('')
  case default
    call err_exit
  end select
endif

end subroutine angle_to_canonical_coords
