!+
! Subroutine offset_coords_m (ele, param, coord, set, set_canonical, set_multi)
!
! Subroutine to effectively offset an element by instead offsetting
! the particle position to correspond to the local element coordinates.
! This subroutine also optionally transforms from cononical coordinates
! to lab coordinates.
! This is generally used with the subroutine tilt_coords.
!
! Subroutine also puts in half the multipole kick of the element.
!
! Also: Changes coord%x%vel and coord%y%vel from cononical momenta (P_x, P_y) 
! to angle (x', y') and back.
!
! Modules Needed:
!   use bmad_struct
!   use bmad_interface
!
! Input:
!   ELE     -- Ele_struct: Element
!     %VALUE(X_OFFSET$) -- Horizontal offset of element.
!     %VALUE(Y_OFFSET$) -- Vertical offset of element.
!     %VALUE(X_PITCH$)  -- Horizontal roll of element.
!     %VALUE(Y_PITCH$)  -- Vertical roll of element.
!   COORD
!     %Z%VEL            -- Used to modify %X%VEL and %Y%VEL
!   SET     -- Logical: 
!               T -> Translate from lab coords to the local element coords.
!               F -> Translate back to lab coords.
!   SET_CANONICAL -- logical: 
!               T -> Convert between (P_x, P_y) and (x', y') also.
!               F -> No conversion between (P_x, P_y) and (x', y').
! 
!
! Output:
!     COORD -- Coord_struct: Coordinates of particles.
!-

!$Id$
!$Log$
!Revision 1.2  2001/09/27 18:31:55  rwh24
!UNIX compatibility updates
!

#include "CESR_platform.inc"


                                                                 
subroutine offset_coords_m(ele, param, coord, set, set_canonical, set_multi)

  use bmad_struct
  implicit none

  type (ele_struct) ele
  type (coord_struct) coord
  type (param_struct) param

  real E_rel, knl(0:n_pole_maxx), tilt(0:n_pole_maxx)

  integer n

  logical set, set_canonical, set_multi

!---------------------------------------------------------------         
! E_rel               

  E_rel = (1 + coord%z%vel)

!----------------------------------------------------------------
! Set...

  if (set) then

! Set: Offset and pitch

    if (ele%key == multipole$ .or. ele%key == ab_multipole$) then
      coord%x%pos = coord%x%pos - ele%value(x_offset$)
      coord%y%pos = coord%y%pos - ele%value(y_offset$)
    else
      coord%x%pos = coord%x%pos - ele%value(x_offset$) +  &
                         ele%value(x_pitch$) * ele%value(l$) / 2
      coord%x%vel = coord%x%vel - ele%value(x_pitch$) * E_rel
      coord%y%pos = coord%y%pos - ele%value(y_offset$) +  &
                         ele%value(y_pitch$) * ele%value(l$) / 2
      coord%y%vel = coord%y%vel - ele%value(y_pitch$) * E_rel
    endif

! Set: Multipoles

    if (set_multi .and. ele%nonzero_multipoles) then
      call multipole_to_vecs(ele, param%particle, knl, tilt)
      knl = knl / 2
      do n = 0, n_pole_maxx
        call multipole_kick (knl(n), tilt(n), n, coord)
      enddo
    endif

! Set: P_x, P_y -> x', y' 

    if (set_canonical .and. coord%z%vel /= 0) then
      coord%x%vel = coord%x%vel / E_rel
      coord%y%vel = coord%y%vel / E_rel
    endif

!----------------------------------------------------------------
! Unset... 

  else

! Unset: x', y' -> P_x, P_y 

    if (set_canonical .and. coord%z%vel /= 0) then
      coord%x%vel = coord%x%vel * E_rel
      coord%y%vel = coord%y%vel * E_rel
    endif

! Unset: Multipoles

    if (set_multi .and. ele%nonzero_multipoles) then
      do n = 0, n_pole_maxx
        call multipole_kick (knl(n), tilt(n), n, coord)
      enddo
    endif

! Unset: Offset and pitch

    if (ele%key == multipole$ .or. ele%key == ab_multipole$) then
      coord%x%pos = coord%x%pos + ele%value(x_offset$)
      coord%y%pos = coord%y%pos + ele%value(y_offset$)
    else
      coord%x%pos = coord%x%pos + ele%value(x_offset$) + &
                         ele%value(x_pitch$) * ele%value(l$) / 2
      coord%x%vel = coord%x%vel + ele%value(x_pitch$) * E_rel
      coord%y%pos = coord%y%pos + ele%value(y_offset$) +  &
                         ele%value(y_pitch$) * ele%value(l$) / 2
      coord%y%vel = coord%y%vel + ele%value(y_pitch$) * E_rel
    endif

  endif

end subroutine
