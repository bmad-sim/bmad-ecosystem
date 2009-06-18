!+
! Subroutine offset_photon (ele, param, coord, set)
!
! Routine to effectively offset an element by instead offsetting
! the photon position to correspond to the local crystal or mirrow coordinates.
!
! Note: The transformation between local and lab coordinates does *not* include 
! a transformation due to a finite graze angle. 
! The reason why the graze angle is ignored is due to the fact that we want 
! the transverse coordinates to be small (so we can do first order optics). 
! However, the graze angle can be large and a transformation with a large angle
! would result in large transverse coords.
! Thus it is up to the calling routine to take care of a finite graze angle.
! [Notice though, that the graze angle is used in the tilt_err transformation.]
!
! Because of the above, there are two sets of local coords: "Incoming" local coords
! and "outgoing" local coords. In both cases the s-axis makes an angle of angle_graze
! with respect to the element surface. 
! The difference is that the angle between the s-axis of the outgoing coords and the 
! s-axis of the incoming coords is 2*angle_graze.
! The x-axis is similarly rotated. The y-axis is the same in both coords.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element
!     %value(x_offset$) -- Horizontal offset of element.
!     %value(x_pitch$)  -- Horizontal roll of element.
!     %value(tilt$)     -- tilt of element.
!     %value(tilt_err$) -- Error tilt of element.
!   coord     -- Coord_struct: Coordinates of the particle.
!     %vec(6)            -- Energy deviation dE/E. 
!                          Used to modify %vec(2) and %vec(4)
!   param     -- lat_param_struct:
!     %particle   -- What kind of particle (for elseparator elements).
!   set       -- Logical: 
!                   T (= set$)   -> Translate from lab coords to incoming local 
!                                     element coords.
!                   F (= unset$) -> Translate from outgoing local coords to lab coords.
!                                               
! Output:
!     coord -- Coord_struct: Coordinates of particle.
!-

subroutine offset_photon (ele, param, coord, set)

use bmad_interface !!!!, except_dummy => offset_photon
use track1_mod, only: track_a_drift

implicit none

type (ele_struct), target :: ele
type (lat_param_struct) :: param
type (coord_struct) :: coord

real(rp) c, s, offset(6), tilt, graze
real(rp), pointer :: v(:), x, xp, y, yp

logical :: set

!----------------------------------------------------------------
! Set...

v => ele%value

if (set) then

  tilt = v(tilt_tot$) + v(tilt_err$)
  graze = v(graze_angle$) + v(graze_angle_err$)

  ! Set: s_offset

  if (v(s_offset_tot$) /= 0) call track_a_drift (coord%vec, v(s_offset_tot$))

  ! Set: Offset and pitch

  if (v(x_offset_tot$) /= 0 .or. v(y_offset_tot$) /= 0 .or. &
            v(x_pitch_tot$) /= 0 .or. v(y_pitch_tot$) /= 0) then
    x => v(x_offset_tot$); xp => v(x_pitch_tot$)
    y => v(y_offset_tot$); yp => v(y_pitch_tot$)
    coord%vec(1) = coord%vec(1) - x
    coord%vec(2) = coord%vec(2) - xp
    coord%vec(3) = coord%vec(3) - y
    coord%vec(4) = coord%vec(4) - yp 
    coord%vec(5) = coord%vec(5) + (-(1 - xp) * x * cos(tilt) + &
                                    (1 - yp) * y * sin(tilt)) / tan(graze) - &
                  (v(g_graze$) * x**2 + v(g_trans$) * y**2) * cos (graze) / 2
  endif

  ! Set: tilt

  call tilt_coords (tilt, coord%vec, set$)

  ! Set: graze_angle_err

  coord%vec(2) = coord%vec(2) - v(graze_angle_err$)

!----------------------------------------------------------------
! Unset... 

else

  tilt = v(tilt_tot$) + v(tilt_err$)
  graze = v(graze_angle$) + v(graze_angle_err$)

  ! unset: graze_angle_err

  coord%vec(2) = coord%vec(2) + v(graze_angle_err$)

  ! Unset: Offset and pitch
  ! Since the offsets and pitches are with respect to the lab input coord system, we have
  ! to translate to the local output coords.

  c = cos(2*v(graze_angle$)) 
  s = sin(2*v(graze_angle$))

  if (v(x_offset_tot$) /= 0 .or. v(y_offset_tot$) /= 0 .or. &
            v(x_pitch_tot$) /= 0 .or. v(y_pitch_tot$) /= 0 .or. v(s_offset_tot$) /= 0) then
    ! First translate to local input coords
    offset = (/ v(x_offset_tot$), v(x_pitch_tot$), &
                v(y_offset_tot$), v(y_pitch_tot$), v(s_offset_tot$), 0.0_rp /)
    call tilt_coords (tilt, offset, set$)

    ! Now apply the offsets
    coord%vec(1) = coord%vec(1) + offset(1) * c + offset(5) * s  ! local x_offset
    coord%vec(2) = coord%vec(2) + offset(2)       ! local x_pitch
    coord%vec(3) = coord%vec(3) + offset(3)       ! local y_offset
    coord%vec(4) = coord%vec(4) + offset(4) * c   ! local y_pitch 

    call track_a_drift (coord%vec, offset(1) * s - offset(5) * c)

  endif

  ! Unset: tilt and tilt_err

  call tilt_coords (v(tilt_tot$) + v(tilt_err$), coord%vec, unset$)

  if (v(tilt_err$) /= 0) then
    coord%vec(2) = coord%vec(2) - s * sin(v(tilt_tot$))
    coord%vec(4) = coord%vec(4) + s * cos(v(tilt_tot$))
  endif



endif

end subroutine
                          

