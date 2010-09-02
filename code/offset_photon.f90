!+
! Subroutine offset_photon (ele, param, coord, set)
!
! Routine to effectively offset an element by instead offsetting
! the photon position to correspond to the local crystal or mirror coordinates.
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
type (coord_struct), target :: coord

real(rp) c2g, s2g, ct, st, offset(6), tilt, graze
real(rp) off(3), rot(3), project_x(3), project_y(3), project_s(3)
real(rp), pointer :: p(:), vec(:)
complex(rp) efield_x, efield_y, efieldout_x, efieldout_y

logical :: set

!----------------------------------------------------------------
! Set...

p   => ele%value  ! parameter
vec => coord%vec

if (set) then

  if (ele%key == mirror$) then
    graze = p(graze_angle$) + p(graze_angle_err$)
  else  ! Crystal
    graze = p(graze_angle_in$) + p(graze_angle_err$)
  endif

  ! Set: s_offset

  if (p(s_offset_tot$) /= 0) then
    vec(1) = vec(1) + vec(2) * p(s_offset_tot$)
    vec(3) = vec(3) + vec(4) * p(s_offset_tot$)
    vec(5) = vec(5) - p(s_offset_tot$)
  endif

  ! Set: x and y offset 
  ! Set: pitch

  vec(1) = vec(1) - p(x_offset_tot$)
  vec(2) = vec(2) - p(x_pitch_tot$)
  vec(3) = vec(3) - p(y_offset_tot$)
  vec(4) = vec(4) - p(y_pitch_tot$) 
  vec(5) = vec(5) + vec(1) * p(x_pitch_tot$) + vec(3) * p(y_pitch_tot$)

  ! Set: tilt

  call tilt_coords (p(tilt_tot$) + p(tilt_err$), vec, set$)

  ! Set: graze_angle_err

  vec(2) = vec(2) + p(graze_angle_err$)
  vec(5) = vec(5) - vec(1) * p(graze_angle_err$)

  ! Set: intensities
  efield_x = sqrt(coord%intensity_x) * cmplx(cos(coord%phase_x), sin(coord%phase_x) )
  efield_y = sqrt(coord%intensity_y) * cmplx(cos(coord%phase_y), sin(coord%phase_y) )
  tilt = p(tilt_tot$) + p(tilt_err$)
  efieldout_x = cos(tilt) * efield_x + sin(tilt)*efield_y
  efieldout_y = -sin(tilt) * efield_x + cos(tilt)*efield_y
  coord%intensity_x = sqrt(abs(efieldout_x))
  coord%phase_x = atan2(aimag(efieldout_x),real(efieldout_x))
  coord%intensity_y = sqrt(abs(efieldout_y))
  coord%phase_y = atan2(aimag(efieldout_y),real(efieldout_y))
  

!----------------------------------------------------------------
! Unset... 

else

  if (ele%key == mirror$) then
    c2g = cos(2*p(graze_angle$)) 
    s2g = sin(2*p(graze_angle$))
  else
    c2g = cos(2*p(graze_angle_out$)) 
    s2g = sin(2*p(graze_angle_out$))
  endif

  ct = cos(p(tilt$)) 
  st = sin(p(tilt$))

  project_x = (/ c2g * ct**2 + st**2, -ct * st + c2g * ct * st, -ct * s2g /)
  project_y = (/ -ct * st + c2g * ct * st, ct**2 + c2g * st**2, -s2g * st /) 
  project_s = (/ ct * s2g, s2g * st, c2g /)

  ! Unset: graze_angle_err

  vec(2) = vec(2) - p(graze_angle_err$)
  vec(5) = vec(5) + vec(1) * p(graze_angle_err$)

  ! Unset: tilt

  call tilt_coords (p(tilt_tot$), vec, unset$)

  ! Unset: tilt_err
  ! The difference between tilt and tilt_err is that tilt also rotates the output 
  ! laboratory coords but tilt_err does not. 
  ! The difference between tilt_err with Set vs Unset is that the tilt_err is 
  ! expressed in terms of the input lab coords.

  if (p(tilt_err$) /= 0) then

    rot = project_s * p(tilt_err$)

    call tilt_coords (rot(3), vec, unset$)

    vec(2) = vec(2) + rot(2) 
    vec(4) = vec(4) - rot(1)
    vec(5) = vec(5) - vec(1) * rot(2) + vec(3) * rot(1)

  endif

  ! Unset: pitch
  ! Since the pitches are with respect to the lab input coord system, we have
  ! to translate to the local output coords.

  rot = project_x * p(y_pitch_tot$) - project_y * p(x_pitch_tot$)
  vec(2) = vec(2) - rot(2)
  vec(4) = vec(4) + rot(1)
  vec(5) = vec(5) + vec(1) * rot(2) - vec(3) * rot(1)

  ! Unset: offset
  ! Translate offsets to the local output coords.

  off = project_x * p(x_offset_tot$) + project_y * p(y_offset_tot$) + project_s * p(s_offset_tot$)

  vec(1) = vec(1) + off(1)
  vec(3) = vec(3) + off(2)
  if (off(3) /= 0) then
    vec(1) = vec(1) - vec(2) * off(3)
    vec(3) = vec(3) - vec(4) * off(3)
    vec(5) = vec(5) + off(3)
  endif

  ! Unset: intensities

  efield_x = cmplx(sqrt(coord%intensity_x)*cos(coord%phase_x),sqrt(coord%intensity_x)*sin(coord%phase_x) )
  efield_y = cmplx(sqrt(coord%intensity_y)*cos(coord%phase_y),sqrt(coord%intensity_y)*sin(coord%phase_y) )
  tilt = p(tilt_tot$) + p(tilt_err$)
  efieldout_x = cos(tilt) * efield_x - sin(tilt)*efield_y
  efieldout_y = sin(tilt) * efield_x + cos(tilt)*efield_y
  coord%intensity_x = sqrt(abs(efieldout_x))
  coord%phase_x = atan2(aimag(efieldout_x),real(efieldout_x))
  coord%intensity_y = sqrt(abs(efieldout_y))
  coord%phase_y = atan2(aimag(efieldout_y),real(efieldout_y))

endif

end subroutine
                          

