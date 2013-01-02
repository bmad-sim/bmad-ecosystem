!+
! Subroutine offset_photon (ele, coord, set, offset_position_only)
!
! Routine to effectively offset an element by instead offsetting
! the photon position and field to correspond to the local crystal or mirror coordinates.
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
!   set       -- Logical: 
!                   T (= set$)   -> Translate from lab coords to incoming local 
!                                     element coords.
!                   F (= unset$) -> Translate from outgoing local coords to lab coords.
!   offset_position_only
!           -- Logical, optional: If present and True then only offset the position
!                coordinates. This is used for example, aperture calculations where
!                offsetting the field is not needed.
!
! Output:
!     coord -- Coord_struct: Coordinates of particle.
!-

subroutine offset_photon (ele, coord, set, offset_position_only)

use track1_mod, dummy => offset_photon

implicit none

type (ele_struct), target :: ele
type (coord_struct), target :: coord

real(rp) c2g, s2g, ct, st, offset(6), tilt, r(3), ds_center
real(rp) off(3), rot(3), project_x(3), project_y(3), project_s(3), rot_mat(3,3)
real(rp), pointer :: p(:), vec(:)
complex(rp) efield_x, efield_y, efieldout_x, efieldout_y

logical :: set
logical, optional :: offset_position_only
logical is_reflective_element

!

select case (ele%key)
case (crystal$, mirror$, multilayer_mirror$)
  is_reflective_element = .true.
case default
  is_reflective_element = .false.
  ds_center = coord%vec(5) - ele%value(l$)/2
end select

!----------------------------------------------------------------
! Set...

p   => ele%value  ! parameter
vec => coord%vec

if (set) then

  ! Set: z_offset

  if (p(z_offset_tot$) /= 0) then
    vec(1) = vec(1) + vec(2) * p(z_offset_tot$) / vec(6)
    vec(3) = vec(3) + vec(4) * p(z_offset_tot$) / vec(6)
    coord%t = coord%t + p(z_offset_tot$)  / vec(6) / c_light 
    coord%s = coord%s + p(z_offset_tot$)
  endif

  ! Set: X and Y offsets

  vec(1) = vec(1) - p(x_offset_tot$)
  vec(3) = vec(3) - p(y_offset_tot$)

  ! Set: pitch

  if (p(x_offset_tot$) /= 0 .or. p(y_offset_tot$) /= 0) then
    call pitches_to_rotation_matrix (p(x_offset_tot$), p(y_offset_tot$), set, rot_mat)
    r = [vec(1:3:2), ds_center]
    vec(1:5:2) = matmul(rot_mat, r)
    vec(5) = vec(5) + ele%value(l$)/2
    vec(2:6:2) = matmul(rot_mat, vec(2:6:2))
  endif

  ! Set: tilt

  tilt = p(tilt_tot$) + p(tilt_err$)
  if (ele%key == crystal$) tilt = tilt - p(tilt_corr$)

  call tilt_coords (tilt, vec)

  ! Set: graze_angle_err

  if (is_reflective_element) then
    vec(2) = vec(2) + p(graze_angle_err$)
    vec(5) = vec(5) - vec(1) * p(graze_angle_err$)
  endif

  ! Set: intensity and phase rotation due to the tilt + tilt_err

  if (logic_option(.false., offset_position_only)) return

  efield_x = coord%e_field_x * cmplx(cos(coord%phase_x), sin(coord%phase_x) )
  efield_y = coord%e_field_y * cmplx(cos(coord%phase_y), sin(coord%phase_y) )
  efieldout_x = cos(tilt) * efield_x + sin(tilt)*efield_y
  efieldout_y = -sin(tilt) * efield_x + cos(tilt)*efield_y
  coord%e_field_x = abs(efieldout_x)
  coord%phase_x = atan2(aimag(efieldout_x),real(efieldout_x))
  coord%e_field_y = abs(efieldout_y)
  coord%phase_y = atan2(aimag(efieldout_y),real(efieldout_y))

!----------------------------------------------------------------
! Unset... 

else

  ! reflective element...

  if (is_reflective_element) then

    ! Unset: graze_angle_err

    vec(2) = vec(2) - p(graze_angle_err$)
    vec(5) = vec(5) + vec(1) * p(graze_angle_err$)

    ! Unset: tilt

    call tilt_coords (-p(tilt_tot$), vec)

    ! Unset: tilt_err
    ! The difference between tilt and tilt_err is that tilt also rotates the output 
    ! laboratory coords but tilt_err does not. 
    ! The difference between tilt_err with Set vs Unset is that the tilt_err is 
    ! expressed in terms of the input lab coords.

    if (ele%key == mirror$) then
      c2g = cos(2*p(graze_angle$)) 
      s2g = sin(2*p(graze_angle$))
    else
      c2g = cos(p(graze_angle_in$)+p(graze_angle_out$)) 
      s2g = sin(p(graze_angle_in$)+p(graze_angle_out$))
    endif

    ct = cos(p(tilt_tot$)) 
    st = sin(p(tilt_tot$))

    ! [project_x, project_y, project_z] is the matrix product of the matrices T.G.T^-1
    ! T rotates x to y by tilt_tot, G rotates z to x by sum of graze angles
    ! project_x is the projection of the x-basis vector in the original basis onto the new basis

    project_x = [c2g * ct**2 + st**2,      -ct * st + c2g * ct * st, -ct * s2g ]
    project_y = [-ct * st + c2g * ct * st, ct**2 + c2g * st**2,      -s2g * st ] 
    project_s = [ct * s2g,                 s2g * st,                 c2g       ]

    if (p(tilt_err$) /= 0) then

      rot = project_s * p(tilt_err$)

      call tilt_coords (-rot(3), vec)

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

    off = project_x * p(x_offset_tot$) + project_y * p(y_offset_tot$) + project_s * p(z_offset_tot$)

    vec(1) = vec(1) + off(1)
    vec(3) = vec(3) + off(2)
    if (off(3) /= 0) then
      vec(1) = vec(1) - vec(2) * off(3)
      vec(3) = vec(3) - vec(4) * off(3)
      coord%t = coord%t - p(z_offset_tot$)  / vec(6) / c_light 
      coord%s = coord%s - p(z_offset_tot$)
    endif

  ! non-reflective element

  else

    ! Unset: tilt

    tilt = p(tilt_tot$) + p(tilt_err$)
    call tilt_coords (-tilt, vec)

    ! Unset: Pitch

    if (p(x_offset_tot$) /= 0 .or. p(y_offset_tot$) /= 0) then
      call pitches_to_rotation_matrix (-p(x_offset_tot$), -p(y_offset_tot$), set, rot_mat)
      r = [vec(1:3:2), ds_center]
      vec(1:5:2) = matmul(rot_mat, r)
      vec(5) = vec(5) - ele%value(l$)/2
      vec(2:6:2) = matmul(rot_mat, vec(2:6:2))
    endif

    ! Unset: X and Y Offsets

    vec(1) = vec(1) - p(x_offset_tot$)
    vec(3) = vec(3) - p(y_offset_tot$)

    ! Unset S Offset

    if (p(z_offset_tot$) /= 0) then
      vec(1) = vec(1) - vec(2) * p(z_offset_tot$) / vec(6)
      vec(3) = vec(3) - vec(4) * p(z_offset_tot$) / vec(6)
      coord%t = coord%t - p(z_offset_tot$)  / vec(6) / c_light 
      coord%s = coord%s - p(z_offset_tot$)
    endif

  endif

  ! Unset: intensities

  efield_x = cmplx(coord%e_field_x*cos(coord%phase_x), coord%e_field_x*sin(coord%phase_x))
  efield_y = cmplx(coord%e_field_y*cos(coord%phase_y), coord%e_field_y*sin(coord%phase_y))
  tilt = p(tilt_tot$) + rot(3)
  efieldout_x = cos(tilt) * efield_x - sin(tilt)*efield_y
  efieldout_y = sin(tilt) * efield_x + cos(tilt)*efield_y
  coord%e_field_x = abs(efieldout_x)
  coord%phase_x = atan2(aimag(efieldout_x),real(efieldout_x))
  coord%e_field_y = abs(efieldout_y)
  coord%phase_y = atan2(aimag(efieldout_y),real(efieldout_y))

endif

end subroutine
                          

