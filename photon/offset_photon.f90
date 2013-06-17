!+
! Subroutine offset_photon (ele, orbit, set, offset_position_only)
!
! Routine to effectively offset an element by instead offsetting
! the photon position and field to correspond to the local crystal or mirror coordinates.
!
! Modules Needed:
!   use bmad
!
! Input:
!   ele       -- Ele_struct: Element
!   orbit     -- Coord_struct: Coordinates of the particle.
!   set       -- Logical: 
!                   T (= set$)   -> Translate from lab coords to local 
!                                     element coords.
!                   F (= unset$) -> Translate from outgoing local coords to lab coords.
!   offset_position_only
!           -- Logical, optional: If present and True then only offset the position
!                coordinates. This is used for example, aperture calculations where
!                offsetting the field is not needed.
!
! Output:
!     orbit -- Coord_struct: Coordinates of particle.
!-

subroutine offset_photon (ele, orbit, set, offset_position_only)

use track1_mod, dummy => offset_photon

implicit none

type (ele_struct), target :: ele
type (coord_struct), target :: orbit

real(rp) graze2, offset(6), tilt, r(3), ds_center, graze, sin_g, cos_g
real(rp) off(3), rot(3), project(3,3), rot_mat(3,3)
real(rp), pointer :: p(:), vec(:)

complex(rp) field(2)

logical :: set
logical, optional :: offset_position_only
logical is_reflective_element

!

select case (ele%key)
case (crystal$, mirror$, multilayer_mirror$)
  is_reflective_element = .true.
  ds_center = 0
case default
  is_reflective_element = .false.
  ds_center = orbit%vec(5) - ele%value(l$)/2
end select

p   => ele%value  ! parameter
vec => orbit%vec

!----------------------------------------------------------------
! Set...

if (set) then

  ! Set: z_offset

  if (p(z_offset_tot$) /= 0) then
    vec(1) = vec(1) + vec(2) * p(z_offset_tot$) / vec(6)
    vec(3) = vec(3) + vec(4) * p(z_offset_tot$) / vec(6)
    orbit%t = orbit%t + p(z_offset_tot$)  / vec(6) / c_light 
    orbit%s = orbit%s + p(z_offset_tot$)
  endif

  ! Set: X and Y offsets

  vec(1) = vec(1) - p(x_offset_tot$)
  vec(3) = vec(3) - p(y_offset_tot$)

  ! Set: pitch

  if (p(x_pitch_tot$) /= 0 .or. p(y_pitch_tot$) /= 0) then
    call pitches_to_rotation_matrix (p(x_pitch_tot$), p(y_pitch_tot$), set, rot_mat)
    r = [vec(1:3:2), ds_center]
    vec(1:5:2) = matmul(rot_mat, r)
    vec(5) = vec(5) + ele%value(l$)/2
    vec(2:6:2) = matmul(rot_mat, vec(2:6:2))
  endif

  ! Set: tilt

  select case (ele%key)
  case (crystal$)
    tilt = p(ref_tilt_tot$) + p(tilt_tot$) + p(tilt_corr$)
  case (mirror$, multilayer_mirror$)
    tilt = p(ref_tilt_tot$) + p(tilt_tot$)
  case default
    tilt = p(tilt_tot$)
  end select

  call tilt_coords (tilt, vec)

  ! Set: intensity and phase rotation due to the tilt

  if (logic_option(.false., offset_position_only)) return

  field = [orbit%field(1) * cmplx(cos(orbit%phase(1)), sin(orbit%phase(1))), &
           orbit%field(2) * cmplx(cos(orbit%phase(2)), sin(orbit%phase(2)))]

  field = [cos(tilt) * field(1) + sin(tilt)*field(2), &
          -sin(tilt) * field(1) + cos(tilt)*field(2)]

  orbit%field(1) = abs(field(1))
  orbit%phase(1) = atan2(aimag(field(1)), real(field(1)))

  orbit%field(2) = abs(field(2))
  orbit%phase(2) = atan2(aimag(field(2)), real(field(2)))

  ! Set: Rotate by the graze angle to body coords. Note vec(5) = 0 initially.

  select case (ele%key)
  case (crystal$)
    graze = p(graze_angle_in$) 
  case (mirror$, multilayer_mirror$)
    graze = p(graze_angle$) 
  case default
    graze = 0
  end select

  if (graze /= 0) then
    sin_g = sin(graze)
    cos_g = cos(graze)

    orbit%vec(2:6:2) =  [cos_g * orbit%vec(2) + sin_g * orbit%vec(6), orbit%vec(4), &
                        -sin_g * orbit%vec(2) + cos_g * orbit%vec(6)]
    orbit%vec(1:5:2) = [cos_g * orbit%vec(1), orbit%vec(3), -sin_g * orbit%vec(1)]
  endif

!----------------------------------------------------------------
! Unset... 

else

  ! reflective element...

  if (is_reflective_element) then

    select case (ele%key)
    case (crystal$)
      graze = p(graze_angle_out$) 
    case (mirror$, multilayer_mirror$)
      graze = p(graze_angle$) 
    end select

    sin_g = sin(graze)
    cos_g = cos(graze)

    ! Translate momentum to laboratory exit coords
    ! and compute position, backpropagating the ray.

    orbit%vec(2:6:2) = [orbit%vec(2) * cos_g + orbit%vec(6) * sin_g, orbit%vec(4), &
                       -orbit%vec(2) * sin_g + orbit%vec(6) * cos_g]

    orbit%vec(1:5:2) = [cos_g * orbit%vec(1) + sin_g * orbit%vec(5), orbit%vec(3), &
                       -sin_g * orbit%vec(1) + cos_g * orbit%vec(5)]

    orbit%vec(1:5:2) = orbit%vec(1:5:2) - orbit%vec(2:6:2) * (orbit%vec(5) / orbit%vec(6))

    ! Unset: tilt

    call tilt_coords (-p(ref_tilt_tot$), vec)

    ! Unset: tilt_tot
    ! The difference between ref_tilt_tot and tilt_tot is that ref_tilt_tot also rotates the output 
    ! laboratory coords but tilt_tot does not. 
    ! The difference between tilt_tot with Set vs Unset is that the tilt_tot is 
    ! expressed in terms of the input lab coords.

    select case (ele%key)
    case (mirror$, multilayer_mirror$)
      graze2 = 2*p(graze_angle$)
      tilt = p(tilt_tot$)
    case (crystal$)
      graze2 = p(graze_angle_in$)+p(graze_angle_out$)
      tilt = p(tilt_tot$) + p(tilt_corr$)
    end select

    ! project is the entrance coords in terms of the exit coords.
    ! EG: project(1,:) is the entrance x-axis in the exit coords.

    rot = [sin(p(ref_tilt_tot$)), -cos(p(ref_tilt_tot$)), 0.0_rp]
    call axis_angle_to_w_mat (rot, graze2, project)

    if (tilt /= 0) then
      call axis_angle_to_w_mat (project(3,:), tilt, rot_mat)
      vec(1:5:2) = matmul(rot_mat, vec(1:5:2))
      vec(2:6:2) = matmul(rot_mat, vec(2:6:2))
    endif

    ! Unset: pitch
    ! Since the pitches are with respect to the lab input coord system, we have
    ! to translate to the local output coords.

    if (p(x_pitch_tot$) /= 0) then
      call axis_angle_to_w_mat (project(3,:), p(x_pitch_tot$), rot_mat)
      vec(1:5:2) = matmul(rot_mat, vec(1:5:2))
      vec(2:6:2) = matmul(rot_mat, vec(2:6:2))
    endif

    if (p(y_pitch_tot$) /= 0) then
      call axis_angle_to_w_mat (-project(2,:), p(y_pitch_tot$), rot_mat)
      vec(1:5:2) = matmul(rot_mat, vec(1:5:2))
      vec(2:6:2) = matmul(rot_mat, vec(2:6:2))
    endif

    ! Unset: offset
    ! Translate offsets to the local output coords.

    off = project(1,:) * p(x_offset_tot$) + project(2,:) * p(y_offset_tot$) + &
          project(3,:) * p(z_offset_tot$)

    vec(1) = vec(1) + off(1)
    vec(3) = vec(3) + off(2)
    if (off(3) /= 0) then
      vec(1) = vec(1) - vec(2) * off(3)
      vec(3) = vec(3) - vec(4) * off(3)
      orbit%t = orbit%t - p(z_offset_tot$)  / vec(6) / c_light 
      orbit%s = orbit%s - p(z_offset_tot$)
    endif

  ! non-reflective element

  else

    ! Unset: tilt

    tilt = p(tilt_tot$)
    call tilt_coords (-tilt, vec)
    rot = 0

    ! Unset: Pitch

    if (p(x_pitch_tot$) /= 0 .or. p(y_pitch_tot$) /= 0) then
      call pitches_to_rotation_matrix (-p(x_pitch_tot$), -p(y_pitch_tot$), set, rot_mat)
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
      orbit%t = orbit%t - p(z_offset_tot$)  / vec(6) / c_light 
      orbit%s = orbit%s - p(z_offset_tot$)
    endif

  endif

  ! Unset: intensities

  field = [cmplx(orbit%field(1)*cos(orbit%phase(1)), orbit%field(1)*sin(orbit%phase(1))), &
           cmplx(orbit%field(2)*cos(orbit%phase(2)), orbit%field(2)*sin(orbit%phase(2)))]

  tilt = tilt + rot(3)
  field = [cos(tilt) * field(1) - sin(tilt)*field(2), &
           sin(tilt) * field(1) + cos(tilt)*field(2)]

  orbit%field(1) = abs(field(1))
  orbit%phase(1) = atan2(aimag(field(1)),real(field(1)))

  orbit%field(2) = abs(field(2))
  orbit%phase(2) = atan2(aimag(field(2)),real(field(2)))

endif

end subroutine
                          

