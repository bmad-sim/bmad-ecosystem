!+
! Subroutine rotate_for_curved_surface (ele, orbit, set, rot_mat)
!
! Routine to rotate just the velocity coords between element body coords and effective 
! body coords ("curved body coords") with respect to the surface at the point of photon impact.
!
! To rotate the photon coords back to the the element body coords, the inverse of the rotation 
! matrix usedto transform from element body coords is needed. Thus rot_mat must be saved between 
! calls to this routine with set = True and set = False.
!
! Input:
!   ele          -- ele_struct: reflecting element
!   orbit        -- coord_struct: Photon position.
!   set          -- Logical: True -> Transform body coords to local curved body coords. 
!                            False -> Transform local curved body to body coords.
!   rot_mat(3,3) -- real(rp): When set = False, rotation matrix calculated from previous call with set = True.
!
! Output:
!   orbit        -- coord_struct: Photon position.
!   rot_mat(3,3) -- real(rp): When set = True, calculated rotation matrix.
!-

subroutine rotate_for_curved_surface (ele, orbit, set, rot_mat)

use photon_utils_mod, dummy => rotate_for_curved_surface

implicit none

type (ele_struct), target :: ele
type (coord_struct) orbit
type (photon_element_struct), pointer :: ph
type (surface_grid_pt_struct), pointer :: pt

real(rp) rot_mat(3,3)
real(rp) rot(3,3), angle
real(rp) dz_dxy(2), x, y, z, zt, g(3), gs
integer ix, iy

logical set, err_flag
character(*), parameter :: r_name = 'rotate_for_curved_surface'

! Transform from local curved to body coords.

if (.not. set) then
  rot = transpose(rot_mat)
  orbit%vec(2:6:2) = matmul(rot, orbit%vec(2:6:2))
  return
endif

! Compute the slope of the surface at that the point of impact.
! curve_rot transforms from standard body element coords to body element coords at point of impact.

ph => ele%photon
x = orbit%vec(1)
y = orbit%vec(3)

if (ph%grid%type == segmented$ .and. ph%grid%active) then
  pt => pointer_to_surface_grid_pt(ele, .true., x, y)
  if (.not. associated(pt)) then
    orbit%state = lost$
    call out_io (s_info$, r_name, 'Photon is outside of grid bounds for: ' // ele%name)
    return
  endif

  dz_dxy = [pt%dz_dx, pt%dz_dy]

else
  if (ph%grid%type == displacement$) then
    call surface_grid_displacement (ele, x, y, err_flag, z, dz_dxy)
    if (err_flag) then
      orbit%state = lost$
      call out_io (s_info$, r_name, 'Photon is outside of grid bounds for: ' // ele%name)
      return
    endif
  else
    dz_dxy = 0
  endif

  do ix = 0, ubound(ph%curvature%xy, 1)
  do iy = 0, ubound(ph%curvature%xy, 2) - ix
    if (ph%curvature%xy(ix, iy) == 0) cycle
    if (ix > 0) dz_dxy(1) = dz_dxy(1) - ix * ph%curvature%xy(ix, iy) * x**(ix-1) * y**iy
    if (iy > 0) dz_dxy(2) = dz_dxy(2) - iy * ph%curvature%xy(ix, iy) * x**ix * y**(iy-1)
  enddo
  enddo

  g = ph%curvature%elliptical
  if (g(3) /= 0) then
    zt = sqrt(1 - (x * g(1))**2 - (y * g(2))**2)
    dz_dxy = dz_dxy - [x*g(1)**2, y*g(2)**2] / (g(3) * zt)
  endif

  gs = ph%curvature%spherical
  if (gs /= 0) then
    zt = sqrt(1 - (x * gs)**2 - (y * gs)**2)
    dz_dxy = dz_dxy - [x, y] * (gs**2 / (gs * zt))
  endif
endif

if (dz_dxy(1) == 0 .and. dz_dxy(2) == 0) then
  call mat_make_unit(rot_mat)
  return
endif

! Compute rotation matrix and goto body element coords at point of photon impact

angle = -atan2(sqrt(dz_dxy(1)**2 + dz_dxy(2)**2), 1.0_rp)
call axis_angle_to_w_mat ([dz_dxy(2), -dz_dxy(1), 0.0_rp], angle, rot_mat)

orbit%vec(2:6:2) = matmul(rot_mat, orbit%vec(2:6:2))

end subroutine rotate_for_curved_surface

