
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Function coords_floor_to_local_curvilinear  (global_position, ele, status, w_mat, use_patch_entrance) result(local_position)
!
! Given a position in global coordinates, return local curvilinear coordinates defined by ele.
!
! If the calculated longitudinal s places global_position outside of ele, the status argument
! will be set appropriately and local_position is not will be defined.
!
! Angular orientation is ignored.
!
! Input:
!   global_position     -- floor_position_struct: %r = [X, Y, Z] position in global coordinates
!   ele                 -- ele_struct: element to find local coordinates of.
!   use_patch_entrance  -- logical, optional: This argument is ignored for non-patch elements. If in a patch: 
!                             True => use entrance coordinates. False (default) => use exit coordinates. 
!                            
! Output:
!   local_position  -- floor_position_struct: %r = [x, y, z] position in local curvilinear coordinates
!                        with z relative to entrance edge of the element.
!   status          -- logical: longitudinal position:
!                               inside$: Inside the element.
!                               upstream_end$: Upstream of element.
!                               downstream_end$: Downstream of element.
!   w_mat(3,3)      -- real(rp) (optional): W matrix at s, to transform vectors. 
!                                  v_global = w_mat.v_local
!                                  v_local = transpose(w_mat).v_global
!-  

function coords_floor_to_local_curvilinear (global_position, ele, status, w_mat, use_patch_entrance) result(local_position)

use nr, only: zbrent
use bmad_interface, dummy => coords_floor_to_local_curvilinear

implicit none

type (floor_position_struct) :: global_position, local_position
type (ele_struct)   :: ele
type (floor_position_struct) :: floor0, floor1
real(rp), optional :: w_mat(3,3)
real(rp) x, y, z, rho, dtheta, tilt, dz0, dz1
integer :: status
logical, optional :: use_patch_entrance

! In all cases remember that ele%floor is the downstream floor coords independent of ele%orientation
! sbend case.

if (ele%key == sbend$ .and. ele%value(g$) /= 0) then
  if (ele%orientation == 1) then
    floor0 = ele%branch%ele(ele%ix_ele-1)%floor        ! Get floor0 from previous element
  else
    floor0 = ele%floor
  endif
  local_position = coords_floor_to_relative (floor0, global_position)
  tilt = ele%value(ref_tilt_tot$)
  if (tilt == 0) then
    x = local_position%r(1)
    y = local_position%r(2)
  else
    x =  local_position%r(1) * cos(tilt) + local_position%r(2) * sin(tilt)
    y = -local_position%r(1) * sin(tilt) + local_position%r(2) * cos(tilt)
  endif
  z = local_position%r(3)
  rho = ele%value(rho$)
  if (rho > 0) then
    dtheta = atan2 (z, x + rho)
  else
    dtheta = atan2 (-z, -(x + rho))
  endif

  if (dtheta < ele%value(angle$)/2 - pi) dtheta = dtheta + twopi
  if (dtheta > ele%value(angle$)/2 + pi) dtheta = dtheta - twopi

  local_position%r(1) = rho * sqrt_one(2*x/rho + (x/rho)**2 + (z/rho)**2)
  local_position%r(2) = y
  local_position%r(3) = dtheta * rho

  if (tilt == 0) then
    call rotate_mat(local_position%w, y_axis$, dtheta)
  else
    call rotate_mat(local_position%w, z_axis$, -tilt)
    call rotate_mat(local_position%w, y_axis$,  dtheta)
    call rotate_mat(local_position%w, z_axis$,  tilt)
  endif

  call update_floor_angles (local_position)

! patch case

elseif (ele%key == patch$) then
  if (ele%orientation == 1) then 
    if (logic_option(.false., use_patch_entrance)) then
      floor0 = ele%branch%ele(ele%ix_ele-1)%floor        ! Get floor0 from previous element
    else
      floor0 = ele%floor
    endif

  else
    if (logic_option(.false., use_patch_entrance)) then
      floor0 = ele%floor
    else
      floor0 = ele%branch%ele(ele%ix_ele+1)%floor        ! Get floor0 from next element
    endif
  endif

  local_position = coords_floor_to_relative (floor0, global_position)

  ! Is the particle inside the patch or outside?

  floor0 = ele%branch%ele(ele%ix_ele-1)%floor
  floor1 = ele%floor

  dz0 = ele%value(upstream_coord_dir$) * dot_product(floor0%w(:,3), (global_position%r - floor0%r)) ! Use w_inv = transpose
  dz1 = ele%value(downstream_coord_dir$) * dot_product(floor1%w(:,3), (global_position%r - floor1%r))

  if (dz0 > 0 .and. dz1 > 0) then
    status = downstream_end$
  elseif (dz0 < 0 .and. dz1 < 0) then
    status = upstream_end$
  else
    status = inside$
  endif

  return

! Straight line case

else
  floor0 = ele%floor
  local_position = coords_floor_to_relative (floor0, global_position)
  if (ele%orientation == 1) local_position%r(3) = local_position%r(3) + ele%value(l$)
endif

! Inside or outside?

if (local_position%r(3) < min(0.0_rp, ele%value(l$))) then
  if (ele%orientation == 1) then
    status = upstream_end$
  else
    status = downstream_end$
  endif
elseif (local_position%r(3) > max(0.0_rp, ele%value(l$))) then
  if (ele%orientation == 1) then
    status = downstream_end$
  else
    status = upstream_end$
  endif
else
  status = inside$
endif

! Optionally return w_mat

if (present(w_mat)) then
  w_mat = local_position%W
endif

end function coords_floor_to_local_curvilinear
