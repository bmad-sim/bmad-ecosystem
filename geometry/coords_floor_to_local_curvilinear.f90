!+
! Function coords_floor_to_local_curvilinear  (global_position, ele, status, w_mat, relative_to_upstream) result(local_position)
!
! Given a position in global coordinates, return local curvilinear (laboratory) coordinates defined by ele.
!
! If the calculated longitudinal s places global_position outside of ele, the status argument
! will be set appropriately and local_position is not will be defined.
!
! Also see: coords_floor_to_curvilinear which will search a lattice branch for the element where
! the global_position is associated with.
!
! Input:
!   global_position     -- floor_position_struct: %r = [X, Y, Z] position in global coordinates
!   ele                 -- ele_struct: element to find local coordinates of.
!   relative_to_upstream  -- logical, optional: Default is False. If True, local_position%r(3) is relative to the 
!                              upstream end which will not be the entrance end if ele%orientation = -1.
!                              For patch elements (which always have ele%orientation = 1), if relative_to_upstream = T, 
!                              local_position%r(1:2) transverse coords will be taken wrt the upstream coords 
!                            
! Output:
!   local_position  -- floor_position_struct: %r = [x, y, z] position in local curvilinear coordinates
!                        with z relative to entrance edge of the element except if relative_to_upstream = T.
!   status          -- logical: longitudinal position:
!                               inside$: Inside the element.
!                               upstream_end$: Upstream of element.
!                               downstream_end$: Downstream of element.
!   w_mat(3,3)      -- real(rp) (optional): W matrix at s, to transform vectors. 
!                                  v_global = w_mat.v_local
!                                  v_local = transpose(w_mat).v_global
!-  

function coords_floor_to_local_curvilinear (global_position, ele, status, w_mat, relative_to_upstream) result(local_position)

use nr, only: zbrent
use bmad_interface, dummy => coords_floor_to_local_curvilinear

implicit none

type (floor_position_struct) :: global_position, local_position
type (ele_struct)   :: ele
type (floor_position_struct) :: floor0, floor1
real(rp), optional :: w_mat(3,3)
real(rp) x, y, z, g, dd, tilt, ds, dz0, dz1, wm(3,3), ref_tilt
integer :: status
logical, optional :: relative_to_upstream

! In all cases remember that ele%floor is the downstream floor coords independent of ele%orientation
! sbend case.

if (ele%key == sbend$) then
  g = ele%value(g$); ref_tilt = ele%value(ref_tilt_tot$)
  local_position = coords_floor_to_relative (ele%floor, global_position)

  if (g /= 0) then
    dd = (local_position%r(1) * cos(ref_tilt) + local_position%r(2) * sin(ref_tilt)) * g
    ds = atan(local_position%r(3) * g / (1 + dd)) / g
    local_position = bend_shift(local_position, g, ds, wm, ref_tilt)
    local_position%r(3) = ds
    if (present(w_mat)) w_mat = matmul(ele%floor%w, transpose(wm))
  else
    if (present(w_mat)) w_mat = ele%floor%w
  endif

  if (.not. (logic_option(.false., relative_to_upstream) .and. ele%orientation == -1)) then
    local_position%r(3) = local_position%r(3) + ele%value(l$)
  endif

  call update_floor_angles (local_position)

! patch case

elseif (ele%key == patch$) then
  if (ele%orientation == 1) then 
    if (logic_option(.false., relative_to_upstream)) then
      floor0 = ele%branch%ele(ele%ix_ele-1)%floor        ! Get floor0 from previous element
    else
      floor0 = ele%floor
    endif

  else
    if (logic_option(.false., relative_to_upstream)) then
      floor0 = ele%floor
    else
      floor0 = ele%branch%ele(ele%ix_ele+1)%floor        ! Get floor0 from next element
    endif
  endif

  local_position = coords_floor_to_relative (floor0, global_position)

  if (present(w_mat)) then
    w_mat = floor0%w
  endif

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

  if (present(w_mat)) then
    w_mat = floor0%w
  endif
endif

!----------------------
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

end function coords_floor_to_local_curvilinear
