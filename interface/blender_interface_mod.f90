module blender_interface_mod

use geometry_mod

implicit none

private skip_ele_blender, write_blender_ele

contains

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine write_blender_lat_layout(iu, lat)

type(lat_struct) :: lat
type(ele_struct), pointer :: ele
integer :: iu, n, ix_ele

! Header

write(iu, '(a)') 'ele_name, ix_ele, x, y, z, theta ,phi, psi, key, L, custom1, custom2, custom3, descrip'

do n = 0, ubound(lat%branch, 1)
  do ix_ele = 0, lat%branch(n)%n_ele_max

    ele => lat%branch(n)%ele(ix_ele)
    if (skip_ele_blender(ele)) cycle
    call write_blender_ele(iu, ele)

  enddo
enddo

end subroutine write_blender_lat_layout

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

function skip_ele_blender (ele) result (skip)

type(ele_struct) :: ele
logical :: skip

! Skip slaves

skip = .true. 

if (ele%slave_status == multipass_slave$ .or. ele%slave_status == super_slave$) return
select case (ele%key)
case (beginning_ele$, patch$, marker$, match$, null_ele$, floor_shift$, group$, overlay$)
  return
end select

skip = .false.

end function skip_ele_blender

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------


subroutine write_blender_ele(iu, ele)

type(ele_struct) :: ele
type(floor_position_struct) :: local_position, p
type (coord_struct) orbit

real(rp) :: w_mat(3,3), w2_mat(3,3)
integer :: iu

character(2), parameter :: c = ', '
character(500) line

!

local_position = floor_position_struct()
local_position%r = [0.0_rp, 0.0_rp, ele%value(L$)/2]

p = coords_local_curvilinear_to_floor (local_position, ele, in_ele_frame=.true., w_mat=w_mat)

select case (ele%key)
case (crystal$, mirror$, multilayer_mirror$)
  orbit%vec = 0
  call init_coord (orbit, orbit%vec, ele, upstream_end$, photon$)
  call offset_photon (ele, orbit, unset$, rot_mat = w2_mat)
  call mat_inverse (w2_mat, w2_mat)
  w_mat = matmul(w_mat, w2_mat)
end select

call floor_w_mat_to_angles (w_mat, 0.0_rp, p%theta, p%phi, p%psi)

write(line, '(2a, i0, 6(es14.6, a), 2a, es14.6, a)') trim(ele%name), c,  ele%ix_ele, c, &
                    p%r(1), c, p%r(2), c, p%r(3), c, p%theta, c, p%phi, c, p%psi, c, &
                    trim(key_name(ele%key)), c, ele%value(l$)

select case (ele%key)
case (sbend$)
  write (line, '(a, 3(a, es14.6))') trim(line), c, ele%value(angle$), c, ele%value(e1$), c, ele%value(e2$)
case (crystal$, detector$, mirror$, multilayer_mirror$, diffraction_plate$)
  write (line, '(a, 3(a, es14.6))') trim(line), c, ele%value(x_limit$), c, ele%value(y_limit$), c, 1e-3_rp
case default
  write (line, '(2a)') trim(line), ', 0.0, 0.0, 0.0'
end select

if (associated(ele%descrip)) then
  write (line, '(4a)') trim(line), ', "', ele%descrip, '"'
else
  write (line, '(2a)') trim(line), ', ""'
endif

!

write (iu, '(a)') trim(line)

end subroutine write_blender_ele

end module

