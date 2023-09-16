module blender_interface_mod

use bmad_interface

implicit none

!private skip_ele_blender, write_blender_ele

contains

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

subroutine write_blender_lat_layout(file_name, lat)

type(lat_struct) :: lat
type(ele_struct), pointer :: ele
integer :: iu, n, ix_ele, ix
character(*) file_name
character(200) basename

! Header

iu = lunget()
open (iu, file = file_name)

write (iu, '(a)') "import os"
write (iu, *)
write (iu, '(a)') "def map_table_dict(line, des):"
write (iu, '(a)') "    d = {}"
write (iu, '(a)') "    vals = line.split(',')[0:14]"
write (iu, '(a)') "    #d['layer']  = vals[0].strip()"
write (iu, '(a)') "    d['name']   = vals[0].strip()"
write (iu, '(a)') "    d['index']  = int(vals[1])"
write (iu, '(a)') "    d['x']      = float(vals[2])"
write (iu, '(a)') "    d['y']      = float(vals[3])"
write (iu, '(a)') "    d['z']      = float(vals[4])"
write (iu, '(a)') "    d['theta']  = float(vals[5])"
write (iu, '(a)') "    d['phi']    = float(vals[6])"
write (iu, '(a)') "    d['psi']    = float(vals[7])"
write (iu, '(a)') "    d['key']    = vals[8].strip()"
write (iu, '(a)') "    d['L']      = float(vals[9])  "
write (iu, '(a)') "    if d['key']=='SBEND':"
write (iu, '(a)') "        d['angle'] = float(vals[10])"
write (iu, '(a)') "        d['e1'] = float(vals[11])"
write (iu, '(a)') "        d['e2'] = float(vals[12])"
write (iu, '(a)') "    d['descrip'] = des"
write (iu, '(a)') "    return d"
write (iu, *)
write (iu, '(a)') 'lat = []'
write (iu, '(a)') '# ele_name, ix_ele, x, y, z, theta ,phi, psi, key, L, custom1, custom2, custom3, descrip'
write (iu, *)

! Elements

do n = 0, ubound(lat%branch, 1)
  do ix_ele = 0, lat%branch(n)%n_ele_max

    ele => lat%branch(n)%ele(ix_ele)
    if (skip_ele_blender(ele)) cycle
    call write_blender_ele(iu, ele)

  enddo
enddo

! Footer

basename = file_name
ix = str_last_in_set (basename, '.')
if (ix /= 0) basename = basename(:ix-1)

write (iu, *)
write (iu, '(a)') "bmad_base_dir = os.getenv('ACC_ROOT_DIR')"
write (iu, '(a)') "basename = '" // trim(basename) // "'"
write (iu, '(a)') "template_file = bmad_base_dir + '/bmad/scripts/blender_base.py'"
write (iu, '(a)') "exec(open(template_file).read())"

close(iu)

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
case (beginning_ele$, patch$, match$, null_ele$, floor_shift$, &
      group$, overlay$, fork$, photon_fork$)
  return
end select

skip = .false.

end function skip_ele_blender

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------

 
subroutine write_blender_ele(iu, ele, old_format)

type(ele_struct) :: ele
type(floor_position_struct) :: local_position, p
type (coord_struct) orbit

real(rp) :: w_mat(3,3), w2_mat(3,3)
integer :: iu

logical, optional :: old_format

character(2), parameter :: c = ', '
character(200) des
character(500) line

!

call mat_make_unit(local_position%w)
local_position%r = [0.0_rp, 0.0_rp, ele%value(L$)/2]

p = coords_local_curvilinear_to_floor (local_position, ele, in_body_frame=.true., w_mat=w_mat)

select case (ele%key)
case (crystal$, mirror$, multilayer_mirror$)
  orbit%vec = 0
  call init_coord (orbit, orbit%vec, ele, upstream_end$, photon$)
  call offset_photon (ele, orbit, unset$, rot_mat = w2_mat)
  call mat_inverse (w2_mat, w2_mat)
  w_mat = matmul(w_mat, w2_mat)
end select

call floor_w_mat_to_angles (w_mat, p%theta, p%phi, p%psi)



write(line, '(2a, i0, a, 6(es16.8, a), 2a, es16.8, a)') trim(ele%name), c,  ele%ix_ele, c, &
                    p%r(1), c, p%r(2), c, p%r(3), c, p%theta, c, p%phi, c, p%psi, c, &
                    trim(key_name(ele%key)), c, ele%value(l$)

select case (ele%key)
case (pipe$)
  write (line, '(a, 3(a, es16.8))') trim(line), c, ele%value(x1_limit$), c, ele%value(y1_limit$), c, 0.002_rp
case (sbend$)
  write (line, '(a, 3(a, es16.8))') trim(line), c, ele%value(angle$), c, ele%value(e1$), c, ele%value(e2$)
case (crystal$, detector$, mirror$, multilayer_mirror$, diffraction_plate$, mask$)
  write (line, '(a, 3(a, es16.8))') trim(line), c, ele%value(x1_limit$), c, ele%value(y1_limit$), c, 1e-3_rp
case (wiggler$)
  write (line, '(a, 3(a, es16.8))') trim(line), c, ele%value(x1_limit$), c, ele%value(y1_limit$), c
case default
  write (line, '(2a)') trim(line), ', 0.0, 0.0, 0.0'
end select

des = ''
if (associated(ele%descrip)) des = ele%descrip

!
if (logic_option(.false., old_format)) then
  write (iu, '(4a)') trim(line), ', "', trim(des), '"'
else
  write (iu, '(5a)') 'lat.append(map_table_dict("', trim(line), '", "', trim(des), '"))'
endif

end subroutine write_blender_ele

end module

