!+
! This module is for creating a Cubit script.
! This module is not used. 
! 3D modeling is now done with the Blender program.
!-

module tao_write_3d_mod

use tao_utils
implicit none

private ele_to_3d_shape, ele_draw_line, write_head_stuff, write_foot_stuff
private ele_draw_limit_box, construct_cubit_name, draw_box, orient_body

integer, private, save :: iu

contains

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!+
! Subroutine tao_write_3d_floor_plan (floor_file, lat)
!
! Routine to write a 3D floor plan file for the 3D_floor_plan.jou script.
!
! Input:
!   floor_file -- character(*): Name of output file.
!   lat        -- lat_struct: Lattice to output.
!-

subroutine tao_write_3d_floor_plan (floor_file, lat)

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (tao_ele_shape_struct), pointer :: ele_shape

character(*) floor_file

integer ib, ie, ios
logical shape_drawn

character(*), parameter :: r_name = 'tao_write_3d_floor_plan'

!

iu = lunget()
open (iu, file = floor_file, iostat = ios)
if (ios /= 0) then
  print *, 'CANNOT OPEN FILE: ', trim(floor_file)
  return
endif

call write_head_stuff()


do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do ie = 1, lat%n_ele_max
    ele => branch%ele(ie)

    call ele_to_3d_shape (ele, shape_drawn)
    if (.not. shape_drawn) call ele_draw_line (ele)

  enddo
enddo

call out_io (s_info$, r_name, 'Written: ' // floor_file)

end subroutine tao_write_3d_floor_plan

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine ele_to_3d_shape (ele, shape_drawn)

type (ele_struct) ele
type (tao_ele_shape_struct), pointer :: ele_shape
logical shape_drawn

!

shape_drawn = .false.

ele_shape => tao_pointer_to_ele_shape (ele, s%plot_page%floor_plan%ele_shape)
if (.not. associated(ele_shape)) return
if (.not. ele_shape%draw) return
if (ele%value(x1_limit$) == 0) return
!select case (ele%key)
!case (diffraction_plate$)

call ele_draw_limit_box (ele, ele_shape)

! end select


shape_drawn = .true.

end subroutine ele_to_3d_shape

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine ele_draw_limit_box (ele, ele_shape)

type (ele_struct) ele
type (tao_ele_shape_struct) ele_shape
type (coord_struct) orbit

real(rp) r0(3), r1(3)
real(rp) w_mat(3,3), axis(3), angle

character(40) cubit_name

! 

r0 = [-ele%value(x1_limit$), -ele%value(y1_limit$), 0.0_rp]
r1 = [ ele%value(x2_limit$),  ele%value(y2_limit$), max(ele%value(l$), 0.001_rp)]

call construct_cubit_name (ele, cubit_name)
call draw_box (r0, r1, ele_shape%color, cubit_name)

if (ele%aperture_at == surface$) then
  orbit%vec = 0
  call init_coord (orbit, orbit%vec, ele, upstream_end$, photon$)
  call offset_photon (ele, orbit, unset$, rot_mat = w_mat)
  call w_mat_to_axis_angle (w_mat, axis, angle)
  if (angle /= 0) then
    write (iu, '(3a, f9.2, a, 3f11.6, a)') &
            'cubit.cmd("rotate ', trim(cubit_name), ' angle ', 180*angle/pi, ' about origin 0 0 0 direction', axis, '")'
  endif
endif

call orient_body (cubit_name, ele%floor)

end subroutine ele_draw_limit_box

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine construct_cubit_name (ele, cubit_name)

type (ele_struct) ele
character(*) cubit_name

!

write (cubit_name, '(2a, i0)') trim(ele%name), '_', ele%ix_ele

end subroutine construct_cubit_name

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine draw_box (r0, r1, color, cubit_name)

real(rp) r0(3), r1(3)

character(*) color, cubit_name

! 

write (iu, '(a)')  '#'
write (iu, '(2a)') '# Draw: ', trim(cubit_name)
write (iu, '(4(a, f12.6))') &
                  'cubit.cmd("brick x', r1(1)-r0(1), ' y', &
                                        r1(2)-r0(2), ' z', r1(3) - r0(3), '")'
write (iu, '(a)')  'id = cubit.get_last_id("volume")'
write (iu, '(5a)') 'cubit.cmd("Volume " + str(id) + " rename ', "'", trim(cubit_name), "'", '")'
write (iu, '(3a, f12.6, a)') &
                  'cubit.cmd("move ', trim(cubit_name), ' x 0 y 0 z ', (r1(3) - r0(3)) / 2, '")'

end subroutine draw_box

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine orient_body (cubit_name, floor)

type (floor_position_struct) floor
character(*) cubit_name

real(rp) w_mat(3,3), axis(3), angle

! Rotate

call floor_angles_to_w_mat (floor%theta, floor%phi, floor%psi, w_mat)
call w_mat_to_axis_angle (w_mat, axis, angle)

if (angle /= 0) then
  write (iu, '(3a, f9.2, a, 3f11.6, a)') &
          'cubit.cmd("rotate ', trim(cubit_name), ' angle ', 180*angle/pi, ' about origin 0 0 0 direction', axis, '")'
endif

! Translate

write (iu, '(3a, 3(f12.6, a))') &
          'cubit.cmd("move ', trim(cubit_name), ' x', floor%r(1), ' y', floor%r(2), ' z', floor%r(3), '")'

end subroutine orient_body

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine ele_draw_line (ele)

type (ele_struct) ele

! 

end subroutine ele_draw_line

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine write_head_stuff ()

write (iu, '(a)') '#!python'
write (iu, '(a)') ''
write (iu, '(a)') 'cubit.cmd("undo on")'
write (iu, '(a)') ''

end subroutine write_head_stuff

!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------

subroutine write_foot_stuff ()

write (iu, '(a)') ''

end subroutine write_foot_stuff

end module
