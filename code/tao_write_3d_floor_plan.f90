!+
! Subroutine tao_write_3d_floor_plan (floor_file, lat)
!
! Routine to write a 3D floor plan file for the 3D_floor_plan.py script.
!
! Input:
!   floor_file -- character(*): Name of output file.
!   lat        -- lat_struct: Lattice to output.
!-

subroutine tao_write_3d_floor_plan (floor_file, lat)

use tao_utils

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch
type (tao_ele_shape_struct), pointer :: ele_shape

character(*) floor_file

integer ib, ie, ios, iu
logical shape_drawn

character(*), parameter :: r_name = 'tao_write_3d_floor_plan'

!

print *, '3D Floor Plan not yet implemented.'
return

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
    if (.not. shape_drawn) call ele_to_line (ele)

  enddo
enddo

call out_io (s_info$, r_name, 'Written: ' // floor_file)

!----------------------------------------------------------------
contains

subroutine ele_to_3d_shape (ele, shape_drawn)

type (ele_struct) ele

logical shape_drawn

!

shape_drawn = .false.

ele_shape => tao_pointer_to_ele_shape (ele, s%plotting%floor_plan%ele_shape)
if (.not. associated(ele_shape)) return
if (.not. ele_shape%draw) return



write (1, '(a)') 'draw_box (, )'

shape_drawn = .true.

end subroutine ele_to_3d_shape

!----------------------------------------------------------------
! contains

subroutine ele_to_line (ele)

type (ele_struct) ele

! 

end subroutine ele_to_line


!----------------------------------------------------------------
! contains

subroutine write_head_stuff ()

write (1, '(a)') 'import matplotlib.pyplot as plt'
write (1, '(a)') 'from mpl_toolkits.mplot3d import Axes3D'
write (1, '(a)') 'import numpy as np'
write (1, '(a)')
write (1, '(a)') 'def draw_box (r0, w0):'
write (1, '(a)') ''
write (1, '(a)') ''
write (1, '(a)') ''
write (1, '(a)') 
write (1, '(a)') 'fig = plt.figure()'
write (1, '(a)') 'ax = plt.Axes3D(fig)'

end subroutine write_head_stuff

!----------------------------------------------------------------
! contains

subroutine write_foot_stuff ()

write (1, '(a)') 'plt.show()'

end subroutine write_foot_stuff

end subroutine tao_write_3d_floor_plan
