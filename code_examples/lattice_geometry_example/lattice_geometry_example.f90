program lattice_geometry_example

use beam_mod
use bmad

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele

real(rp) :: ds, sc

integer :: i_dim, i, j, iu
integer :: namelist_file, n_char

character(100) :: lat_name, lat_path, base_name, in_file
character(30), parameter :: r_name = 'lattice_geometry_example'

namelist / lattice_geometry_example_params / &
    lat_name

!------------------------------------------
!Defaults for namelist
lat_name = 'lat.bmad'

!Read namelist
in_file = 'lattice_geometry_example.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)

namelist_file = lunget()
!print *, 'Opening: ', trim(in_file)
open (namelist_file, file = in_file, status = "old")
read (namelist_file, nml = lattice_geometry_example_params)
close (namelist_file)

!Trim filename
n_char= SplitFileName(lat_name, lat_path, base_name) 

!Parse Lattice
call bmad_parser (lat_name, lat)
!branch => lat%branch(0)

!--------------------
sc = 0.1

iu=lunget()
call type_header()
open(iu, file='reference_frame.dat')
call type_lat_floor(iu, 0.0_rp, 0.0_rp, in_body_frame=.false.)
close(iu)
call type_header()
open(iu, file='element_frame.dat')
call type_lat_floor(iu, 0.0_rp, 0.0_rp, in_body_frame=.true.)
close(iu)


open(iu, file='reference_boxes.dat')
call type_header()
call type_lat_floor(iu, 1*sc,  1*sc, in_body_frame=.false.)
call type_lat_floor(iu, -1*sc, 1*sc, in_body_frame=.false.)
call type_lat_floor(iu, -1*sc,-1*sc, in_body_frame=.false.)
call type_lat_floor(iu,  1*sc,-1*sc, in_body_frame=.false.)
close(iu)

open(iu, file='element_boxes.dat')
call type_header()
call type_lat_floor(iu, 1*sc,  1*sc, in_body_frame=.true.)
call type_lat_floor(iu, -1*sc, 1*sc, in_body_frame=.true.)
call type_lat_floor(iu, -1*sc,-1*sc, in_body_frame=.true.)
call type_lat_floor(iu,  1*sc,-1*sc, in_body_frame=.true.)
close(iu)


contains

subroutine type_lat_floor(iu, x, y, in_body_frame)
real(rp) :: x, y
integer :: i, iu
logical, optional :: in_body_frame
do i=1, lat%n_ele_track
  ele=>lat%ele(i)
  call type_floor(iu, ele, x,  y, in_body_frame = logic_option(.false., in_body_frame))
enddo
end subroutine


subroutine type_header()
write(iu,'(12a10)') 'x', 'y', 'z', 'Xn_x', 'Xn_y', 'Xn_z', 'Yn_x', 'Yn_y', 'Yn_z', 'Zn_x', 'Zn_y', 'Zn_z' 
end subroutine
                                 
 
subroutine type_floor(iu, ele, x, y, in_body_frame)
type(floor_position_struct) :: f0, f1, f2
type(ele_struct) :: ele
real(rp) :: ds, s, w_mat(3,3), x, y
integer i, iu, n_steps
logical, optional :: in_body_frame

n_steps = 20
ds = ele%value(L$)/(n_steps-1)
if (ds==0) n_steps=0


do i=0, n_steps-1
  s = i*ds
  f0 = floor_position_struct(vec3_zero$, mat3_unit$, 0.0_rp, 0.0_rp, 0.0_rp)
  f0%r = [x, y, s]
  f1= coords_local_curvilinear_to_floor (f0, ele, in_body_frame = logic_option(.false., in_body_frame), w_mat = w_mat)
  write(iu,'(12f10.5)') f1%r, w_mat(:,1), w_mat(:,2), w_mat(:,3)
enddo

end subroutine





end program
