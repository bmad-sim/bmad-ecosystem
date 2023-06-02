!+
! Program photon_surface_data
!
! Program to output data about an element that can be used with gnuplot 
! to generate a 3D surface plot.
!
! Syntax:
!   > output_surface_data <param_file>
! Where
!   <param_file>  -- Name of the parameter file.
!-
 
program photon_surface_data

use photon_utils_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (photon_element_struct), pointer :: ph
type (ele_pointer_struct), allocatable :: eles(:)

real(rp) dr(2), r0(2), x, y, z, del_r(2), z0x, z0y, z1x, z1y

integer ie, n_arg, n_loc, ix, iy, granularity
integer ix_bounds(2), iy_bounds(2)

logical err

character(100) lat_file_name, param_file_name, data_file_name
character(40) ele_name

namelist / params / lat_file_name, ele_name, data_file_name, dr, r0, ix_bounds, iy_bounds, granularity

!

n_arg = command_argument_count()
if (n_arg /= 1) then
  print *, 'Syntax:'
  print *, '  > output_surface_data <param_file>'
  print *, ' Where'
  print *, '   <param_file>  -- Name of the parameter file.'
  stop
endif

!

r0 = real_garbage$
dr = real_garbage$
ix_bounds = int_garbage$
iy_bounds = int_garbage$
granularity = 1
data_file_name = ''
ele_name = ''
lat_file_name = ''

call get_command_argument(1, param_file_name)
open (1, file = param_file_name, status = 'old')
read (1, nml = params)
close (1)

if (data_file_name == '') data_file_name = 'surface.dat'

!

call bmad_parser(lat_file_name, lat)

if (ele_name == '') then
  ele => null()
  do ie = 1, lat%n_ele_max
    ele => lat%ele(ie)
    if (ele%key == crystal$) exit
  enddo

  if (.not. associated(ele)) then
    print *, 'Cannot find crystal element!'
    stop
  endif

else
  call lat_ele_locator(ele_name, lat, eles, n_loc)
  if (n_loc == 0) then
    print *, 'No elements found with name: ' // trim(ele_name)
    stop
  elseif (n_loc > 1) then
    print *, 'Note: Multiple elements matching name. Will use first one.'
  endif

  ele => eles(1)%ele
endif 

!

if (.not. associated(ele%photon)) then
  print *, 'No surface associated with element!'
  stop
endif

ph => ele%photon

if (any(ph%grid%dr /= 0)) then
  if (all(r0 == real_garbage$)) r0 =  ph%grid%r0
  if (all(dr == real_garbage$)) dr =  ph%grid%dr
  if (allocated(ph%grid%pt)) then
    if (all(ix_bounds == int_garbage$)) ix_bounds = [lbound(ph%grid%pt,1), ubound(ph%grid%pt,1)]
    if (all(iy_bounds == int_garbage$)) iy_bounds = [lbound(ph%grid%pt,2), ubound(ph%grid%pt,2)]
  endif
endif

if (r0(1) == int_garbage$) r0(1) = 0
if (r0(2) == int_garbage$) r0(2) = 0

if (any(dr == real_garbage$)) then
  print *, 'dr not set!'
  stop
endif

if (any(ix_bounds == int_garbage$) .or. any(iy_bounds == int_garbage$)) then
  print *, 'ix_bounds or iy_bounds not set!'
  stop
endif

!

ix_bounds = ix_bounds * granularity
iy_bounds = iy_bounds * granularity
dr = dr / granularity
del_r = 1d-6 * dr

!

open (1, file = data_file_name)
write (1, '(a)') '#     Ix      Iy          X                Y                Z              dZ/dX         dZ/dY'

do ix = ix_bounds(1), ix_bounds(2)
  write (1, *)
  x = r0(1) + ix * dr(1)
  do iy = iy_bounds(1), iy_bounds(2)
    y = r0(2) + iy * dr(2)
    z = z_at_surface (ele, x, y, err, .true.)
    z0x = z_at_surface (ele, x-del_r(1), y, err, .true.)
    z0y = z_at_surface (ele, x, y-del_r(2), err, .true.)
    z1x = z_at_surface (ele, x+del_r(1), y, err, .true.)
    z1y = z_at_surface (ele, x, y+del_r(2), err, .true.)
    write (1, '(2i8, 4x, 3es17.8, 2es14.6)') ix, iy, x, y, z, (z1x-z0x)/(2*del_r(1)), (z1y-z0y)/(2*del_r(2))
  enddo
enddo

close (1)
print '(2a)', 'Data file: ', trim(data_file_name)

end program
