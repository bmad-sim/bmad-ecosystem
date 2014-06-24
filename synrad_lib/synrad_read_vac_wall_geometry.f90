!+
! Subroutine synrad_read_vac_wall_geometry (wall_file, seg_len_max, s_lat, geometry, walls)
!
! Routine to read the vacuum wall geometry from two files: A file specifying the outline
! of the components used in constructing the machine and a file specifying where the components
! are located in the machine.
!
! Input:
!   wall_file      -- Character(*): Name of the wall file.
!   seg_len_max    -- Real(rp): Maximum length of wall segments.
!   s_lat          -- Real(rp): Lattice length
!   geometry       -- Integer: Type of lattice. open$ or closed$
!
! Output:
!   walls -- Walls_struct: wall structure.
!-

subroutine synrad_read_vac_wall_geometry (wall_file, seg_len_max, s_lat, geometry, walls)

use synrad_mod, except => synrad_read_vac_wall_geometry
use synrad3d_utils
use filename_mod

implicit none

type (walls_struct), target :: walls
type (wall_struct), pointer :: inside, outside
type (sr3d_wall_struct) wall3d

real(rp) s_lat, seg_len_max
real(rp) s, x_in, x_out

integer geometry, i, lun, ios, n_in, n_out

character(*) wall_file
character(40) name
character(200) file
character(*), parameter :: r_name = 'synrad_read_vac_wall_geometry'

logical phantom

namelist / wall_pt / s, x_in, x_out, name, phantom
 
! init

outside => walls%positive_x_wall
inside  => walls%negative_x_wall

outside%side = positive_x$
inside%side = negative_x$

! If a synrad3d file...

if (wall_file(1:10) == 'synrad3d::') then
  call sr3d_read_wall_file (wall_file(11:), s_lat, geometry, wall3d)
  call synrad3d_wall_to_synrad_walls (wall3d, seg_len_max, s_lat, geometry, walls)
  return
endif

! open file

lun = lunget()
call fullfilename (wall_file, file)
open (lun, file = file, status = 'old', iostat = ios)
if (ios /= 0) then
  call out_io (s_fatal$, r_name, 'CANNOT FIND FILE: ' // trim(wall_file))
  return
endif

! count number of wall points

n_in = -1; n_out = -1
do 
  x_in = real_garbage$; x_out = real_garbage$
  read (lun, nml = wall_pt, iostat = ios)
  if (ios > 0) then
    call out_io (s_fatal$, r_name, 'ERROR READING SYNRAD WALL FILE WALL_PT NAMELIST.')
    rewind (lun)
    do
      read (lun, nml = wall_pt) ! To generatge error message
    enddo
  endif
  if (ios < 0) exit
  if (x_in /= real_garbage$) n_in = n_in + 1
  if (x_out /= real_garbage$) n_out = n_out + 1
enddo

if (n_in == 0 .or. n_out == 0) then
  call out_io (s_fatal$, r_name, 'NO INSIDE AND/OR OUTSIDE WALL POINTS FOUND IN SYNRAD WALL FILE: ' // trim(wall_file))
  return
endif

allocate (inside%pt(0:n_in), outside%pt(0:n_out))
inside%n_pt_tot = n_in
outside%n_pt_tot = n_out

! Read wall file

rewind (lun)
n_in = -1; n_out = -1

do 
  x_in = real_garbage$; x_out = real_garbage$
  phantom = .false.
  read (lun, nml = wall_pt, iostat = ios)
  if (ios < 0) exit

  if (x_in /= real_garbage$) then
    n_in = n_in + 1
    inside%pt(n_in)%s = s
    inside%pt(n_in)%x = x_in
    inside%pt(n_in)%name = name
    inside%pt(n_in)%phantom = .false.
    inside%pt(n_in)%type = no_alley$
  endif

  if (x_out /= real_garbage$) then
    n_out = n_out + 1
    outside%pt(n_out)%s = s
    outside%pt(n_out)%x = x_out
    outside%pt(n_out)%name = name
    outside%pt(n_out)%phantom = phantom
    outside%pt(n_out)%type = possible_alley$
  endif
enddo

!

call delete_overlapping_wall_points (outside)
call delete_overlapping_wall_points (inside)

forall (i = 0:n_out) outside%pt(i)%ix_pt = i
forall (i = 0:n_in)  inside%pt(i)%ix_pt = i

call create_alley (inside)
call create_alley (outside)

! check that endpoints are correct

if (abs(outside%pt(outside%n_pt_tot)%s - s_lat) > 0.01) then
  print *, 'Note: Outside wall ends at:', outside%pt(outside%n_pt_tot)%s
  print *, '      And not at lattice end of:', s_lat
  print *, '      [But last point is always adjusted to have s = s_lat]'
endif

if (abs(inside%pt(inside%n_pt_tot)%s - s_lat) > 0.01) then
  print *, 'Note: Inside wall ends at:', inside%pt(inside%n_pt_tot)%s
  print *, '      And not at lattice end of:', s_lat
  print *, '      [But last point is always adjusted to have s = s_lat]'
endif

outside%pt(outside%n_pt_tot)%s = s_lat
inside%pt(inside%n_pt_tot)%s = s_lat

! do some checking

call check_wall (inside, s_lat, geometry)
call check_wall (outside, s_lat, geometry)

! segment wall

call break_wall_into_segments (inside, seg_len_max)
call break_wall_into_segments (outside, seg_len_max)

end subroutine
