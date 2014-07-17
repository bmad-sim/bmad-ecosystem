!+
! Subroutine synrad_read_vac_wall_geometry (wall_file, seg_len_max, branch, walls)
!
! Routine to read the vacuum wall geometry from two files: A file specifying the outline
! of the components used in constructing the machine and a file specifying where the components
! are located in the machine.
!
! Input:
!   wall_file      -- Character(*): Name of the wall file.
!   seg_len_max    -- Real(rp): Maximum length of wall segments.
!   branch         -- branch_struct: lattice branch to use
!
! Output:
!   walls -- Walls_struct: wall structure.
!-

subroutine synrad_read_vac_wall_geometry (wall_file, seg_len_max, branch, walls)

use synrad_mod, except => synrad_read_vac_wall_geometry
use synrad3d_utils
use filename_mod

implicit none

type (walls_struct), target :: walls
type (branch_struct) branch
type (wall_struct), pointer :: inside, outside
type (sr3d_wall_struct) wall3d

real(rp) seg_len_max
real(rp) s, x_in, x_out

integer i, lun, ios, n_in, n_out

character(*) wall_file
character(40) name
character(200) file
character(*), parameter :: r_name = 'synrad_read_vac_wall_geometry'

logical phantom

namelist / wall_pt / s, x_in, x_out, name, phantom
 
! If a synrad3d file...

if (wall_file(1:10) == 'synrad3d::') then
  call sr3d_read_wall_file (wall_file(11:), branch%param%total_length, branch%param%geometry, wall3d)
  call synrad3d_wall_to_synrad_walls (wall3d, seg_len_max, branch, walls)
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

outside => walls%positive_x_wall
inside  => walls%negative_x_wall

allocate (inside%pt(0:n_in), outside%pt(0:n_out))
inside%n_pt_max = n_in
outside%n_pt_max = n_out

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
  endif

  if (x_out /= real_garbage$) then
    n_out = n_out + 1
    outside%pt(n_out)%s = s
    outside%pt(n_out)%x = x_out
    outside%pt(n_out)%name = name
    outside%pt(n_out)%phantom = phantom
  endif
enddo

!

call synrad_setup_walls (walls, branch, seg_len_max)

end subroutine
