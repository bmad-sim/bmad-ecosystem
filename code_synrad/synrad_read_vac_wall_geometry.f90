!+
! Subroutine synrad_read_vac_wall_geometry (wall_file, seg_len_max, branch, walls, err_flag, seg_len_phantom_max)
!
! Routine to read the vacuum wall geometry from two files: A file specifying the outline
! of the components used in constructing the machine and a file specifying where the components
! are located in the machine.
!
! Input:
!   wall_file           -- character(*): Name of the wall file.
!   seg_len_max         -- real(rp): Maximum length of wall segments.
!   branch              -- branch_struct: lattice branch to use
!   seg_len_phantom_max -- real(rp), optional: If present then use this number for phantom segments
!                             instead of seg_len_max.
!
! Output:
!   walls     -- walls_struct: wall structure.
!   err_flag  -- logical, optional: Set true if there is a problem
!-

subroutine synrad_read_vac_wall_geometry (wall_file, seg_len_max, branch, walls, err_flag, seg_len_phantom_max)

use synrad_mod, except => synrad_read_vac_wall_geometry
use synrad3d_parse_wall
use synrad3d_wall_to_synrad_walls_mod, except2 => synrad_read_vac_wall_geometry
use filename_mod

implicit none

type (walls_struct), target :: walls
type (branch_struct) branch
type (wall_struct), pointer :: minus_side, plus_side

real(rp), optional :: seg_len_phantom_max
real(rp) seg_len_max
real(rp) s, x_in, x_out, x_plus, x_minus

integer i, lun, ios, n_minus, n_plus

character(*) wall_file
character(40) name
character(200) file
character(*), parameter :: r_name = 'synrad_read_vac_wall_geometry'

logical, optional :: err_flag
logical phantom, err

namelist / wall_pt / s, x_in, x_out, x_plus, x_minus, name, phantom
 
!

if (present(err_flag)) err_flag = .true.

! If a synrad3d file...

if (wall_file(1:10) == 'synrad3d::') then
  call sr3d_read_wall_file (wall_file(11:), branch%lat, err)
  if (err) return
  call synrad3d_wall_to_synrad_walls (branch, seg_len_max, walls)
  if (present(err_flag)) err_flag = .false.
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

n_minus = -1; n_plus = -1
do 
  x_in = real_garbage$; x_out = real_garbage$
  x_minus = real_garbage$; x_plus = real_garbage$
  read (lun, nml = wall_pt, iostat = ios)
  if (ios > 0) then
    call out_io (s_fatal$, r_name, 'ERROR READING SYNRAD WALL FILE WALL_PT NAMELIST.')
    rewind (lun)
    do
      read (lun, nml = wall_pt) ! To generatge error message
    enddo
  endif
  if (ios < 0) exit
  if (x_in /= real_garbage$) x_minus = x_in
  if (x_out /= real_garbage$) x_plus = x_out
  if (x_minus /= real_garbage$) n_minus = n_minus + 1
  if (x_plus /= real_garbage$) n_plus = n_plus + 1
enddo

if (n_minus == 0 .or. n_plus == 0) then
  call out_io (s_fatal$, r_name, 'NO X_PLUS AND/OR X_MINUS WALL POINTS FOUND IN SYNRAD WALL FILE: ' // trim(wall_file))
  return
endif

plus_side => walls%positive_x_wall
minus_side  => walls%negative_x_wall

allocate (minus_side%pt(0:n_minus), plus_side%pt(0:n_plus))
minus_side%n_pt_max = n_minus
plus_side%n_pt_max = n_plus

! Read wall file

rewind (lun)
n_minus = -1; n_plus = -1

do 
  x_in = real_garbage$; x_out = real_garbage$
  x_minus = real_garbage$; x_plus = real_garbage$
  phantom = .false.
  read (lun, nml = wall_pt, iostat = ios)
  if (ios < 0) exit

  if (x_in /= real_garbage$) x_minus = x_in
  if (x_out /= real_garbage$) x_plus = x_out

  if (x_minus /= real_garbage$) then
    n_minus = n_minus + 1
    minus_side%pt(n_minus)%s = s
    minus_side%pt(n_minus)%x = x_minus
    minus_side%pt(n_minus)%name = name
    minus_side%pt(n_minus)%phantom = .false.
  endif

  if (x_plus /= real_garbage$) then
    n_plus = n_plus + 1
    plus_side%pt(n_plus)%s = s
    plus_side%pt(n_plus)%x = x_plus
    plus_side%pt(n_plus)%name = name
    plus_side%pt(n_plus)%phantom = phantom
  endif
enddo

!

call synrad_setup_walls (walls, branch, seg_len_max, seg_len_phantom_max)

if (present(err_flag)) err_flag = .false.

end subroutine
