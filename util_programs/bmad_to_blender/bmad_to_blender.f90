program bmad_to_blender

use blender_interface_mod
use beam_mod



type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) :: orb0
type (coord_struct), allocatable :: closed_orb(:)

integer :: i_dim, ios, n, ix_ele
integer :: namelist_file, outfile, n_char

character(100) :: lat_name, lat_path, base_name, in_file, outfile_name
character(30) :: suffix
character(30), parameter :: r_name = 'bmad_to_blender'


namelist / bmad_to_blender_params / &
    lat_name, suffix

!------------------------------------------
!Defaults for namelist
lat_name = 'lat.bmad'
suffix = 'layout_table'


!Read namelist
in_file = 'bmad_to_blender.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)

namelist_file = lunget()
print *, 'Opening: ', trim(in_file)
open (namelist_file, file = in_file, status = "old")
read (namelist_file, nml = bmad_to_blender_params, iostat=ios)
close (namelist_file)
if (ios > 0) then
  ! Error in namelist
  close (namelist_file)
else if (ios < 0) then
  ! No namelist found. Use 
  lat_name = in_file
endif





!Trim filename
n_char= SplitFileName(lat_name, lat_path, base_name) 

! Prepare outfile_name
call file_suffixer (lat_name, outfile_name, trim(suffix), .true.)

!Parse Lattice
call bmad_parser (lat_name, lat)
!branch => lat%branch(0)

outfile = lunget()
open (outfile, file = outfile_name)
!call write_blender_lat_layout(outfile, lat)
!call write_blender_lat_layout(outfile_name, lat)

! Old format
write (outfile, '(a)') '# ele_name, ix_ele, x, y, z, theta ,phi, psi, key, L, custom1, custom2, custom3, descrip'

! Elements

do n = 0, ubound(lat%branch, 1)
  do ix_ele = 0, lat%branch(n)%n_ele_max
    ele => lat%branch(n)%ele(ix_ele)
    if (skip_ele_blender(ele)) cycle
    call write_blender_ele(outfile, ele, old_format = .true.)

  enddo
enddo


write(*, '(2a)') trim(outfile_name), ' written'
close(outfile)

end program
