!+
! Program construct_taylor_map
!
! Program to construct the truncated Taylor series map for a lattice.
! Optionally, including spin.
!
! Usage:
!   <path_to_bin_dir>/construct_taylor_map {-spin} {-from <ele_name_or_num>} <lattice_file_name>
!
! Where:
!   -spin                    -> Include spin part of map.
!   -from <ele_name_or_num>  -> Specify start/stop location of map.
!-

program construct_taylor_map

use bmad_interface

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_pointer_struct), allocatable :: eles(:)
type (taylor_struct) xfer_map(6), spin_map(4)
type (coord_struct), allocatable :: orb0(:)

integer i, ix, n_loc, ixe

logical ok, include_spin, err

character(40) from_ele
character(100) arg, lat_file, out_file

! Read in command line

include_spin = .false.
from_ele = 'beginning'
lat_file = ''
out_file = ''
ok = .true.

i = 0
do while (i < command_argument_count())
  i=i+1; call get_command_argument(i, arg)
  select case (arg)
  case ('-spin')
    include_spin = .true.
  case ('-from')
    i=i+1; i=i+1; call get_command_argument(i, from_ele)
  case ('-out')
    i=i+1; i=i+1; call get_command_argument(i, out_file)
  case default
    if (arg(1:1) == '-') then
      print *, 'I DO NOT UNDERSTAND: ', trim(arg)
      print *
      ok = .false.
    endif
    lat_file = arg
  end select
enddo

!

if (.not. ok .or. lat_file == '') then
  print '(a)', 'Usage:'
  print '(a)', '  <bin_dir>/construct_taylor_map {-spin} {-from <ele_name_or_num>} {-out <out_file>} <lat_file>'
  print '(a)', ''
  print '(a)', 'Where:'
  print '(a)', '  -spin                    -> Include spin part of map.'
  print '(a)', '  -from <ele_name_or_num>  -> Exit end of element is start/stop location of map.'
  print '(a)', '  -out <out_file>          -> Specify output data file. Default is <lat_file>.map'
  print '(a)', ''
  print '(a)', 'Output:'
  stop
endif

!

call bmad_parser (lat_file, lat)

call lat_ele_locator (from_ele, lat, eles, n_loc)
if (n_loc > 1) then
  print '(a)', 'More than one element matches: ' // trim(from_ele)
  print '(a)', 'Nothing done.'
  stop
endif

if (n_loc == 0) then
  print '(a)', 'No element match: ' // trim(from_ele)
  print '(a)', 'Nothing done.'
  stop
endif

!

ele => eles(1)%ele
ixe = ele%ix_ele

call closed_orbit_calc (lat, orb0, ix_branch = ele%ix_branch)

!

if (out_file == '') out_file = trim(lat_file) // '.map'
open (1, file = out_file)

if (include_spin) then
  call ptc_transfer_map_with_spin (ele%branch, xfer_map, spin_map, orb0(ixe), err, ixe, ixe, .true.)
  call type_taylors (xfer_map, file_id = 1)
  call type_taylors(spin_map, file_id = 1)

else
  ix = ele%ix_ele
  call transfer_map_calc (lat, xfer_map, err, ix, ix, orb0(ix), ele%ix_branch, .true., .true.)
  call type_taylors (xfer_map, file_id = 1)
endif

print '(a)', 'Output file: ', trim(out_file)

end program
