program synrad3d_test

use synrad3d_track_mod
use photon_reflection_mod

implicit none

type (lat_struct), target :: lat
type (sr3d_wall_struct), target :: wall
type (sr3d_photon_track_struct) :: photon
type (sr3d_photon_coord_struct) p
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)

real(rp) vel
integer ios, num_ignored, n_photon

logical is_inside, err, absorbed

character(100) old_wall_file, wall_file

namelist / in / p, wall_file

! Init

call bmad_parser('lat.bmad', lat)

sr3d_params%specular_reflection_only = .true.
sr3d_params%allow_absorption = .false.
num_ignored = 0
allocate (wall_hit(0:1))

! Specular reflection test...

old_wall_file = 'xxx' 
open (1, file = 'output.now')
open (2, file = 'specular.input')

n_photon = 0
do 

  read (2, nml = in, iostat = ios)
  if (ios < 0) exit 
  if (ios > 0) then
    print *, 'Error reading photon starting position at photon index:', n_photon
    call err_exit
  endif

  if (wall_file /= old_wall_file) then
    if (allocated(wall%section)) deallocate (wall%section)
    if (allocated(wall%gen_shape)) deallocate (wall%gen_shape)
    call sr3d_read_wall_file (wall_file, lat%ele(lat%n_ele_track)%s, lat%param%geometry, wall)
    old_wall_file = wall_file
    cycle
  endif

  vel = sqrt(p%vec(2)**2 + p%vec(4)**2 + p%vec(6)**2)
  if (abs(vel - 1) > 0.1) then
    print *, 'ERROR: PHOTON VELOCITY NOT PROPERLY NORMALIZED TO 1 FOR PHOTON:', n_photon
    stop
  endif
  p%vec(2:6:2) = p%vec(2:6:2) / vel

  p%energy = 1000             ! Arbitrary
  p%ix_ele = element_at_s(lat, p%vec(5), .true.)
  photon%start = p
  photon%n_wall_hit = 0

  n_photon = n_photon + 1
  photon%ix_photon = n_photon
  photon%ix_photon_generated = n_photon

  call sr3d_check_if_photon_init_coords_outside_wall (p, wall, is_inside, num_ignored)

  call sr3d_track_photon (photon, lat, wall, wall_hit, err, .true.)
  write (1, '(a, i0, a, 3f20.16)') '"Photon:',  n_photon, '"    ABS   1.0E-14', &
                                                            wall_hit(1)%after_reflect%vec(2:6:2)
enddo

end program
