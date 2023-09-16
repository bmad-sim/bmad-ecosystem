program synrad3d_test

use synrad3d_track_mod
use synrad3d_parse_wall
use photon_reflection_mod
use synrad_mod
use random_mod

implicit none

type (lat_struct), target :: lat
type (sr3d_photon_track_struct) :: photon
type (coord_struct) p
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (ele_struct), pointer :: ele
type (wall3d_section_struct), pointer :: section(:)
type (walls_struct), target :: walls
type (wall_struct), pointer :: wall
type (photon_reflect_surface_struct) surface
type (diffuse_param_struct) d_param

real(rp) vel, graze_angle_in, energy, theta_out, phi_out
integer i, ios, num_ignored, n_photon

logical is_inside, err, absorbed

character(100) old_wall_file, wall_file

namelist / in / p, wall_file

!

open (1, file = 'output.now')

! Diffuse scattering test

call photon_reflection_std_surface_init (surface)
graze_angle_in = 1
energy = 1
call ran_seed_put (1234)
call photon_diffuse_scattering (graze_angle_in, energy, surface, theta_out, phi_out, d_param)

write (1, '(a, 3f14.9)') '"Diffuse" ABS 1E-8', theta_out, phi_out

! Conversion to synrad wall test

call bmad_parser('lat0.bmad', lat)
call synrad_read_vac_wall_geometry ('synrad3d::lat0.wall', 0.1_rp, lat%branch(0), walls, err)

wall => walls%negative_x_wall
do i = 0, wall%n_pt_max
  write (1, '(a, i0, a, 2f16.10)') '"Synrad-Wall', i, '"  ABS 1E-12', wall%pt(i)%s, wall%pt(i)%x 
enddo

! Init lattice

call bmad_parser('lat.bmad', lat)

sr3d_params%specular_reflection_only = .true.
sr3d_params%allow_absorption = .false.
num_ignored = 0
allocate (wall_hit(0:1))

! Specular reflection test...

old_wall_file = 'xxx' 
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
    call sr3d_read_wall_file (wall_file, lat)
    old_wall_file = wall_file
    cycle
  endif

  vel = sqrt(p%vec(2)**2 + p%vec(4)**2 + p%vec(6)**2)
  if (abs(vel - 1) > 0.1) then
    print *, 'ERROR: PHOTON VELOCITY NOT PROPERLY NORMALIZED TO 1 FOR PHOTON:', n_photon
    stop
  endif
  p%vec(2:6:2) = p%vec(2:6:2) / vel
  p%direction = sign_of(p%vec(6))

  p%p0c = 1000             ! Arbitrary
  p%s = p%vec(5)
  p%ix_ele = element_at_s(lat, p%s, .true.)
  p%location = inside$
  ele => lat%ele(p%ix_ele)
  p%vec(5) = p%s - (ele%s - ele%value(l$))
  photon%start%orb = p
  photon%start%ix_branch = 0
  photon%n_wall_hit = 0

  n_photon = n_photon + 1
  photon%ix_photon = n_photon
  photon%ix_photon_generated = n_photon

  call sr3d_check_if_photon_init_coords_outside_wall (photon%start, lat, is_inside, num_ignored)

  call sr3d_track_photon (photon, lat, wall_hit, err, .true.)
  write (1, '(a, i0, a, 3f20.16)') '"Photon:',  n_photon, '"    ABS   1.0E-14', wall_hit(1)%after_reflect%vec(2:6:2)
enddo

! Multi-section test

call sr3d_read_wall_file ('multi_sec.wall', lat)
section => lat%branch(0)%wall3d(1)%section
write (1, '(3a)')       '"Multi-sec22-name"  STR  "', trim(section(22)%name), '"'
write (1, '(3a)')       '"Multi-sec23-name"  STR  "', trim(section(23)%name), '"'
write (1, '(a, f12.4)') '"Multi-sec22-s"   ABS 0 ', section(22)%s
write (1, '(a, f12.4)') '"Multi-sec23-s"   ABS 0 ', section(23)%s
write (1, '(a, i3)')    '"Multi-sec-size"  ABS 0  ', size(section)

! Branch test

call bmad_parser ('branch.bmad', lat)
call sr3d_read_wall_file ('branch.wall', lat)

ele => lat%ele(1)
p%s = 1.2_rp - 1d-12
p%vec = 0
p%vec(6) = 1
p%ix_ele = 1
p%location = inside$
p%vec(5) = p%s - (ele%s - ele%value(l$))
photon%start%orb = p
photon%start%ix_branch = 0
photon%n_wall_hit = 0

call sr3d_track_photon (photon, lat, wall_hit, err, .true.)
write (1, '(a, 3f20.16)') '"Branch-Photon"    ABS   1.0E-14', wall_hit(1)%before_reflect%vec(2:6:2)
write (1, '(a, i0)')      '"Ix_Branch"        ABS   0      ', photon%now%ix_branch

!

close (1)

end program
