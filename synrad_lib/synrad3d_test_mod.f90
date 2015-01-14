module synrad3d_test_mod

use synrad3d_track_mod
use synrad3d_output_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_diffuse_reflection_test (param_file)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- Character(*): Input parameter file.
!-

subroutine sr3d_diffuse_reflection_test (param_file)

implicit none

type (photon_reflect_surface_struct) surface

real(rp) graze_angle_in, energy, angle_radians
real(rp) p_reflect_rough, p_reflect_smooth, theta_out, phi_out
real(rp) surface_roughness_rms, roughness_correlation_len

integer n_photons
integer i, ix, ios, random_seed

character(*) param_file
character(200) output_file, surface_reflection_file

namelist / diffuse_reflection_test / graze_angle_in, energy, n_photons, surface_roughness_rms, &
            roughness_correlation_len, surface_reflection_file, output_file, random_seed

!

open (1, file = param_file)
output_file = 'test_diffuse_reflection.dat'

random_seed = 0
read (1, nml = diffuse_reflection_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING DIFFUSE_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND DIFFUSE_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
close (1)

!

call ran_seed_put (random_seed)

if (surface_reflection_file == '') then
  call photon_reflection_std_surface_init (surface)
else
  call read_surface_reflection_file (surface_reflection_file, surface)
endif

if (surface_roughness_rms > 0) surface%surface_roughness_rms = surface_roughness_rms
if (roughness_correlation_len > 0) surface%roughness_correlation_len = roughness_correlation_len

angle_radians = graze_angle_in * pi / 180
call photon_reflectivity (angle_radians, energy, surface, p_reflect_rough, p_reflect_smooth)

!

open (2, file = output_file)

write (2, *) 'Grazing angle in (deg):    ', graze_angle_in
write (2, *) 'Energy (eV):               ', energy
write (2, *) 'surface_roughness_rms:     ', surface_roughness_rms
write (2, *) 'roughness_correlation_len: ', roughness_correlation_len
write (2, *) 'Rough surface reflection probability: ', p_reflect_rough
write (2, *) 'Smooth surface reflection probability:', p_reflect_smooth
write (2, *) 'surface_reflection_file: "', trim(surface_reflection_file), '"'
write (2, *) 'random_seed:                 "', random_seed

write (2, *) '          #          theta_out                     phi_out'
do i = 1, n_photons
  call photon_diffuse_scattering (angle_radians, energy, surface, theta_out, phi_out)
  write (2, *) i, theta_out, phi_out
enddo

close (2)
print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_diffuse_reflection_test

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_specular_reflection_test (param_file)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- Character(*): Input parameter file.
!-

subroutine sr3d_specular_reflection_test (param_file)

implicit none

type (lat_struct), target :: lat
type (sr3d_wall_struct), target :: wall
type (sr3d_photon_track_struct) :: photon
type (coord_struct) p
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (random_state_struct) ran_state
type (branch_struct), pointer :: branch

real(rp) vel
integer i, ios, num_ignored, random_seed, n_photon

logical is_inside, err, absorbed

character(*) param_file
character(100) photon_start_input_file, output_file, lattice_file, wall_file

namelist / specular_reflection_test / photon_start_input_file, output_file, lattice_file, wall_file
namelist / start / p, ran_state, random_seed

! Read parameters

open (1, file = param_file)
output_file = ''
read (1, nml = specular_reflection_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING SPECULAR_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND SPECULAR_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
close (1)

if (output_file == '') output_file = 'test_specular_reflection.dat'

! Get lattice

if (lattice_file(1:6) == 'xsif::') then
  call xsif_parser(lattice_file(7:), lat)
else
  call bmad_parser (lattice_file, lat)
endif

branch => lat%branch(0)

! Init wall

call sr3d_read_wall_file (wall_file, branch)

! Open photon start input file and count the number of photons

print *, 'Opening photon starting position input file: ', trim(photon_start_input_file)
open (1, file = photon_start_input_file, status = 'old')
open (2, file = output_file)

allocate (wall_hit(0:10))
sr3d_params%specular_reflection_only = .true.
sr3d_params%allow_absorption = .false.
num_ignored = 0
n_photon = 0

do

  read (1, nml = start, iostat = ios)
  if (ios < 0) exit 
  if (ios > 0) then
    print *, 'Error reading photon starting position at photon index:', n_photon
    call err_exit
  endif

  vel = sqrt(p%vec(2)**2 + p%vec(4)**2 + p%vec(6)**2)
  if (abs(vel - 1) > 0.1) then
    print *, 'ERROR: PHOTON VELOCITY NOT PROPERLY NORMALIZED TO 1 FOR PHOTON:', n_photon
    stop
  endif
  p%vec(2:6:2) = p%vec(2:6:2) / vel

  p%p0c = 1000             ! Arbitrary
  p%ix_ele = element_at_s(lat, p%s, .true.)
  photon%start%orb = p
  photon%n_wall_hit = 0

  call sr3d_check_if_photon_init_coords_outside_wall (photon%start, branch, is_inside, num_ignored)

  n_photon = n_photon + 1
  photon%ix_photon_generated = n_photon
  photon%ix_photon = n_photon

  call sr3d_track_photon (photon, branch, wall_hit, err, .true.)
  call print_hit_points (2, photon, wall_hit)  

enddo

print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_specular_reflection_test 

end module

