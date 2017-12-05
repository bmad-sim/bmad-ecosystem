module synrad3d_test_mod

use synrad3d_track_mod
use synrad3d_output_mod
use synrad3d_parse_wall

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_roughness_scan_test (param_file)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- character(*): Input parameter file.
!-

subroutine sr3d_roughness_scan_test (param_file)

implicit none

type (photon_reflect_surface_struct) surface
type (diffuse_param_struct) d_param

real(rp) graze_angle_in, energy, surface_roughness_min, surface_roughness_max
real(rp) roughness_correlation_len, graze_angle_out_min, graze_angle_out_max, graze_angles_out(100)
real(rp) roughness, theta_out, phi_out, prob, graze_angle

integer i, j, surface_roughness_n_pts, graze_angle_out_n_pts, random_seed, ios, n_angle

logical ok

character(*) param_file
character(200) output_file, surface_reflection_file
character(1000) :: line

namelist / roughness_scan_test / &
    graze_angle_in, energy, surface_roughness_min, surface_roughness_max, surface_roughness_n_pts, &
    roughness_correlation_len, graze_angle_out_min, graze_angle_out_max, graze_angle_out_n_pts, graze_angles_out, &
    random_seed, output_file, surface_reflection_file

!

graze_angles_out = -1
graze_angle_out_min = -1 
graze_angle_out_max = -1 
graze_angle_out_n_pts = -1 
output_file = ''
random_seed = 0
surface_reflection_file = ''
output_file = 'roughness_scan.dat'

open (1, file = param_file, status = 'old')
read (1, nml = roughness_scan_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
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

if (roughness_correlation_len > 0) surface%roughness_correlation_len = roughness_correlation_len

if (graze_angles_out(1) > 0) then
  do j = 1, size(graze_angles_out)
    if (graze_angles_out(j) > 0) cycle
    n_angle = j
    exit
  enddo
else
  n_angle = graze_angle_out_n_pts
endif

!

open (2, file = output_file, recl = 1000)

write (2, '(a)')    '#                     |               Graze_Angle_Out'
line =              '# Index     Roughness |'

do j = 1, n_angle
  if (graze_angles_out(1) > 0) then
    graze_angle = graze_angles_out(i)
  else
    graze_angle = graze_angle_out_min + (j - 1) * (graze_angle_out_max - graze_angle_out_min) / max(graze_angle_out_n_pts - 1, 1)
  endif
  write (line(16+j*10:), '(f10.6)') graze_angle
enddo

write (2, '(a)') trim(line)

do i = 1, surface_roughness_n_pts
  roughness = exp(log(surface_roughness_min) + (i - 1) * &
          (log(surface_roughness_max) - log(surface_roughness_min)) / max(surface_roughness_n_pts - 1, 1))
  surface%surface_roughness_rms = roughness
  call photon_diffuse_scattering (graze_angle_in, energy, surface, theta_out, phi_out, d_param)

  write (line, '(i7, es14.6)') i, roughness

  do j = 1, n_angle
    if (graze_angles_out(1) > 0) then
      graze_angle = graze_angles_out(i)
    else
      graze_angle = graze_angle_out_min + (j - 1) * (graze_angle_out_max - graze_angle_out_min) / max(graze_angle_out_n_pts - 1, 1)
    endif
    call spline_evaluate(d_param%prob_spline(0:d_param%n_pt_spline), graze_angle, ok, prob)
    write (line(16+j*10:), '(f10.6)') prob / d_param%chx_norm
  enddo

  write (2, '(a)') trim(line)
enddo
  

close (2)
print *, 'Output file: ', trim(output_file)

end subroutine sr3d_roughness_scan_test

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_reflection_test (param_file, who)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- character(*): Input parameter file.
!   who           -- character(*): "reflection" or "diffuse_reflection"
!-

subroutine sr3d_reflection_test (param_file, who)

implicit none

type (photon_reflect_surface_struct) surface

real(rp) graze_angle_in, energy
real(rp) p_reflect, rel_p_specular, theta_out, phi_out
real(rp) surface_roughness_rms, roughness_correlation_len

integer n_photons
integer i, ix, ios, random_seed

character(*) param_file, who
character(200) output_file, surface_reflection_file

namelist / reflection_test / graze_angle_in, energy, n_photons, surface_roughness_rms, &
            roughness_correlation_len, surface_reflection_file, output_file, random_seed

! Set defaults

select case (who)
case ('reflection');          output_file = 'test_reflection.dat'
case ('diffuse_reflection');  output_file = 'test_diffuse_reflection.dat'
case default;                 call err_exit
end select

random_seed = 0

! Read parameters

open (1, file = param_file, status = 'old')
read (1, nml = reflection_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
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

call photon_reflectivity (graze_angle_in, energy, surface, p_reflect, rel_p_specular)

!

open (2, file = output_file)

write (2, *) 'Grazing angle in (rad):    ', graze_angle_in
write (2, *) 'Energy (eV):               ', energy
write (2, *) 'surface_roughness_rms:     ', surface_roughness_rms
write (2, *) 'roughness_correlation_len: ', roughness_correlation_len
write (2, *) 'Reflection probability:    ', p_reflect
write (2, *) 'Specular Reflection / Reflection probability ratio:', rel_p_specular
write (2, *) 'surface_reflection_file: "', trim(surface_reflection_file), '"'
write (2, *) 'random_seed:                 "', random_seed

write (2, *) '          #          theta_out                     phi_out'
do i = 1, n_photons
  select case (who)
  case ('reflection')
    call photon_reflection (graze_angle_in, energy, surface, theta_out, phi_out)
  case ('diffuse_reflection')
    call photon_diffuse_scattering (graze_angle_in, energy, surface, theta_out, phi_out)
  end select
  write (2, *) i, pi/2-theta_out, phi_out
enddo

close (2)
print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_reflection_test

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

! Init wall

call sr3d_read_wall_file (wall_file, lat)

! Open photon start input file and count the number of photons

print *, 'Opening photon starting position input file: ', trim(photon_start_input_file)
open (1, file = photon_start_input_file, status = 'old')
open (2, file = output_file)

allocate (wall_hit(0:10))
sr3d_params%specular_reflection_only = .true.
sr3d_params%allow_absorption = .false.
num_ignored = 0
n_photon = 0

branch => lat%branch(0)

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

  call sr3d_check_if_photon_init_coords_outside_wall (photon%start, lat, is_inside, num_ignored)

  n_photon = n_photon + 1
  photon%ix_photon_generated = n_photon
  photon%ix_photon = n_photon

  call sr3d_track_photon (photon, lat, wall_hit, err, .true.)
  call sr3d_print_hit_points (2, photon, wall_hit, branch)

enddo

print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_specular_reflection_test 

end module

