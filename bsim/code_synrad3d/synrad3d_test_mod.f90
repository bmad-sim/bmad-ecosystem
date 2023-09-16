module synrad3d_test_mod

use synrad3d_track_mod
use synrad3d_output_mod
use synrad3d_parse_wall
use photon_reflection_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_diffuse_probability_test (param_file)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- character(*): Input parameter file.
!-

subroutine sr3d_diffuse_probability_test (param_file)

implicit none

type (photon_reflect_surface_struct) surface
type (diffuse_param_struct) d_param

real(rp) surface_roughness_rms, graze_angle_in, energy, row_min, row_max, row_values(400)
real(rp) roughness_correlation_len, graze_angle_out_min, graze_angle_out_max, graze_angles_out(100)
real(rp) row_value, theta_out, phi_out, prob, graze_angle_out, p_reflect, rel_p_specular

integer i, j, row_n_pts, graze_angle_out_n_pts, ios, n_angle

logical ok, row_log_scale

character(*) param_file
character(20) row_type, prob_normalization
character(200) output_file, surface_reflection_file
character(1000) :: line
character(*), parameter :: r_name = 'sr3d_diffuse_probability_test'

namelist / diffuse_probability_test / &
    graze_angle_in, energy, surface_roughness_rms, row_min, row_max, row_n_pts, row_log_scale, &
    roughness_correlation_len, graze_angle_out_min, graze_angle_out_max, graze_angle_out_n_pts, graze_angles_out, &
    output_file, surface_reflection_file, row_type, prob_normalization

!

row_values = real_garbage$
row_n_pts = -1
row_log_scale = .false.
graze_angles_out = real_garbage$
graze_angle_out_n_pts = -1
prob_normalization = ''
surface_roughness_rms = -1
roughness_correlation_len = -1
output_file = ''
surface_reflection_file = ''
output_file = 'diffuse_probability.dat'

open (1, file = param_file, status = 'old')
read (1, nml = diffuse_probability_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING DIFFUSE_PROBABILITY_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND DIFFUSE_PROBABILITY_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
close (1)

!

if (surface_reflection_file == '') then
  call photon_reflection_std_surface_init (surface)
else
  call read_surface_reflection_file (surface_reflection_file, surface)
endif

if (roughness_correlation_len > 0) surface%roughness_correlation_len = roughness_correlation_len
if (surface_roughness_rms > 0)     surface%surface_roughness_rms = surface_roughness_rms

if (graze_angle_out_n_pts > 0) then
  n_angle = graze_angle_out_n_pts
else
  do j = 1, size(graze_angles_out)
    if (graze_angles_out(j) /= real_garbage$) cycle
    n_angle = j
    exit
  enddo
endif

!

open (2, file = output_file, recl = 1000)

write (2, '(a)')    '#                         |               Graze_Angle_Out'
line =              '# Index ' // row_type // ' |'

do j = 1, n_angle
  if (graze_angle_out_n_pts > 0) then
    graze_angle_out = graze_angle_out_min + (j - 1) * (graze_angle_out_max - graze_angle_out_min) / max(graze_angle_out_n_pts - 1, 1)
  else
    graze_angle_out = graze_angles_out(i)
  endif
  write (line(20+j*10:), '(f10.6)') graze_angle_out
enddo

write (2, '(a)') trim(line)

do i = 1, row_n_pts
  if (row_log_scale) then
    row_value = exp(log(row_min) + (i - 1) * (log(row_max) - log(row_min)) / max(row_n_pts - 1, 1))
  else
    row_value = row_min + (i - 1) * (row_max - row_min) / max(row_n_pts - 1, 1)
  endif

  select case (row_type)
  case ('energy');            energy = row_value
  case ('roughness');         surface%surface_roughness_rms = row_value
  case ('correlation');       surface%roughness_correlation_len = row_value
  case ('graze_angle_in');    graze_angle_in = row_value
  case ('azimuth_angle_out'); phi_out = row_value
  case default
    call out_io (s_fatal$, r_name, 'BAD ROW_TYPE PARAMETER: ' // row_type)
    stop
  end select

  call photon_reflectivity (graze_angle_in, energy, surface, p_reflect, rel_p_specular)
  call photon_diffuse_scattering (graze_angle_in, energy, surface, theta_out, phi_out, d_param)

  write (line, '(i7, es14.6)') i, row_value

  do j = 1, n_angle
    if (graze_angle_out_n_pts > 0) then
      graze_angle_out = graze_angle_out_min + (j - 1) * (graze_angle_out_max - graze_angle_out_min) / max(graze_angle_out_n_pts - 1, 1)
    else
      graze_angle_out = graze_angles_out(i)
    endif
    call spline_evaluate(d_param%prob_spline(0:d_param%n_pt_spline), sin(graze_angle_out), ok, prob)
    prob = prob / d_param%chx_norm

    if (row_type == 'azimuth_angle_out') then
      prob = prob * ptwo(surface%surface_roughness_rms, surface%roughness_correlation_len, phi_out, d_param) / d_param%c_norm
    endif

    select case (prob_normalization)
    case ('1');               ! Nothing to do
    case ('reflect');         prob = prob * (1 - rel_p_specular)
    case ('all');             prob = prob * (1 - rel_p_specular) * p_reflect
    case default
    call out_io (s_fatal$, r_name, 'BAD PROB_NORMALIZATION PARAMETER: ' // prob_normalization)
    stop
  end select

    write (line(20+j*10:), '(f10.6)') prob
  enddo

  write (2, '(a)') trim(line)
enddo
  

close (2)
print *, 'Output file: ', trim(output_file)

end subroutine sr3d_diffuse_probability_test

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_monte_carlo_reflection_test (param_file)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- character(*): Input parameter file.
!-

subroutine sr3d_monte_carlo_reflection_test (param_file)

implicit none

type (photon_reflect_surface_struct) surface

real(rp) graze_angle_in, energy
real(rp) p_reflect, rel_p_specular, theta_out, phi_out
real(rp) surface_roughness_rms, roughness_correlation_len

integer n_photons
integer i, ix, ios, random_seed

logical include_specular_reflections

character(*) param_file
character(200) output_file, surface_reflection_file

namelist / reflection_test / graze_angle_in, energy, n_photons, surface_roughness_rms, &
            roughness_correlation_len, surface_reflection_file, output_file, random_seed, include_specular_reflections

! Read parameters

random_seed = 0
include_specular_reflections = .true.
output_file = 'test_monte_carlo_reflection.dat'

open (1, file = param_file, status = 'old')
read (1, nml = reflection_test, iostat = ios)
if (ios > 0) then
  print *, 'ERROR READING MONTE_CARLO_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
  stop
endif
if (ios < 0) then
  print *, 'CANNOT FIND MONTE_CARLO_REFLECTION_TEST NAMELIST IN FILE: ' // trim(param_file)
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
  if (include_specular_reflections) then
    call photon_reflection (graze_angle_in, energy, surface, theta_out, phi_out)
  else
    call photon_diffuse_scattering (graze_angle_in, energy, surface, theta_out, phi_out)
  endif
  write (2, *) i, theta_out, phi_out
enddo

close (2)
print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_monte_carlo_reflection_test

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

call bmad_parser (lattice_file, lat)

! Init wall

call sr3d_read_wall_file (wall_file, lat)

! Open photon start input file and count the number of photons

print *, 'Opening photon starting position input file: ', trim(photon_start_input_file)
open (1, file = photon_start_input_file, status = 'old')

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

  call sr3d_track_photon (photon, lat, wall_hit, err, one_reflection_only = .true.)
  call sr3d_write_hit_points (output_file, photon, wall_hit, lat)

enddo

print *, 'Output file: ' // trim(output_file)

end subroutine sr3d_specular_reflection_test 

end module

