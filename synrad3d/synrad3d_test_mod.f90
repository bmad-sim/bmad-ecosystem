module synrad3d_test_mod

use synrad3d_utils
use synrad3d_output_mod
use photon_reflection_mod

contains

!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------
!+
! Subroutine sr3d_reflection_test (param_file)
! 
! Routine to proform the reflection test.
!
! Input:
!   param_file    -- Character(*): Input parameter file.
!-

subroutine sr3d_reflection_test (param_file)

implicit none

real(rp) graze_angle_in, energy, angle_radians
real(rp) p_reflect_rough, p_reflect_smooth, theta_out, phi_out
real(rp) surface_roughness_rms, roughness_correlation_len

integer n_photons
integer i, ix

character(*) param_file
character(200) output_file, reflection_probability_file

namelist / reflection_test / graze_angle_in, energy, n_photons, surface_roughness_rms, &
            roughness_correlation_len, reflection_probability_file, output_file

!

open (1, file = param_file)
read (1, nml = reflection_test)
close (1)

!

if (reflection_probability_file /= '') call read_surface_reflection_file (reflection_probability_file, ix)
call set_surface_roughness (surface_roughness_rms, roughness_correlation_len)

angle_radians = graze_angle_in * pi / 180
call photon_reflectivity (angle_radians, energy, p_reflect_rough, p_reflect_smooth)

!

open (2, file = 'photons.dat')

write (2, *) 'Grazing angle in (deg):    ', graze_angle_in
write (2, *) 'Energy (eV):               ', energy
write (2, *) 'surface_roughness_rms:     ', surface_roughness_rms
write (2, *) 'roughness_correlation_len: ', roughness_correlation_len
write (2, *) 'Rough surface reflection probability: ', p_reflect_rough
write (2, *) 'Smooth surface reflection probability:', p_reflect_smooth
write (2, *) 'reflection_probability_file: "', trim(reflection_probability_file), '"'

do i = 1, n_photons
  call photon_diffuse_scattering (angle_radians, energy, theta_out, phi_out)
  write (2, *) 'theta_out      phi_out'
  write (2, *) theta_out, phi_out
enddo

close (2)

end subroutine sr3d_reflection_test

end module

