!+
! Program s3drefl
!
! Output photon reflection function
!-

program s3drefl

use synrad3d_plot_mod
use synrad3d_output_mod
use bookkeeper_mod
use photon_reflection_mod

implicit none

type (sr3d_wall_struct), target :: wall
type (photon_reflect_surface_struct), pointer :: surface

character(200) surface_reflection_file
real(rp) roughness_rms_min, roughness_rms_max, &
         correlation_len_min, correlation_len_max, &
         angle_in_min, angle_in_max, &
         energy_min, energy_max

integer roughness_rms_nsteps, correlation_len_nsteps, angle_in_nsteps, energy_nsteps

character(200) param_file, output_file

integer iu_output_file

real(rp) roughness_rms, correlation_len, rms_set, correlation_set

real(rp) dr, dc, da, de

integer i,j, k, l

integer n_surface

real(rp) angle_in, energy, theta_out, phi_out

namelist / s3drefl_parameters / surface_reflection_file, &
                                roughness_rms_min, roughness_rms_max, roughness_rms_nsteps, &
                                correlation_len_min, correlation_len_max, correlation_len_nsteps, &
                                angle_in_min, angle_in_max, angle_in_nsteps, &
                                energy_min, energy_max, energy_nsteps, &
                                output_file

! Default input parameters. Specular reflection only.

 surface_reflection_file = '' ! For custom reflection parameters.
 roughness_rms_min = 0        ! Specular only
 roughness_rms_max = 0
 roughness_rms_min = 100e-9        
 roughness_rms_max = 200e-9
 roughness_rms_nsteps = 2     ! Ignored if min>max
 correlation_len_min = 0      ! Specular only
 correlation_len_max = 0
 correlation_len_min = 5000e-9
 correlation_len_max = 10000e-9
 correlation_len_nsteps = 2   ! Ignored if min>max
 angle_in_min = 0             ! Min incident angle (degrees)
 angle_in_max = 90.           ! Max incident angle (degrees)
 angle_in_nsteps = 100        ! Ignored if min>max
 energy_min = 100.
 energy_max = 9000.
 energy_nsteps = 89           ! Ignored if min>max
 output_file = 's3drefl.out'

! Read input file

param_file = 's3drefl.init'
print *, 'Input parameter file: ', trim(param_file)
open (1, file = param_file, status = 'old')
read (1, nml = s3drefl_parameters)
close (1)


! Convert incident angle to radians
angle_in_min = angle_in_min * pi / 2. / 90.
angle_in_max = angle_in_max * pi / 2. / 90.

! Open output file

if (output_file /= '') then
  iu_output_file = lunget()
  open (iu_output_file, file = output_file)
  print *, 'Creating output file: ', trim(output_file)

! Write header
  write(iu_output_file,'(4(e15.7,2x,e15.7,i6))') &
         roughness_rms_min, roughness_rms_max, roughness_rms_nsteps, &
         correlation_len_min, correlation_len_max, correlation_len_nsteps, &
         angle_in_min, angle_in_max, angle_in_nsteps, &
         energy_min, energy_max, energy_nsteps

endif

! Load different surface reflection parameters if wanted

if (surface_reflection_file /= '') then
  call read_surface_reflection_file (surface_reflection_file, wall%surface(i))
else
  n_surface =1 
  allocate (wall%surface(n_surface))
!  call photon_reflection_init()
endif

! Calculate step sizes
dr=0.
if ( roughness_rms_nsteps >0 ) then
  if ( roughness_rms_min < roughness_rms_max ) dr = ( roughness_rms_max - roughness_rms_min ) / roughness_rms_nsteps
endif

dc=0.
if ( correlation_len_nsteps > 0 ) then
  if ( correlation_len_min < correlation_len_max ) dc = ( correlation_len_max - correlation_len_min ) / correlation_len_nsteps
endif

da=0.
if ( angle_in_nsteps > 0 ) then
  if ( angle_in_min < angle_in_max ) da = ( angle_in_max - angle_in_min ) / angle_in_nsteps
endif

de=0.
if ( energy_nsteps > 0 ) then
  if ( energy_min < energy_max ) de = ( energy_max - energy_min ) / energy_nsteps
endif

do i = 0, roughness_rms_nsteps
  roughness_rms = roughness_rms_min + i * dr

  do j = 0, correlation_len_nsteps
    correlation_len = correlation_len_min + j * dc

!    call set_surface_roughness (roughness_rms, correlation_len, rms_set, correlation_set)

    wall%surface(1)%surface_roughness_rms = roughness_rms
    wall%surface(1)%roughness_correlation_len = correlation_len
    call photon_reflection_std_surface_init (wall%surface(1))

    write(*,'(a,e22.7,3x,e22.7)')'Roughness, correlation length:', roughness_rms, correlation_len

    write(iu_output_file,'(2(e15.7,2x))') roughness_rms, correlation_len

    do k = 0, angle_in_nsteps
      angle_in = angle_in_min + k * da

      do l = 0, energy_nsteps
        energy = energy_min + l * de

        surface => wall%surface(1)
        call photon_diffuse_scattering ( angle_in, energy, surface, theta_out, phi_out )

!        write (*,'(a,f6.2,2x,e10.5,a,f6.2,2x,f6.2)') &
!              ' Incident angle (deg), energy (eV): ', 180*angle_in/pi, 1e9*energy, &
!              ' Outgoing polar and azimuthal angles: ', 180*theta_out/pi, 180*phi_out/pi

!       Write output file entry
        write(iu_output_file,'(4(e15.7,2x))') angle_in, energy, theta_out, phi_out


      end do
    end do
  end do
end do


if ( iu_output_file > 0 ) close(iu_output_file)

end program
