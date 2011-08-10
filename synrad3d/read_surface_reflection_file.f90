!+ 
! Subroutine read_surface_reflection_file (file_name)
!
! Routine to read the reflection parameters from a file for the photon_reflectivity routine. 
! 
! Input:
!   file_name -- Character(*): Name of the file
!-

subroutine read_surface_reflection_file (file_name)

use photon_reflection_mod
use sim_utils

implicit none

type (photon_reflect_table_struct), pointer :: prt

character(*) file_name

integer i, j, iu, n_table, it, n_angles, n_energy

real(rp) angles(100), energy_min, energy_max, energy_delta, row(200)

namelist / surface_general / n_table
namelist / table_general / angles, energy_min, energy_max, energy_delta
namelist / table_row / row

! Init standard table

call photon_reflection_init()

! Open file

iu = lunget()
open (iu, file = file_name, status = 'old')

! Allocate the tables

read (iu, nml = surface_general)

reflect_surface => reflect_surface_array(2)
allocate (reflect_surface%table(n_table))

! Fill in each table

do it = 1, n_table
  prt => reflect_surface%table(it)

  angles = -1
  read (iu, nml = table_general)
  do i = 1, size(angles)
    if (angles(i) < 0) then
      n_angles = i - 1
      exit
    endif
  enddo

  n_energy = 1 + (energy_max - energy_min) / energy_delta
  allocate(prt%angle(n_angles), prt%energy(n_energy), prt%p_reflect(n_angles,n_energy), prt%reflect_prob(n_angles), prt%int1(n_energy))

  prt%angle = angles(1:n_angles)
  prt%energy = [(energy_min + (i-1) * energy_delta, i = 1, n_energy)]
  prt%max_energy = prt%energy(n_energy)

  do i = 1, n_energy
    read (iu, nml = table_row)
    if (nint(row(1)) /= i) then
      print *, 'ERROR READING SURFACE FILE: ROW MISMATCH:', i, nint(row(1))
      call err_exit
    endif
    prt%p_reflect(:, i) = row(2:n_angles+1)
  enddo

enddo

close (iu)

call finalize_reflectivity_tables (reflect_surface)

end subroutine read_surface_reflection_file
