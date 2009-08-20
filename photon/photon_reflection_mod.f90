module photon_reflection_mod

use precision_def
use physical_constants

integer, private, parameter :: N_energy=97, N_angles=16
real(rp), private, save :: ang_vec(N_angles) = (/ 0,1,2,3,4,5,6,7,10,15,20,30,45,60,75,90 /)
real(rp), private, save :: energy_vec(N_energy), reflect_vec(N_angles,N_energy)
real(rp), private, save :: y2(N_angles,N_energy)
logical, save, private :: init_needed = .true.

contains

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

! initialization routine for photrefp
! read in tables for interpolation

Subroutine photon_reflection_init

  use nr

  implicit none

  character(100) datafile
  character(100) line

  integer i, j

  !

  do i = 1, N_angles
    write (datafile, '(a, i0, a)') 'PhotoReflTables/xray-', nint(ang_vec(i)), '.dat'
    open(10, file=datafile, status='old')
    read(10, *) line
    read(10, *) line
    do j = 1, N_energy
      read(10, *) energy_vec(j), reflect_vec(i, j)
    enddo
  enddo

  !  set up for spine interpolation

  call splie2(ang_vec, energy_vec, reflect_vec, y2)

end subroutine

!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

! evaluate using spline interpolation
! angle is in degrees,  energy is in eV,  
! ref is the reflectivity, a pure number between 0 and 1.

subroutine photon_reflectivity (angle, energy,  ref)

  use nr

  implicit none

  real(rp) angle, angle_deg, energy, energym, ref

  !

  if (init_needed) then
    call photon_reflection_init
    init_needed = .false.
  endif

  energym = max(30.0_rp, energy)
  angle_deg = angle * 180 / pi
  ref = max(0.0_rp, splin2 (ang_vec, energy_vec, reflect_vec, y2, angle_deg, energym))

end subroutine

end module
