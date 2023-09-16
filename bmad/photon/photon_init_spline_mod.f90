!+
! Module photon_init_spline_mod
!
! Module for initializing phtons given the appropriate splines fits
! to the photon probability distributions.
!
! !!! NOTE: Currently this module is not used and is not being developed. !!!
!-

module photon_init_spline_mod

use bmad_routine_interface
use spline_mod

type photon_init_x_angle_spline_struct
  type (spline_struct), allocatable :: prob(:), pl(:), pc(:), pl45(:)
end type

type photon_init_y_angle_spline_struct
  type (spline_struct), allocatable :: prob(:), pl(:), pc(:), pl45(:)
  type (photon_init_x_angle_spline_struct), allocatable :: x_angle(:)
end type

type photon_init_splines_struct
  character(16) source_type                ! 'bend', 'wiggler', 'undulator'
  integer spline_space_dimensions          ! Dimensions: [energy, y_angle, x_angle, x, y]
  type (spline_struct), allocatable :: energy_prob(:)
  type (photon_init_y_angle_spline_struct), allocatable :: y_angle(:)
end type

contains

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
!+
! Subroutine photon_read_spline (spline_dir, splines)
!
! Routine to initialize a photon using a set of spline fits.
!
! Input:
!   spline_dir -- character(*): Root directory for the spline fits.
!
! Output:
!   splines   -- photon_init_splines_struct: Spline structure
!-

subroutine photon_read_spline (spline_dir, splines)

implicit none

type (photon_init_splines_struct) splines
type (spline_struct), allocatable, target :: prob_spline(:), pl_spline(:), pc_spline(:), pl45_spline(:)

character(*) spline_dir
character(len(spline_dir)+1) s_dir
character(20) source_type, basename
character(200) angle_spline_file

real(rp) dE_spline_max, dP_spline_max

integer i, j, n, ix, num_rows_energy, num_rows_y_angle, num_rows_x_angle, iu
integer spline_space_dimensions

namelist / master_params / source_type, dE_spline_max, dP_spline_max, num_rows_energy, &
            num_rows_y_angle, num_rows_x_angle, spline_space_dimensions
namelist / spline / prob_spline, pl_spline, pc_spline, pl45_spline

! Add '/' suffix if needed

s_dir = spline_dir
n = len_trim(spline_dir)
if (spline_dir(n:n) /= '/') s_dir = trim(s_dir) // '/'

! Read general parameters

iu = lunget()
open (iu, file = trim(s_dir) // 'spline.params', status = 'old')
read (iu, nml = master_params)
close(iu)

splines%source_type = source_type
splines%spline_space_dimensions = spline_space_dimensions

! Read energy spline

iu = lunget()
open (iu, file = trim(s_dir) // 'spline/energy.spline', status = 'old')
allocate (prob_spline(num_rows_energy))
read(iu, nml = spline)
close (iu)

call move_alloc (prob_spline, splines%energy_prob)

!---------------------------------
! Read bend vertical angle splines

if (spline_space_dimensions == 2) then

  allocate (splines%y_angle(num_rows_energy))

  do i = 1, num_rows_energy
    write (angle_spline_file, '(2a, i0, a)') trim(s_dir), 'spline/e', i, '_y.spline'
    open (iu, file = angle_spline_file, status = 'old')
    allocate (prob_spline(num_rows_y_angle))
    if (spline_space_dimensions == 2) allocate (pl_spline(num_rows_y_angle),  &
                                 pc_spline(num_rows_y_angle), pl45_spline(num_rows_y_angle))
    read(iu, nml = spline)
    close (iu)

    call move_alloc (prob_spline, splines%y_angle(i)%prob)
    if (spline_space_dimensions == 2) then
      call move_alloc (pl_spline, splines%y_angle(i)%pl)
      call move_alloc (pc_spline, splines%y_angle(i)%pc)
      call move_alloc (pl45_spline, splines%y_angle(i)%pl45)
    endif
  enddo
  return

endif

!---------------------------------
! 3D


end subroutine photon_read_spline


end module
