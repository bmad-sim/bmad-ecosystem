!+ 
! Subroutine write_lattice_in_foreign_format (out_type, out_file_name, lat, ref_orbit, &
!        use_matrix_model, include_apertures, dr12_drift_max, ix_branch, err)
!
! Subroutine to write a Elegant, MAD-8, MAD-X, OPAL, SAD, or JULIA, lattice file using the 
! information in a lat_struct. Optionally, only part of the lattice can be generated.
!
! To write a Bmad lattice file, use: write_bmad_lattice_file
!
! Note: When translating to MAD8: sad_mult and patch element are translated
!  to a MAD8 matrix element (which is a 2nd order map). In this case, the ref_orbit orbit is
!  used as the reference orbit for construction of the 2nd order map.
!
! If a sad_mult or patch element is translated to a matrix element, and the referece orbit
! is non-zero, the calculation must use 2nd order maps thourghout in order to avoid "feed down".
! If the PTC map order is different from 2, PTC will be temperarily switched to 2. 
!
! The MAD drift model is approximate and this can be a problem if the reference orbit is large.
! For a drift, the value of transfer matrix element R12 is equal to L/(1+pz) for small
! deviations of the ref_orbit from zero. dr12_drift_max sets the maximum deviation of R12 beyound 
! which an extra matrix element is inserted to make the MAD model better agree with Bmad.
!
! Note: sol_quad elements are replaced by a drift-matrix-drift or solenoid-quad model.
! Note: wiggler elements are replaced by a drift-matrix-drift or drift-bend model.
!
! Input:
!   out_type          -- character(*): Either 'ELEGANT', 'MAD-8', 'MAD-X', 'SAD', or 'OPAL-T', 'JULIA'.
!   out_file_name     -- character(*): Name of the mad output lattice file.
!   lat               -- lat_struct: Holds the lattice information.
!   ref_orbit(0:)     -- coord_struct, allocatable, optional: Referece orbit for sad_mult and patch elements.
!                          This argument must be present if the lattice has sad_mult or patch elements and is
!                          being translated to MAD-8 or SAD.
!   use_matrix_model  -- logical, optional: Use a drift-matrix_drift model for wigglers/undulators?
!                           [A MAD "matrix" is a 2nd order Taylor map.] This switch is ignored for SAD conversion.
!                           Default is False -> Use a bend-drift-bend model. 
!                           Note: sol_quad elements always use a drift-matrix-drift model.
!   include_apertures -- logical, optional: If True (the default), add to the output lattice a zero length
!                           collimator element next to any non-collimator element that has an aperture.
!                           Note: MADX translations for non-drift elements can handle non-collimator elements 
!                           with an aperture so in this case this argument is ignored.
!   dr12_drift_max    -- real(rp), optional: Max deviation for drifts allowed before a correction matrix element
!                           is added. Default value is 1d-5. A negative number means use default.
!   ix_branch         -- Integer, optional: Index of lattice branch to use. Default = 0.
!
! Output:
!   err               -- logical, optional: Set True if, say a file could not be opened.
!-

subroutine write_lattice_in_foreign_format (out_type, out_file_name, lat, ref_orbit, &
                    use_matrix_model, include_apertures, dr12_drift_max, ix_branch, err)

use bmad_interface, dummy => write_lattice_in_foreign_format
use opal_interface_mod, only: write_opal_lattice_file

implicit none

type (lat_struct), target :: lat
type (coord_struct), allocatable, optional :: ref_orbit(:)

real(rp), optional :: dr12_drift_max

integer, optional :: ix_branch
integer iu, ios

character(*), parameter :: r_name = "write_lattice_in_foreign_format"
character(*) out_type, out_file_name
character(300) line

logical, optional :: use_matrix_model, include_apertures, err

! 

select case (out_type)

case ('JULIA')
  call write_lattice_in_julia (out_file_name, lat)

case ('ELEGANT')
  call write_lattice_in_elegant_format (out_file_name, lat, ref_orbit, use_matrix_model, &
                                     include_apertures, dr12_drift_max, ix_branch, err)

case ('MAD-8', 'MAD-X')
  call write_lattice_in_mad_format (out_type, out_file_name, lat, ref_orbit, use_matrix_model, &
                                     include_apertures, dr12_drift_max, ix_branch, err)

case ('SAD')
  call write_lattice_in_sad_format (out_file_name, lat, include_apertures, ix_branch, err)

case ('OPAL-T')
  iu = lunget()
  call fullfilename (out_file_name, line)
  open (iu, file = line, iostat = ios)
  if (ios /= 0) then
    call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(out_file_name))
    return
  endif

  call write_opal_lattice_file (iu, lat, err)
  close (iu)

case default
  call out_io (s_error$, r_name, 'BAD OUT_TYPE: ' // out_type)
  return

end select

end subroutine write_lattice_in_foreign_format
