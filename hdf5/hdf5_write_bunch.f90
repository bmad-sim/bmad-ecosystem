module hdf5_bunch_mod

use openpmd_interface
use bmad_interface

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_write_beam (file_name, bunches, lat)

type (bunch_struct), target :: bunches(:)
type (bunch_struct), pointer :: bunch
type (lat_struct) lat

integer(HID_T) f_id, g_id, z_id, r_id, b_id
integer ib, ix, h5_err

character(*) file_name
character(20) date_time, root_path, bunch_path, particle_path
character(100) this_path

! Open a new file using default properties.

call h5open_f(h5_err)  ! Init Fortran interface
call h5fcreate_f (file_name, H5F_ACC_TRUNC_F, f_id, h5_err)

! Write some header stuff

call date_and_time_stamp (date_time, .true.)
root_path = '/data/'
bunch_path = 'bunch%T/'
particle_path = 'particles/'

call pmd_write_string_attrib(f_id, 'openPMD', '2.0.0')
call pmd_write_string_attrib(f_id, 'openPMDextension', 'BeamPhysics; SpeciesType')
call pmd_write_string_attrib(f_id, 'basePath', trim(root_path) // trim(bunch_path))
call pmd_write_string_attrib(f_id, 'particlesPath', trim(particle_path))
call pmd_write_string_attrib(f_id, 'software', 'Bmad')
call pmd_write_string_attrib(f_id, 'softwareVersion', '1.0')
call pmd_write_string_attrib(f_id, 'date', date_time)
call pmd_write_string_attrib(f_id, 'latticeFile', lat%input_file_name)
if (lat%lattice /= '') then
  call pmd_write_string_attrib(f_id, 'latticeName', lat%lattice)
endif

! Loop over bunches

call h5gcreate_f(f_id, trim(root_path), r_id, h5_err)

do ib = 1, size(bunches)
  bunch => bunches(ib)
  ! Create bunch group
  ix = index(bunch_path, '%T')
  write (this_path, '(a, i0, 2a)') bunch_path(1:ix-1), ib, trim(bunch_path(ix+2:))
  call h5gcreate_f(r_id, trim(this_path), b_id, h5_err)

  call h5gcreate_f(b_id, 'position', z_id, h5_err)
  call pmd_write_real_vector_to_dataset(z_id, 'x', unit_m, bunch%particle(:)%vec(1), h5_err)
  call pmd_write_real_vector_to_dataset(z_id, 'y', unit_m, bunch%particle(:)%vec(3), h5_err)
  call pmd_write_real_to_pseudo_dataset(z_id, 'z', unit_m, 0.0_rp, size(bunch%particle), h5_err)
  call h5gclose_f(z_id, h5_err)

  call h5gclose_f(b_id, h5_err)
enddo

call h5gclose_f(r_id, h5_err)
call h5fclose_f(f_id, h5_err)

end subroutine hdf5_write_beam 

end module hdf5_bunch_mod
