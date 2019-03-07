module hdf5_bunch_mod

use hdf5_interface
use bmad_interface

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_write_beam (file_name, bunches, lat)

type (bunch_struct), target :: bunches(:)
type (bunch_struct), pointer :: bunch
type (lat_struct) lat

integer(HID_T) fid, gid, zid
integer ib, ix, h5_err

character(*) file_name
character(20) date_time, base_path, particle_path
character(100) bunch_path

! Open a new file using default properties.

call h5open_f(h5_err)  ! Init Fortran interface
call h5fcreate_f (file_name, H5F_ACC_TRUNC_F, fid, h5_err)

! Write some header stuff

call date_and_time_stamp (date_time, .true.)
base_path = '/data/bunch%T/'
particle_path = 'particles/'

call hdf5_write_string_attrib(fid, 'openPMD', '2.0.0')
call hdf5_write_string_attrib(fid, 'openPMDextension', 'BeamPhysics; SpeciesType')
call hdf5_write_string_attrib(fid, 'basePath', base_path)
call hdf5_write_string_attrib(fid, 'particlesPath', particle_path)
call hdf5_write_string_attrib(fid, 'software', 'Bmad')
call hdf5_write_string_attrib(fid, 'softwareVersion', '1.0')
call hdf5_write_string_attrib(fid, 'date', date_time)
call hdf5_write_string_attrib(fid, 'latticeFile', lat%input_file_name)
if (lat%lattice /= '') then
  call hdf5_write_string_attrib(fid, 'latticeName', lat%lattice)
endif

! Loop over bunches

do ib = 1, size(bunches)
  bunch => bunches(ib)
  ! Create bunch group
  ix = index(base_path, '%T')
  write (bunch_path, '(a, i0, 2a)') base_path(1:ix-1), ib, trim(base_path(ix+2:)), trim(particle_path)
  call h5gcreate_f(fid, bunch_path, gid, h5_err)

  call h5gcreate(gid, 'position', zid, h5_err)
  call hdf5_write_real_vector_to_group(zid, 'x', bunch%particle(:)%vec(1))
  call hdf5_write_real_vector_to_group(zid, 'y', bunch%particle(:)%vec(3))
  call hdf5_write_real_to_group(zid, 'z', 0.0_rp, size(bunch%particle))
  call h5gclose_f(zid, h5_err)
enddo

end subroutine hdf5_write_beam 

end module hdf5_bunch_mod
