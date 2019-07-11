!+
! Subroutine hdf5_write_grid_field (file_name, ele, g_field, err_flag)
!
! Routine to write a binary grid_field structure.
! Note: The file name should have a ".h5" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!   g_field       -- grid_field_struct: Grid field.
!
! Ouput:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_grid_field (file_name, ele, g_field, err_flag)

use hdf5_openpmd_mod
use bmad_interface

implicit none

type (grid_field_struct), target :: g_field
type (ele_struct) ele

integer(hid_t) f_id, r_id, b_id, b2_id, z_id
integer i, j, k, n, ix, ig, h5_err
logical err_flag, err

character(*) file_name
character(*), parameter :: r_name = 'dhf5_write_grid_field'
character(40) this_path
character(28) date_time, root_path, bunch_path, mesh_path
character(2), parameter :: b_name(3) = ['Bx', 'By', 'Bz'], e_name(3) = ['Ex', 'Ey', 'Ez']

!

err_flag = .true.
call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call date_and_time_stamp (date_time, .true., .true.)
root_path = '/data/'
bunch_path = '%T/'
mesh_path = 'ExternalFieldMesh/'

call hdf5_write_attribute_string(f_id, 'dataType',               'Bmad:grid_field', err)
call hdf5_write_attribute_string(f_id, 'openPMD', '2.0.0', err)
call hdf5_write_attribute_string(f_id, 'openPMDextension', 'BeamPhysics;SpeciesType', err)
call hdf5_write_attribute_string(f_id, 'basePath', trim(root_path) // trim(bunch_path), err)
call hdf5_write_attribute_string(f_id, 'meshesPath', trim(mesh_path), err)
call hdf5_write_attribute_string(f_id, 'software', 'Bmad', err)
call hdf5_write_attribute_string(f_id, 'softwareVersion', '1.0', err)
call hdf5_write_attribute_string(f_id, 'date', date_time, err)

ix = index(bunch_path, '%T')
ig = 1   ! Only one grid to write
write (this_path, '(a, i0, a)') bunch_path(1:ix-1), ig, trim(bunch_path(ix+2:))
call h5gcreate_f(r_id, trim(this_path), b_id, h5_err)
call h5gcreate_f(b_id, mesh_path, b2_id, h5_err)

call hdf5_write_attribute_string(b2_id,  'master_parameter',     attribute_name(ele, g_field%master_parameter), err)
call hdf5_write_attribute_string(b2_id,  'geometry',             grid_field_geometry_name(g_field%geometry), err)
call hdf5_write_attribute_string(b2_id,  'fieldType',            em_field_type_name(g_field%field_type), err)
call hdf5_write_attribute_string(b2_id,  'eleAnchorPt',          anchor_pt_name(g_field%ele_anchor_pt), err)
call hdf5_write_attribute_real(b2_id,    'field_Scale',          g_field%field_scale, err)
call hdf5_write_attribute_real(b2_id,    'phi0Fieldmap',         g_field%phi0_fieldmap, err)
call hdf5_write_attribute_real(b2_id,    'r0',                   g_field%r0, err)
call hdf5_write_attribute_real(b2_id,    'dr',                   g_field%dr, err)
call hdf5_write_attribute_int(b2_id,     'harmonic',             g_field%harmonic, err)
call hdf5_write_attribute_int(b2_id,     'interpolationOrder',   g_field%interpolation_order, err)
call hdf5_write_attribute_int(b2_id,     'lbound',               lbound(g_field%ptr%pt), err)
call hdf5_write_attribute_int(b2_id,     'ubound',               ubound(g_field%ptr%pt), err)
call hdf5_write_attribute_int(b2_id,     'curvedRefFrame',       int_logic(g_field%curved_ref_frame), err)

if (g_field%field_type == magnetic$ .or. g_field%field_type == mixed$) then
  do i = 1, 3
    call h5gcreate_f(b2_id, B_name(i), z_id, h5_err)
    call hdf5_write_real_tensor_to_dataset (z_id, 'real', 'real', unit_tesla, real(g_field%ptr%pt%B(i)), err)
    call hdf5_write_real_tensor_to_dataset (z_id, 'imaginary', 'imaginary', unit_tesla, real(g_field%ptr%pt%B(i)), err)
    call h5gclose_f(z_id, h5_err)
  enddo
endif

if (g_field%field_type == electric$ .or. g_field%field_type == mixed$) then
  do i = 1, 3
    call h5gcreate_f(b2_id, E_name(i), z_id, h5_err)
    call hdf5_write_real_tensor_to_dataset (z_id, 'real', 'real', unit_tesla, real(g_field%ptr%pt%E(i)), err)
    call hdf5_write_real_tensor_to_dataset (z_id, 'imaginary', 'imaginary', unit_tesla, real(g_field%ptr%pt%E(i)), err)
    call h5gclose_f(z_id, h5_err)
  enddo
endif

call h5gclose_f(b2_id, h5_err)
call h5gclose_f(b_id, h5_err)

call h5fclose_f(f_id, h5_err)
err_flag = .false.

end subroutine hdf5_write_grid_field

