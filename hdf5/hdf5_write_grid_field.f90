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

integer(hid_t) f_id, r_id, b_id, z_id
integer i, j, k, n, ix, ig, h5_err
logical err_flag, err

character(*) file_name
character(*), parameter :: r_name = 'dhf5_write_grid_field'
character(40) this_path
character(28) date_time, root_path
character(2), parameter :: b_name(3) = ['Bx', 'By', 'Bz'], e_name(3) = ['Ex', 'Ey', 'Ez']

!

err_flag = .true.
call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call date_and_time_stamp (date_time, .true., .true.)
root_path = '/ExternalFieldMesh/'

call hdf5_write_attribute_string(f_id, 'dataType',          'Bmad:grid_field', err)
call hdf5_write_attribute_string(f_id, 'openPMD',           '2.0.0', err)
call hdf5_write_attribute_string(f_id, 'openPMDextension',  'BeamPhysics;SpeciesType', err)
call hdf5_write_attribute_string(f_id, 'externalFieldPath',  trim(root_path), err)
call hdf5_write_attribute_string(f_id, 'software',          'Bmad', err)
call hdf5_write_attribute_string(f_id, 'softwareVersion',   '1.0', err)
call hdf5_write_attribute_string(f_id, 'date',              date_time, err)

call h5gcreate_f(r_id, trim(root_path), b_id, h5_err)

call hdf5_write_attribute_string(b_id,  'masterParameter',      attribute_name(ele, g_field%master_parameter), err)
call hdf5_write_attribute_string(b_id,  'geometry',             grid_field_geometry_name(g_field%geometry), err)
call hdf5_write_attribute_string(b_id,  'fieldType',            em_field_type_name(g_field%field_type), err)
call hdf5_write_attribute_string(b_id,  'eleAnchorPt',          anchor_pt_name(g_field%ele_anchor_pt), err)
call hdf5_write_attribute_real(b_id,    'field_Scale',          g_field%field_scale, err)
call hdf5_write_attribute_real(b_id,    'phi0Fieldmap',         g_field%phi0_fieldmap, err)
call hdf5_write_attribute_real(b_id,    'r0',                   g_field%r0, err)
call hdf5_write_attribute_real(b_id,    'dr',                   g_field%dr, err)
call hdf5_write_attribute_int(b_id,     'harmonic',             g_field%harmonic, err)
call hdf5_write_attribute_int(b_id,     'interpolationOrder',   g_field%interpolation_order, err)
call hdf5_write_attribute_int(b_id,     'lbound',               lbound(g_field%ptr%pt), err)
call hdf5_write_attribute_int(b_id,     'ubound',               ubound(g_field%ptr%pt), err)
call hdf5_write_attribute_int(b_id,     'curvedRefFrame',       int_logic(g_field%curved_ref_frame), err)

if (g_field%field_type == magnetic$ .or. g_field%field_type == mixed$) then
  do i = 1, 3
    call h5gcreate_f(b_id, B_name(i), z_id, h5_err)
    call pmd_write_real_tensor_to_dataset (z_id, 'real', 'real', unit_tesla, real(g_field%ptr%pt%B(i)), err)
    call pmd_write_real_tensor_to_dataset (z_id, 'imaginary', 'imaginary', unit_tesla, real(g_field%ptr%pt%B(i)), err)
    call h5gclose_f(z_id, h5_err)
  enddo
endif

if (g_field%field_type == electric$ .or. g_field%field_type == mixed$) then
  do i = 1, 3
    call h5gcreate_f(b_id, E_name(i), z_id, h5_err)
    call pmd_write_real_tensor_to_dataset (z_id, 'real', 'real', unit_tesla, real(g_field%ptr%pt%E(i)), err)
    call pmd_write_real_tensor_to_dataset (z_id, 'imaginary', 'imaginary', unit_tesla, real(g_field%ptr%pt%E(i)), err)
    call h5gclose_f(z_id, h5_err)
  enddo
endif

call h5gclose_f(b_id, h5_err)
call h5fclose_f(f_id, h5_err)
err_flag = .false.

end subroutine hdf5_write_grid_field

