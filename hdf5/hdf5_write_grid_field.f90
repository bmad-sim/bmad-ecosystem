!+
! Subroutine hdf5_write_grid_field (file_name, ele, g_field, err_flag)
!
! Routine to write a binary grid_field structure.
! Note: The file name should have a ".h5" suffix.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!   g_field       -- grid_field_struct: Cylindrical map.
!
! Ouput:
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_grid_field (file_name, ele, g_field, err_flag)

use hdf5_interface
use bmad_interface

implicit none

type (grid_field_struct), target :: g_field
type (ele_struct) ele

integer(hid_t) f_id
integer i, j, k, n, h5_err
logical err_flag, err

character(*) file_name
character(*), parameter :: r_name = 'dhf5_write_grid_field'
character(200) f_name

!

err_flag = .true.
call hdf5_open_file (file_name, 'WRITE', f_id, err);  if (err) return

call hdf5_write_attribute_string(f_id, 'fileType',               'Bmad:grid_field', err)
call hdf5_write_attribute_string(f_id, 'file_name',              file_name, err)
call hdf5_write_attribute_string(f_id, 'master_parameter',       attribute_name(ele, g_field%master_parameter), err)
call hdf5_write_attribute_string(f_id, 'geometry',               grid_field_geometry_name(g_field%geometry), err)
call hdf5_write_attribute_string(f_id, 'field_type',             em_field_type_name(g_field%field_type), err)
call hdf5_write_attribute_string(f_id, 'ele_anchor_pt',          anchor_pt_name(g_field%ele_anchor_pt), err)
call hdf5_write_attribute_real(f_id,   'field_scale',            g_field%field_scale, err)
call hdf5_write_attribute_real(f_id,   'phi0_fieldmap',          g_field%phi0_fieldmap, err)
call hdf5_write_attribute_real(f_id,   'r0',                     g_field%r0, err)
call hdf5_write_attribute_real(f_id,   'dr',                     g_field%dr, err)
call hdf5_write_attribute_int(f_id,    'harmonic',               g_field%harmonic, err)
call hdf5_write_attribute_int(f_id,    'interpolation_order',    g_field%interpolation_order, err)
call hdf5_write_attribute_int(f_id,    'lbound',                 lbound(g_field%ptr%pt), err)
call hdf5_write_attribute_int(f_id,    'ubound',                 ubound(g_field%ptr%pt), err)


call h5fclose_f(f_id, h5_err)
err_flag = .false.

end subroutine hdf5_write_grid_field

