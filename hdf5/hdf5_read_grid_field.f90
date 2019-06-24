!+
! Subroutine hdf5_read_grid_field (file_name, ele, g_field, err_flag)
!
! Routine to read a binary grid_field structure.
!
! Input:
!   file_name     -- character(*): File to create.
!   ele           -- ele_struct: Element associated with the map.
!
! Ouput:
!   g_field       -- grid_field_struct, cylindrical map.
!   err_flag      -- logical: Set True if there is an error. False otherwise.
!-

subroutine hdf5_read_grid_field (file_name, ele, g_field, err_flag)

use hdf5_interface
use bmad_interface

implicit none

type (grid_field_struct), target :: g_field
type (ele_struct) ele

integer(HID_T) f_id
integer i, j, k, n0(3), n1(3), iver, h5_err
logical err_flag, err

character(*) file_name
character(40) master_name
character(*), parameter :: r_name = 'hdf5_read_grid_field'

!

err_flag = .true.
call hdf5_open_file (file_name, 'READ', f_id, err);  if (err) return

call h5fclose_f(f_id, h5_err)
err_flag = .false.

end subroutine hdf5_read_grid_field
