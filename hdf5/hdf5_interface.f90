module hdf5_interface

use h5lt
use hdf5
use sim_utils

implicit none

! Misc

integer(hsize_t), parameter :: h5_size_1 = 1
integer(hsize_t), parameter :: h5_size_7 = 7

type hdf5_info_struct
  integer :: type_class = -1
  integer(HSIZE_T) :: dimension(3) = 0   ! Dimensions. EG: Scaler attributes are [1, 0, 0].
  integer(SIZE_T) :: size
end type

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_write_string_attrib(id, name, string)
integer(HID_T) :: id
character(*) :: name, string
integer error
!
call H5LTset_attribute_string_f(id, '.', name, trim(string), error)
end subroutine hdf5_write_string_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_write_integer_attrib(id, name, i)
integer(HID_T) :: id
character(*) :: name
integer :: i
integer error
!
call H5LTset_attribute_int_f(id, '.', name, [i], h5_size_1, error)
end subroutine hdf5_write_integer_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_write_real_attrib(id, name, x)
integer(HID_T) :: id
character(*) :: name
real(rp) :: x
integer error
call H5LTset_attribute_double_f(id, '.', name, [x], h5_size_1, error)
end subroutine hdf5_write_real_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_find_group(f_id, group_name, g_id, error, print_error)

integer(HID_T) f_id, g_id
integer h5_err

logical error, print_error
logical exists

character(*) group_name
character(*), parameter :: r_name = 'hdf5_find_group'

!

error = .true.
call H5Lexists_f(f_id, group_name, exists, h5_err, H5P_DEFAULT_F)
if (.not. exists) then
  if (print_error) then
    call out_io (s_error$, r_name, 'GROUP DOES NOT EXIST: ' // quote(group_name))
  endif
  return
endif
 
call H5Gopen_f (f_id, group_name, g_id, h5_err, H5P_DEFAULT_F)
if (h5_err == -1) return
error = .false.

end subroutine hdf5_find_group

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function hdf5_num_attributes(root_id) result (num)

integer(HID_T) :: root_id
integer num, h5_err

!

call H5Aget_num_attrs_f (root_id, num, h5_err)

end function hdf5_num_attributes

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_get_attribute_by_index(root_id, attrib_indx, attrib_id, attrib_name)

integer(HID_T) root_id, attrib_id
integer(SIZE_T) nam_len
integer attrib_indx, h5_err

character(*) attrib_name

!

call H5Aopen_by_idx_f (root_id, ".", H5_INDEX_CRT_ORDER_F, H5_ITER_INC_F, int(attrib_indx-1, HSIZE_T), &
                                                                 attrib_id, h5_err, aapl_id=H5P_DEFAULT_F)
nam_len = len(attrib_name)
call H5Aget_name_f(attrib_id, nam_len, attrib_name, h5_err)
call H5Aclose_f(attrib_id, h5_err)

end subroutine hdf5_get_attribute_by_index

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function hdf5_attribute_info(root_id, attrib_name, error, print_error) result (info)

type (hdf5_info_struct) info

integer(HID_T) root_id, a_id
integer h5_err

logical error, print_error, exists

character(*) attrib_name
character(*), parameter :: r_name = 'hdf5_attribute_info'

!

error = .true.

call H5Aexists_f (root_id, attrib_name, exists, h5_err)
if (.not. exists .or. h5_err == -1) then
  if (print_error) then
    call out_io (s_error$, r_name, 'ATTRIBUTE IS NOT PRESENT: ' // attrib_name)
  endif
  return
endif

call H5LTget_attribute_info_f(root_id, '.', attrib_name, info%dimension, info%type_class, info%size, h5_err)

if (h5_err < 0) then
  if (print_error) call out_io (s_error$, r_name, 'CANNOT FILE ATTRIBUTE: ' // attrib_name)
  return
endif

error = .false.

end function hdf5_attribute_info

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function hdf5_get_attribute_int(root_id, attrib_name, error, print_error) result (attrib_value)

type (hdf5_info_struct) info

integer(HID_T) root_id, a_id
integer attrib_value
integer h5_err, type_class, a_val(10)

logical error, print_error

character(*) attrib_name
character(*), parameter :: r_name = 'hdf5_get_attribute_int'

!

attrib_value = 0

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)

if (info%type_class == H5T_INTEGER_F) then
  call H5LTget_attribute_int_f(root_id, '.', attrib_name, a_val, h5_err)
  attrib_value = a_val(1)
else
  if (print_error) call out_io (s_error$, r_name, 'ATTRIBUTE IS NOT OF INTEGER TYPE: ' // attrib_name)
  return
endif

error = .false.

end function hdf5_get_attribute_int

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function hdf5_get_attribute_real(root_id, attrib_name, error, print_error) result (attrib_value)

type (hdf5_info_struct) info

integer(HID_T) root_id, a_id
integer h5_err

real(rp) attrib_value, a_val(10)

logical error, print_error

character(*) attrib_name
character(*), parameter :: r_name = 'hdf5_get_attribute_real'

!

attrib_value = 0
error = .true.

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)

if (info%type_class == H5T_INTEGER_F .or. info%type_class == H5T_FLOAT_F) then
  call H5LTget_attribute_double_f(root_id, '.', attrib_name, a_val, h5_err)
  attrib_value = a_val(1)
else
  if (print_error) call out_io (s_error$, r_name, 'ATTRIBUTE IS NOT OF REAL TYPE: ' // attrib_name)
  return
endif

error = .false.

end function hdf5_get_attribute_real

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_get_attribute_string(root_id, attrib_name, string, error, print_error)

type (hdf5_info_struct) info

integer(HID_T) root_id, a_id
integer attrib_value
integer h5_err

logical error, print_error

character(*) attrib_name
character(:), allocatable :: string
character(*), parameter :: r_name = 'hdf5_get_attribute_string'

!

attrib_value = 0

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)

if (info%type_class /= H5T_STRING_F) then
  if (print_error) then
    call out_io (s_error$, r_name, 'ATTRIBUTE: ' // attrib_name, 'IS NOT A STRING!')
  endif
  return
endif

allocate(character(info%size) :: string)
call H5LTget_attribute_string_f(root_id, '.', attrib_name, string, h5_err)
if (h5_err < 0) then
  if (print_error) then
    call out_io (s_error$, r_name, 'CANNOT READ ATTRIBUTE: ' // attrib_name)
  endif
  return
endif

error = .false.

end subroutine hdf5_get_attribute_string

end module
