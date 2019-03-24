module hdf5_interface

use h5lt
use hdf5
use sim_utils

implicit none

! Misc

integer(hsize_t), parameter :: h5_size_1 = 1
integer(hsize_t), parameter :: h5_size_7 = 7
integer, parameter :: H5O_TYPE_ATTRIBUTE_F = 123

! %element_type identifies the type of element (group, dataset or attribute) can be:
!   H5O_TYPE_GROUP_F
!   H5O_TYPE_DATASET_F
!   H5O_TYPE_ATTRIBUTE_F   ! Defined by bmad. Not by HDF5.
!   anything else is not useful.
!
! %data_type identifies the type of the underlying data. Not relavent for groups. can be:
!   H5T_FLOAT_F
!   H5T_INTEGER_F
!   H5T_STRING_F
!   anything else is not useful.

type hdf5_info_struct
  integer :: element_type = -1         ! Type of the element. See above.
  integer :: data_type = -1            ! Type of associated data. Not used for groups. See above.
  integer(HSIZE_T) :: data_dim(3) = 0  ! Dimensions. Not used for groups. EG: Scaler data has [1, 0, 0].
  integer(SIZE_T) :: data_size = -1    ! Size of datums. Not used for groups. For strings size = # of characters.
  integer :: num_attributes = -1       ! Number of associated attributes. Used for groups and datasets only.
end type

interface hdf5_read_dataset_int
  module procedure h5lt_read_dataset_int_kind_4_rank_0
  module procedure h5lt_read_dataset_int_kind_4_rank_1
  module procedure h5lt_read_dataset_int_kind_4_rank_2
  module procedure h5lt_read_dataset_int_kind_4_rank_3
end interface

interface hdf5_read_dataset_double
  module procedure h5lt_read_dataset_double_kind_8_rank_0
  module procedure h5lt_read_dataset_double_kind_8_rank_1
  module procedure h5lt_read_dataset_double_kind_8_rank_2
  module procedure h5lt_read_dataset_double_kind_8_rank_3
end interface

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
! Note: Use H5Gclose_f to close the group.

function hdf5_open_group (f_id, group_name, error, print_error) result (g_id)

integer(HID_T) f_id, g_id
integer h5_err

logical error, print_error
logical exists

character(*) group_name
character(*), parameter :: r_name = 'hdf5_open_group'

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

end function hdf5_open_group

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function hdf5_open_object(f_id, object_name, info, error, print_error) result (obj_id)

type (hdf5_info_struct) info

integer(HID_T) f_id, obj_id
integer h5_err

logical error, print_error

character(*) object_name
character(*), parameter :: r_name = 'hdf5_open_object'

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  obj_id = hdf5_open_dataset (f_id, object_name, error, print_error) 
elseif (info%element_type == H5O_TYPE_GROUP_F) then
  obj_id = hdf5_open_group(f_id, object_name, error, print_error) 
endif

end function hdf5_open_object

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_close_object(obj_id, info)

type (hdf5_info_struct) info

integer(HID_T) obj_id
integer h5_err

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  call H5Dclose_f(obj_id, h5_err)
elseif (info%element_type == H5O_TYPE_GROUP_F) then
  call H5Gclose_f(obj_id, h5_err)
endif

end subroutine hdf5_close_object

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! Note: Use H5Dclose_f to close the group.

function hdf5_open_dataset (f_id, dataset_name, error, print_error) result (ds_id)

integer(HID_T) f_id, ds_id
integer h5_err

logical error, print_error
logical exists

character(*) dataset_name
character(*), parameter :: r_name = 'hdf5_open_dataset'

!

error = .true.
call H5Lexists_f(f_id, dataset_name, exists, h5_err, H5P_DEFAULT_F)
if (.not. exists) then
  if (print_error) then
    call out_io (s_error$, r_name, 'DATASET DOES NOT EXIST: ' // quote(dataset_name))
  endif
  return
endif
 
call H5Dopen_f (f_id, dataset_name, ds_id, h5_err, H5P_DEFAULT_F)
if (h5_err == -1) return
error = .false.

end function hdf5_open_dataset

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

call H5LTget_attribute_info_f(root_id, '.', attrib_name, info%data_dim, info%data_type, info%data_size, h5_err)
info%element_type = H5O_TYPE_ATTRIBUTE_F

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

if (info%data_type == H5T_INTEGER_F) then
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

if (info%data_type == H5T_INTEGER_F .or. info%data_type == H5T_FLOAT_F) then
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

if (info%data_type /= H5T_STRING_F) then
  if (print_error) then
    call out_io (s_error$, r_name, 'ATTRIBUTE: ' // attrib_name, 'IS NOT A STRING!')
  endif
  return
endif

allocate(character(info%data_size) :: string)
call H5LTget_attribute_string_f(root_id, '.', attrib_name, string, h5_err)
if (h5_err < 0) then
  if (print_error) then
    call out_io (s_error$, r_name, 'CANNOT READ ATTRIBUTE: ' // attrib_name)
  endif
  return
endif

error = .false.

end subroutine hdf5_get_attribute_string

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

function hdf5_get_object_info (root_id, name, error, print_error) result (info)

type (hdf5_info_struct) info
type (H5O_info_t) :: infobuf 

integer(hid_t), value :: root_id
integer stat, h5_err

character(*) name

logical error, print_error

!

error = .true.

call H5Oget_info_by_name_f(root_id, name, infobuf, h5_err)
info%element_type = infobuf%type
info%num_attributes = infobuf%num_attrs

if (info%element_type == H5O_TYPE_DATASET_F) then
  call H5LTget_dataset_info_f(root_id, name, info%data_dim, info%data_type, info%data_size, h5_err)
endif

error = .false.

end function hdf5_get_object_info

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

  SUBROUTINE h5lt_read_dataset_double_kind_8_rank_0(loc_id,dset_name,buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    REAL(KIND=8),INTENT(INout), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    f_ptr = C_LOC(buf               )
    namelen = LEN(dset_name)
    errcode = h5ltread_dataset_c(loc_id,namelen,dset_name,H5T_NATIVE_DOUBLE,f_ptr)
  END SUBROUTINE h5lt_read_dataset_double_kind_8_rank_0

  SUBROUTINE h5lt_read_dataset_double_kind_8_rank_1(loc_id,dset_name,buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    REAL(KIND=8),INTENT(INout), DIMENSION(dims(1)), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    f_ptr = C_LOC(buf(1)            )
    namelen = LEN(dset_name)
    errcode = h5ltread_dataset_c(loc_id,namelen,dset_name,H5T_NATIVE_DOUBLE,f_ptr)
  END SUBROUTINE h5lt_read_dataset_double_kind_8_rank_1

  SUBROUTINE h5lt_read_dataset_double_kind_8_rank_2(loc_id,dset_name,buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    REAL(KIND=8),INTENT(INout), DIMENSION(dims(1),dims(2)), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    f_ptr = C_LOC(buf(1,1)          )
    namelen = LEN(dset_name)
    errcode = h5ltread_dataset_c(loc_id,namelen,dset_name,H5T_NATIVE_DOUBLE,f_ptr)
  END SUBROUTINE h5lt_read_dataset_double_kind_8_rank_2

  SUBROUTINE h5lt_read_dataset_double_kind_8_rank_3(loc_id,dset_name,buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    REAL(KIND=8),INTENT(INout), DIMENSION(dims(1),dims(2),dims(3)), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    f_ptr = C_LOC(buf(1,1,1)        )
    namelen = LEN(dset_name)
    errcode = h5ltread_dataset_c(loc_id,namelen,dset_name,H5T_NATIVE_DOUBLE,f_ptr)
  END SUBROUTINE h5lt_read_dataset_double_kind_8_rank_3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

  SUBROUTINE h5lt_read_dataset_int_kind_4_rank_0(loc_id,dset_name, buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    INTEGER(KIND=4),INTENT(INout), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    INTEGER(hid_t) :: type_id
    f_ptr = C_LOC(buf               )
    namelen = LEN(dset_name)
    type_id = h5kind_to_type(KIND(buf               ), H5_INTEGER_KIND)
    errcode = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  END SUBROUTINE h5lt_read_dataset_int_kind_4_rank_0

  SUBROUTINE h5lt_read_dataset_int_kind_4_rank_1(loc_id,dset_name, buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    INTEGER(KIND=4),INTENT(INout), DIMENSION(dims(1)), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    INTEGER(hid_t) :: type_id
    f_ptr = C_LOC(buf(1)            )
    namelen = LEN(dset_name)
    type_id = h5kind_to_type(KIND(buf(1)            ), H5_INTEGER_KIND)
    errcode = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  END SUBROUTINE h5lt_read_dataset_int_kind_4_rank_1

  SUBROUTINE h5lt_read_dataset_int_kind_4_rank_2(loc_id,dset_name, buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    INTEGER(KIND=4),INTENT(INout), DIMENSION(dims(1),dims(2)), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    INTEGER(hid_t) :: type_id
    f_ptr = C_LOC(buf(1,1)          )
    namelen = LEN(dset_name)
    type_id = h5kind_to_type(KIND(buf(1,1)          ), H5_INTEGER_KIND)
    errcode = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  END SUBROUTINE h5lt_read_dataset_int_kind_4_rank_2

  SUBROUTINE h5lt_read_dataset_int_kind_4_rank_3(loc_id,dset_name, buf,dims,errcode)
    IMPLICIT NONE
    INTEGER(hid_t)  , INTENT(IN) :: loc_id
    CHARACTER(LEN=*), INTENT(IN) :: dset_name
    INTEGER(hsize_t), DIMENSION(*), INTENT(in) :: dims
    INTEGER(KIND=4),INTENT(INout), DIMENSION(dims(1),dims(2),dims(3)), TARGET :: buf
    INTEGER :: errcode 
    TYPE(C_PTR) :: f_ptr
    INTEGER(size_t) :: namelen
    INTEGER(hid_t) :: type_id
    f_ptr = C_LOC(buf(1,1,1)        )
    namelen = LEN(dset_name)
    type_id = h5kind_to_type(KIND(buf(1,1,1)        ), H5_INTEGER_KIND)
    errcode = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  END SUBROUTINE h5lt_read_dataset_int_kind_4_rank_3

end module
