!+
! Module hdf5_interface
!
! Interface routines for HDF5.
!
! HDF5 (Hierarchical Data Format version 5) is a set of file format design to store large amounts of data.
! See the web documentation on HDF5 for more info.
!-

module hdf5_interface

use h5lt
use hdf5

use sim_utils
use fortran_cpp_utils

use iso_fortran_env
USE iso_c_binding

implicit none

! Misc

integer, parameter :: H5O_TYPE_ATTRIBUTE_F = 123

! %element_type identifies the type of element (group, dataset or attribute) can be:
!   H5O_TYPE_GROUP_F
!   H5O_TYPE_DATASET_F
!   H5O_TYPE_ATTRIBUTE_F   ! Defined by Bmad. Not by HDF5.
!   Anything else is not useful.
!
! %data_class_type identifies the type of the underlying data. Not relavent for groups. can be:
!   H5T_FLOAT_F
!   H5T_INTEGER_F
!   H5T_STRING_F
!   H5T_COMPOUND_F      ! A compound type is used for storing complex numbers.
!   Anything else is not useful.
! For further info see the HDF5 "Datatype Interface API" help.

type hdf5_info_struct
  integer :: element_type = -1         ! Type of the element. See above.
  integer :: data_class_type = -1      ! Class type of associated data. Not used for groups. See above.
  integer(hsize_t) :: data_dim(3) = 0  ! Dimensions. Not used for groups. EG: Scaler data has [1, 0, 0].
  integer(size_t) :: data_size = -1    ! Size of datums. Not used for groups. For strings size = # of characters.
  integer :: num_attributes = -1       ! Number of associated attributes. Used for groups and datasets only.
end type

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_dataset_int (group_id, dset_name, int_buf, error, err_str)
!
! Routine to read a dataset of integers.
! This overloads hdf5_read_dataset_int_rankN where N is 0, 1, 2, or 3.
! Note: The int_buf array size and shape must be correct for the dataset read.
!
! Input:
!   loc_id          -- integer(hid_t): Id of group containing the dataset.
!   dset_name       -- character(*): Name of the dataset.
!   err_str         -- character(*), optional: String to use with error message.
!
! Output:
!   int_buf         -- integer: For datasets storing a single value.
!   int_buf(:)      -- integer: For datasets storing a 1D array.
!   int_buf(:,:)    -- integer: For datasets storing a 2D array.
!   int_buf(:,:,:)  -- integer: For datasets storing a 3D array.
!   error           -- logical: Set True if there is an error.
!-

interface hdf5_read_dataset_int
  module procedure hdf5_read_dataset_int_rank0
  module procedure hdf5_read_dataset_int_rank1
  module procedure hdf5_read_dataset_int_rank2
  module procedure hdf5_read_dataset_int_rank3
end interface

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_dataset_real (group_id, dset_name, real_buf, error, err_str)
!
! Routine to read a dataset of reals.
! This overloads hdf5_read_dataset_real_rankN where N is 0, 1, 2, or 3.
! Note: The real_buf array size and shape must be correct for the dataset read.
!
! Input:
!   loc_id          -- integer(hid_t): Id of group containing the dataset.
!   dset_name       -- character(*): Name of the dataset.
!   err_str         -- character(*), optional: String to use with error message.
!
! Output:
!   real_buf         -- real(rp): For datasets storing a single value.
!   real_buf(:)      -- real(rp): For datasets storing a 1D array.
!   real_buf(:,:)    -- real(rp): For datasets storing a 2D array.
!   real_buf(:,:,:)  -- real(rp): For datasets storing a 3D array.
!   error            -- logical: Set True if there is an error.
!-

interface hdf5_read_dataset_real
  module procedure hdf5_read_dataset_real_rank0
  module procedure hdf5_read_dataset_real_rank1
  module procedure hdf5_read_dataset_real_rank2
  module procedure hdf5_read_dataset_real_rank3
end interface

!------------------------------------------------------------------------------------------

interface hdf5_read_attribute_real
  module procedure hdf5_read_attribute_real_rank0
  module procedure hdf5_read_attribute_real_rank1
end interface

interface hdf5_read_attribute_int
  module procedure hdf5_read_attribute_int_rank0
  module procedure hdf5_read_attribute_int_rank1
end interface

interface hdf5_read_attribute_string
  module procedure hdf5_read_attribute_string_rank0
  module procedure hdf5_read_attribute_string_rank1
end interface

interface hdf5_read_attribute_alloc_string
  module procedure hdf5_read_attribute_alloc_string_rank0
end interface

interface hdf5_write_dataset_int
  module procedure hdf5_write_dataset_int_rank0
  module procedure hdf5_write_dataset_int_rank1
  module procedure hdf5_write_dataset_int_rank2
  module procedure hdf5_write_dataset_int_rank3
end interface

interface hdf5_write_dataset_real
  module procedure hdf5_write_dataset_real_rank0
  module procedure hdf5_write_dataset_real_rank1
  module procedure hdf5_write_dataset_real_rank2
  module procedure hdf5_write_dataset_real_rank3
end interface

interface hdf5_write_attribute_real
  module procedure hdf5_write_attribute_real_rank0
  module procedure hdf5_write_attribute_real_rank1
end interface

interface hdf5_write_attribute_int
  module procedure hdf5_write_attribute_int_rank0
  module procedure hdf5_write_attribute_int_rank1
end interface

interface hdf5_write_attribute_string
  module procedure hdf5_write_attribute_string_rank0
  module procedure hdf5_write_attribute_string_rank1
end interface

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_attribute_alloc_string_rank0(root_id, attrib_name, string, error, print_error, err_str)
!
! Routine to read a string attribute.
! Also see: hdf5_read_attribute_string
!
! Input:
!   root_id       -- integer(hid_t): ID of group or dataset containing the attribute.
!   attrib_name   -- character(*): Name of the attribute.
!   print_error   -- logical: If true, print an error message if there is a problem.
!   err_str       -- character(*), optional: String to use with error message.
!
! Output:
!   error         -- logical: Set true if there is an error. False otherwise.
!   string        -- character(:), allocatable: Variable length string to hold the attribute value.
!-

subroutine hdf5_read_attribute_alloc_string_rank0(root_id, attrib_name, string, error, print_error, err_str)

type (hdf5_info_struct) info

integer(hid_t) root_id, a_id
integer h5_err

logical error, print_error

character(*) attrib_name
character(:), allocatable :: string
character(*), optional :: err_str
character(*), parameter :: r_name = 'hdf5_read_attribute_alloc_string_rank0'

!

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)

if (info%data_class_type /= H5T_STRING_F) then
  if (print_error) then
    call out_io (s_error$, r_name, 'ATTRIBUTE: ' // string_option(attrib_name, err_str), 'IS NOT A STRING!')
  endif
  allocate(character(0) :: string)
  return
endif

allocate(character(info%data_size) :: string)
call H5LTget_attribute_string_f(root_id, '.', attrib_name, string, h5_err)
if (h5_err < 0) then
  if (print_error) then
    call out_io (s_error$, r_name, 'CANNOT READ ATTRIBUTE: ' // string_option(attrib_name, err_str))
  endif
  string = ''
  return
endif

error = .false.

end subroutine hdf5_read_attribute_alloc_string_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_attribute_string_rank0(root_id, attrib_name, string, error, print_error, err_str)
!
! Routine to read a string attribute.
! Also see: hdf5_read_attribute_alloc_string
!
! Input:
!   root_id       -- integer(hid_t): ID of group or dataset containing the attribute.
!   attrib_name   -- character(*): Name of the attribute.
!   print_error   -- logical: If true, print an error message if there is a problem.
!   err_str       -- character(*), optional: String to use with error message.
!
! Output:
!   error         -- logical: Set true if there is an error. False otherwise.
!   string        -- character(*): String to hold the attribute value. Set to blank if there is an error.
!-

subroutine hdf5_read_attribute_string_rank0(root_id, attrib_name, string, error, print_error, err_str)

type (hdf5_info_struct) info

integer(hid_t) root_id, a_id
integer h5_err, n

logical error, print_error

character(*) attrib_name
character(*) :: string
character(*), optional :: err_str
character(*), parameter :: r_name = 'hdf5_read_attribute_string'

!

string = ''

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)

if (info%data_class_type /= H5T_STRING_F) then
  if (print_error) then
    call out_io (s_error$, r_name, 'ATTRIBUTE: ' // string_option(attrib_name, err_str), 'IS NOT A STRING!')
  endif
  return
endif

call H5LTget_attribute_string_f(root_id, '.', attrib_name, string, h5_err)
if (h5_err < 0) then
  if (print_error) then
    call out_io (s_error$, r_name, 'CANNOT READ ATTRIBUTE: ' // string_option(attrib_name, err_str))
  endif
  return
endif

! This is to get around an HDF5 (V1.10.4) bug where extra garbage characters can be present.

n = min(len(string), info%data_size)
string = string(1:n)

error = .false.

end subroutine hdf5_read_attribute_string_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_attribute_string_rank1(root_id, attrib_name, string, error, print_error, err_str)
!
! Routine to create an HDF5 string array attribute.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group or dataset the attribute is to be put in.
!   attrib_name   -- character(*): Name of the attribute.
!   print_error   -- logical: If true, print an error message if there is a problem.
!   err_str       -- character(*), optional: String to use with error message.
!
! Output
!   strings(:)    -- character(*), allocatable: String array.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_read_attribute_string_rank1(root_id, attrib_name, strings, error, print_error, err_str)

type (hdf5_info_struct) info
type(c_ptr) :: f_ptr

integer(hid_t) :: root_id, attrib_id, attrib_type, space, mem_type
integer(size_t) size, string_len
integer(hsize_t) dims(1), maxdims(1)
integer h5_err

character(*) :: attrib_name
character(*), optional :: err_str
character(*), allocatable, target :: strings(:)
character(*), parameter :: r_name = 'hdf5_read_attribute_string_rank1'
logical error, print_error

!

error = .true.

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)

if (info%data_class_type /= H5T_STRING_F) then
  if (print_error) then
    call out_io (s_error$, r_name, 'ATTRIBUTE: ' // string_option(attrib_name, err_str), 'IS NOT A STRING ARRAY!')
  endif
  return
endif

call H5Aopen_f(root_id, attrib_name, attrib_id, h5_err); if (h5_err < 0) return
 
call H5Aget_type_f(attrib_id, attrib_type, h5_err)
call H5Tget_size_f(attrib_type, size, h5_err)

! Get dataspace and allocate memory for read buffer.

call H5Aget_space_f(attrib_id, space, h5_err)
call H5Sget_simple_extent_dims_f(space, dims, maxdims, h5_err)

if (allocated(strings)) deallocate(strings)
allocate(strings(1:dims(1)))

string_len = len(strings(1))
if (size > string_len+1) then
   call out_io (s_error$, r_name, 'STRINGS ARGUMENT LENGTH IS TO SMALL.')
   return
endif

! Create the memory datatype.

call H5Tcopy_f(H5T_FORTRAN_S1, mem_type, h5_err)
call H5Tset_size_f(mem_type, string_len, h5_err)

! Read the data.

f_ptr = C_LOC(strings(1)(1:1))
call H5Aread_f(attrib_id, mem_type, f_ptr, h5_err)

! Close and release resources.

call H5Aclose_f(attrib_id, h5_err)
call H5Sclose_f(space, h5_err)
call H5Tclose_f(mem_type, h5_err)

error = .false.

end subroutine hdf5_read_attribute_string_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_attribute_string_rank0 (root_id, attrib_name, string, error)
!
! Routine to create an HDF5 attribute whose value is a string.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group or dataset the attribute is to be put in.
!   attrib_name   -- character(*): Name of the attribute.
!   string        -- character(*): String attribute value.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_attribute_string_rank0 (root_id, attrib_name, string, error)

integer(hid_t) :: root_id
character(*) :: attrib_name, string
integer h5_err
logical error
!
error = .true.
call H5LTset_attribute_string_f(root_id, '.', attrib_name, trim(string), h5_err); if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_attribute_string_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_attribute_string_rank1(root_id, attrib_name, string, error)
!
! Routine to create an HDF5 string array attribute.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group or dataset the attribute is to be put in.
!   attrib_name   -- character(*): Name of the attribute.
!   string(:)     -- character(*): String array.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_attribute_string_rank1(root_id, attrib_name, string, error)

integer(hid_t) :: root_id, dspace_id, atype, attr_id
character(*) :: attrib_name, string(:)
character(:), allocatable :: str
integer h5_err, i, j, nn, n_len, n0, n1
integer(hsize_t) dim(1)
integer(size_t) sze
logical error

!

error = .true.
n_len = maxval(len_trim(string))
nn = size(string) * n_len
allocate(character(nn):: str)

do i = 1, size(string)
  n0 = (i - 1) * n_len
  n1 = len_trim(string(i))
  str(n0+1:n0+n1) = trim(string(i))
  do j = n0+n1+1, n0+n_len
    str(j:j) = c_null_char
  enddo
enddo

dim(1) = size(string)
call H5Screate_simple_f(1, dim, dspace_id, h5_err)
call H5Tcopy_f(H5T_C_S1, atype, h5_err)
sze = n_len
call H5Tset_size_f(atype, sze, h5_err)
call H5Acreate_f(root_id, attrib_name, atype, dspace_id, attr_id, h5_err, H5P_DEFAULT_F, H5P_DEFAULT_F)
call H5Awrite_f(attr_id, atype, str, dim, h5_err)

call H5Aclose_f(attr_id, h5_err)
call H5Tclose_f(atype, h5_err)
call H5Sclose_f(dspace_id, h5_err)

error = .false.

end subroutine hdf5_write_attribute_string_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_attribute_int_rank0(root_id, attrib_name, ival, error)
!
! Routine to create an attribute with a scalar integer value.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group or dataset the attribute is to be put in.
!   attrib_name   -- character(*): Name of the attribute.
!   ival          -- integer: Integer value of the attribute.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_attribute_int_rank0(root_id, attrib_name, ival, error)

integer(hid_t) :: root_id
character(*) :: attrib_name
integer :: ival
integer h5_err
logical error
!
error = .true.
call H5LTset_attribute_int_f(root_id, '.', attrib_name, [ival], 1_size_t, h5_err); if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_attribute_int_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_attribute_int_rank1(root_id, attrib_name, ival, error)
!
! Routine to create an attribute with a vector integer value.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group or dataset the attribute is to be put in.
!   attrib_name   -- character(*): Name of the attribute.
!   ival(:)       -- integer: Integer array attribute value.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_attribute_int_rank1(root_id, attrib_name, ival, error)

integer(hid_t) :: root_id
integer(size_t) iz 
character(*) :: attrib_name
integer :: ival(:)
integer h5_err
logical error
!
error = .true.
iz = size(ival)
call H5LTset_attribute_int_f(root_id, '.', attrib_name, ival, iz, h5_err); if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_attribute_int_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_attribute_real_rank0(root_id, attrib_name, rval, error)
!
! Routine to create an attribute with a scalar real value.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group or dataset the attribute is to be put in.
!   attrib_name   -- character(*): Name of the attribute.
!   rval          -- real(rp): real value of the attribute.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_attribute_real_rank0(root_id, attrib_name, rval, error)

integer(hid_t) :: root_id
character(*) :: attrib_name
real(rp) :: rval
integer h5_err
logical error
!
error = .true.
call H5LTset_attribute_double_f(root_id, '.', attrib_name, [rval], 1_size_t, h5_err); if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_attribute_real_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_attribute_real_rank1(root_id, attrib_name, rval, error)
!
! Routine to create an attribute with a real vector value.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group or dataset the attribute is to be put in.
!   attrib_name   -- character(*): Name of the attribute.
!   rval(:)       -- real(rp): real vector value of the attribute.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_attribute_real_rank1(root_id, attrib_name, rval, error)

integer(hid_t) :: root_id
integer(size_t) iz 
character(*) :: attrib_name
real(rp) :: rval(:)
integer h5_err
logical error
!
error = .true.
iz = size(rval)
call H5LTset_attribute_double_f(root_id, '.', attrib_name, rval, iz, h5_err); if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_attribute_real_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_open_file (file_name, action, file_id, error, verbose)
!
! Routine to open an HDF5 file.
!
! Note: To close the file When finished use:
!   call h5fclose_f(file_name, h5_err)  ! h5_err is an integer
!
! Input:
!   file_name   -- character(*): Name of the file
!   action      -- character(*): Possibilities are:
!                     'READ'    -- Read only.
!                     'WRITE'   -- New file for writing to.
!                     'APPEND'  -- If file exists, open file for reading/writing. 
!                                  If file does not exist, create new file.
!   verbose     -- logical, optional: Default False. If set True, toggle verbose output on.
!
! Output:
!   file_id     -- integer(hid_t): File handle.
!   error       -- logical: Set True if there is an error. False otherwise.
!-

subroutine hdf5_open_file (file_name, action, file_id, error, verbose)

integer(hid_t) file_id
integer h5_err, h_err

logical error, exist
logical, optional :: verbose

character(*) file_name, action
character(*), parameter :: r_name = 'hdf5_open_file'
character(200) full_name

!

error = .true.

call fullfilename(file_name, full_name, exist)
if (.not. exist) then
  call out_io (s_error$, r_name, 'MALFORMED FILE NAME OR ENVIRONMENT VARIABLE NOT DEFINED: ' // file_name)
  return
endif

call h5open_f(h5_err)         ! Init Fortran interface.

if (logic_option(.false., verbose)) then
  call H5Eset_auto_f(1, h5_err)   ! Verbose
else
  call H5Eset_auto_f(0, h5_err)   ! Run silent
endif

select case (action)
case ('READ')
  call h5fopen_f(full_name, H5F_ACC_RDONLY_F, file_id, h5_err)

case ('WRITE')
  call h5fcreate_f (full_name, H5F_ACC_TRUNC_F, file_id, h5_err)

case ('APPEND')
  inquire (file = full_name, exist = exist)
  if (exist) then
    call h5fopen_f(full_name, H5F_ACC_RDWR_F, file_id, h5_err)
  else
    call h5fcreate_f (full_name, H5F_ACC_TRUNC_F, file_id, h5_err)
  endif

case default
  call out_io(s_fatal$, r_name, 'BAD ACTION ARGUMENT! ' // quote(action))
  stop
end select

call H5Eset_auto_f(1, h_err)    ! Reset
call h5eclear_f(h_err)

if (h5_err < 0) then
  select case (action)
  case ('READ');    call out_io (s_error$, r_name, 'CANNOT OPEN FILE FOR READING: ' // file_name)
  case ('WRITE');   call out_io (s_error$, r_name, 'CANNOT CREATE FILE FOR WRITING: ' // file_name)
  case ('APPEND');  call out_io (s_error$, r_name, 'CANNOT OPEN FILE FOR APPENDING: ' // file_name)
  end select
  return
endif

error = .false.

end subroutine hdf5_open_file

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Function hdf5_open_object(root_id, object_name, info, error, print_error) result (obj_id)
!
! Routine to open an existing group or dataset.
!
! Note: Use hdf5_close_object to close the object.
! Also see: hdf5_open_group and hdf5_open_dataset.
!
! Input:
!   root_id     -- integer(hid_t): ID of the group containing the object to be opened.
!   object_name -- character(*): Name of the object to be opened
!   info        -- hdf5_info_struct: Information on the object.
!   print_error -- logical: Print an error message if there is an error?
!
! Output:
!   error       -- logical: Set True if there is an error. False otherwise.
!   obj_id      -- integer(hid_t): Object ID.
!-

function hdf5_open_object(root_id, object_name, info, error, print_error) result (obj_id)

type (hdf5_info_struct) info

integer(hid_t) root_id, obj_id
integer h5_err

logical error, print_error

character(*) object_name
character(*), parameter :: r_name = 'hdf5_open_object'

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  obj_id = hdf5_open_dataset (root_id, object_name, error, print_error) 
elseif (info%element_type == H5O_TYPE_GROUP_F) then
  obj_id = hdf5_open_group(root_id, object_name, error, print_error) 
else
  error = .true.
  obj_id = 0
  if (print_error) then
    call out_io (s_error$, r_name, 'UNKNOWN OBJECT TYPE FOR: ' // quote(object_name))
  endif
endif

end function hdf5_open_object

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_close_object(obj_id, info)
!
! Routine to close a group or dataset.
!
! Note: Use hdf5_open_object to open the object.
!
! Input:
!   obj_id      -- integer(hid_t): Object ID.
!   info        -- hdf5_info_struct: Information on the object. 
!                     Obtained when hdf5_open_object was called.
!-

subroutine hdf5_close_object(obj_id, info)

type (hdf5_info_struct) info

integer(hid_t) obj_id
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
!+
! Function hdf5_exists (root_id, object_name, error, print_error) result (exists)
!
! Routine to check if a object with object_name exists relative to root_id.
!
! Input:
!   root_id     -- integer(hid_t): ID of the base grroup.
!   object_name -- character(*): Path of the object.
!   print_error   -- logical: If true, print an error message if there is a problem.
!
! Output:
!   error         -- logical: Set true if there is an error. For example, if any element in the path 
!                     of object_name, except for the target, does not exist.
!   exists        -- logical: Object exists.
!-

function hdf5_exists (root_id, object_name, error, print_error) result (exists)

integer(hid_t) root_id
integer h5_err

logical error, print_error
logical exists

character(*) object_name
character(*), parameter :: r_name = 'hdf5_exists'

!

call H5Lexists_f(root_id, trim(object_name), exists, h5_err, H5P_DEFAULT_F)
error = (h5_err /= 0)
if (error .and. print_error) then
  call out_io (s_error$, r_name, 'CANNOT QUERY EXISTANCE: ' // quote(object_name))
endif

end function hdf5_exists

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Function hdf5_open_group (root_id, group_name, error, print_error) result (g_id)
!
! Rouine to open an existing group.
!
! Notes: 
!   Use H5Gclose_f to close the group.
!   Use H5Gcreate_f to create a new group.
!
! Also see: hdf5_open_object and hdf5_open_dataset.
!
! Input:
!   root_id     -- integer(hid_t): ID of the Parent group containing the group to be opened.
!   group_name  -- character(*): Name of the group to be opened.
!   print_error -- logical: Print an error message if there is an error?
!
! Output:
!   error       -- logical: Set True if there is an error. False otherwise.
!   g_id        -- integer(hid_t): Group ID.
!-

function hdf5_open_group (root_id, group_name, error, print_error) result (g_id)

integer(hid_t) root_id, g_id
integer h5_err

logical error, print_error
logical exists

character(*) group_name
character(*), parameter :: r_name = 'hdf5_open_group'

! H5Lexists_f is not smart enough to recognize that group_name = "." or "./" is 
! equivalent to the root group.

error = .false.

if (group_name == '.' .or. group_name == './') then
  g_id = root_id
  return
endif

!

call H5Lexists_f(root_id, group_name, exists, h5_err, H5P_DEFAULT_F)
if (.not. exists) then
  if (print_error) then
    call out_io (s_error$, r_name, 'GROUP DOES NOT EXIST: ' // quote(group_name))
  endif
  error = .true.
  g_id = 0
  return
endif
 
call H5Gopen_f (root_id, group_name, g_id, h5_err, H5P_DEFAULT_F)
if (h5_err == -1) return

end function hdf5_open_group

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Function hdf5_group_n_links (root_id, error, storage_type) result (n_links)
!
! Routine to return the number of objects in a group.
!
! Use this with hdf5_get_object by index to interate over all objects in a group.
!
! Input:
!   root_id         -- integer(hid_t): ID of the group.
!
! Output:
!   error           -- logical: Set True if there is an error. False otherwise.
!   storage_type    -- integer, optional: Type of storage for the object links.
!                           H5G_STORAGE_TYPE_COMPACT_F: Compact storage
!                           H5G_STORAGE_TYPE_DENSE_F: Indexed storage
!                           H5G_STORAGE_TYPE_SYMBOL_TABLE_F: Symbol tables
!   n_links         -- integer: number of links.
!-

function hdf5_group_n_links (root_id, error, storage_type) result (n_links)

integer (hid_t) root_id
integer n_links, store_typ, h5_err, max_corder
integer, optional :: storage_type

logical error

!

call H5Gget_info_f(root_id, store_typ, n_links, max_corder, h5_err)
error = (h5_err /= 0)
if (present(storage_type)) storage_type = store_typ

end function hdf5_group_n_links

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_get_object_by_index (root_id, idx, obj_name, obj_info, error)
!
! Routine get an object's id, name, and other info given the group the object is in and its index number.
! This routine can be used to iterate over all the objects of a group.
!
! Use hdf5_group_n_links to get the number of objects. Objects in the group have indexes 0 through n_links-1.
! When iterating over the objects in the group, the objects will be sorted by object name in increasing name order.
!
! Use hdf5_open_object to open the object.
!
! Input:
!   root_id   -- integer(hid_t): ID of the group containing the object
!   idx       -- integer: Index number for the object.
!
! Output:
!   ob_name   -- character(*): Name of the object.
!   obj_info  -- hdf5_info_struct: Object information.
!   error     -- logical: Set true if there is an error getting the object information.
!-

subroutine hdf5_get_object_by_index (root_id, idx, obj_name, obj_info, error)

type (hdf5_info_struct) obj_info

integer(hid_t) root_id, obj_id
integer(size_t) g_size
integer(hsize_t) idxh
integer idx, h5_err

logical error
character(*) obj_name
character(200) c_name

!

idxh = idx
call H5Lget_name_by_idx_f (root_id, '.', H5_INDEX_NAME_F, H5_ITER_INC_F, idxh, c_name, h5_err, g_size)
error = (h5_err /= 0)
if (error) return

call to_f_str(c_name, obj_name)
obj_info = hdf5_object_info (root_id, obj_name, error, .true.)

end subroutine hdf5_get_object_by_index

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Function hdf5_open_dataset(root_id, dataset_name, error, print_error) result (ds_id)
!
! Routine to open an existing group or dataset.
!
! Note: Use H5Dclose_f close the dataset.
! Also see: hdf5_open_group and hdf5_open_object.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group containing the dataset to be opened.
!   dataset_name  -- character(*): Name of the dataset to be opened.
!   print_error   -- logical: Print an error message if there is an error?
!
! Output:
!   error         -- logical: Set True if there is an error. False otherwise.
!   ds_id         -- integer(hid_t): Dataset ID.
!-

function hdf5_open_dataset (root_id, dataset_name, error, print_error) result (ds_id)

integer(hid_t) root_id, ds_id
integer h5_err

logical error, print_error
logical exists

character(*) dataset_name
character(*), parameter :: r_name = 'hdf5_open_dataset'

!

error = .true.
call H5Lexists_f(root_id, dataset_name, exists, h5_err, H5P_DEFAULT_F)
if (.not. exists) then
  if (print_error) then
    call out_io (s_error$, r_name, 'DATASET DOES NOT EXIST: ' // quote(dataset_name))
  endif
  ds_id = 0
  return
endif
 
call H5Dopen_f (root_id, dataset_name, ds_id, h5_err, H5P_DEFAULT_F)
if (h5_err == -1) return
error = .false.

end function hdf5_open_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Function hdf5_num_attributes(root_id) result (num)
!
! Routine to return the number of attributes associated with a group or dataset.
!
! Also see: hdf5_get_attribute_by_index
!
! Input:
!   root_id     -- integer(hid_t): Group or dataset ID number.
!
! Output:
!   num         -- integer: Number of attributes in the group or dataset.
!-

function hdf5_num_attributes(root_id) result (num)

integer(hid_t) :: root_id
integer num, h5_err

!

call H5Aget_num_attrs_f (root_id, num, h5_err)

end function hdf5_num_attributes

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_get_attribute_by_index(root_id, attrib_indx, attrib_id, attrib_name)
!
! Routine to return the ID and name of an attribute given the attribute's index number.
! This routine is useful for looping over all the attributes in a group or dataset.
!
! Input:
!   root_id       -- integer(hid_t): ID number of the group or dataset containing the attribute.
!   attrib_indx   -- integer: Attribute index. Will be in the range 1 to hdf5_num_attributes.
!
! Output:
!   attrib_id     -- integer(hid_t): ID number of the attribute.
!   attrib_name   -- character(*): Name of the attribute.
!-

subroutine hdf5_get_attribute_by_index(root_id, attrib_indx, attrib_id, attrib_name)

integer(hid_t) root_id, attrib_id
integer(size_t) nam_len
integer attrib_indx, h5_err

character(*) attrib_name

!

call H5Aopen_by_idx_f (root_id, ".", H5_INDEX_CRT_ORDER_F, H5_ITER_INC_F, int(attrib_indx-1, hsize_t), &
                                                                 attrib_id, h5_err, aapl_id=H5P_DEFAULT_F)
nam_len = len(attrib_name)
call H5Aget_name_f(attrib_id, nam_len, attrib_name, h5_err)
call H5Aclose_f(attrib_id, h5_err)

end subroutine hdf5_get_attribute_by_index

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Function hdf5_attribute_info(root_id, attrib_name, error, print_error) result (info)
!
! Routine to return information on an attribute given the attribute name and encomposing group.
!
! Input:
!   root_id       -- integer(hid_t): ID of group or dataset containing the attribute.
!   attrib_name   -- character(*): Name of the attribute.
!   print_error   -- logical: If true, print an error message if there is a problem.
!
! Output:
!   error         -- logical: Set true if there is an error. False otherwise.
!   info          -- hdf5_info_struct: Information on the attribute.
!-

function hdf5_attribute_info(root_id, attrib_name, error, print_error) result (info)

type (hdf5_info_struct) info

integer(hid_t) root_id, a_id
integer h5_err

logical error, print_error, exists

character(*) attrib_name
character(*), parameter :: r_name = 'hdf5_attribute_info'

!

error = .true.
info = hdf5_info_struct()

call H5Aexists_f (root_id, attrib_name, exists, h5_err)
if (.not. exists .or. h5_err == -1) then
  if (print_error) then
    call out_io (s_error$, r_name, 'ATTRIBUTE IS NOT PRESENT: ' // attrib_name)
  endif
  return
endif

call H5LTget_attribute_info_f(root_id, '.', attrib_name, info%data_dim, info%data_class_type, info%data_size, h5_err)
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
!+
! Function hdf5_object_info (root_id, obj_name, error, print_error) result (info)
!
! Routine to get information on an object (group or dataset).
!
! Input:
!   root_id       -- integer(hid_t): ID of group containing the object in question.
!   obj_name      -- character(*): Name of the object.
!   print_error   -- logical: If true, print an error message if there is a problem.
!
! Output:
!   error         -- logical: Set true if there is an error. False otherwise.
!   info          -- hdf5_info_struct: Information on the object.
!-

function hdf5_object_info (root_id, obj_name, error, print_error) result (info)

type (hdf5_info_struct) info
type (H5O_info_t) :: infobuf 

integer(hid_t), value :: root_id
integer stat, h5_err

character(*) obj_name

logical error, print_error

!

error = .true.

call H5Oget_info_by_name_f(root_id, trim(obj_name), infobuf, h5_err)
info%element_type = infobuf%type
info%num_attributes = infobuf%num_attrs

if (info%element_type == H5O_TYPE_DATASET_F) then
  call H5LTget_dataset_info_f(root_id, trim(obj_name), info%data_dim, info%data_class_type, info%data_size, h5_err)
endif

error = .false.

end function hdf5_object_info

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_attribute_int_rank0(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)
!
! Routine to read an scaler (rank 0) integer attribute value.
! Overloaded by: hdf5_read_attribute_int
!
! Input:
!   root_id       -- integer(hid_t): ID of group or dataset containing the attribute.
!   attrib_name   -- character(*): Name of the attribute.
!   print_error   -- logical: If true, print an error message if there is a problem.
!   err_str       -- character(*), optional: String to use with error message.
!   dflt_value    -- integer, optional: Default value if there is an error. 
!                         If not present, the default is 0.
!
! Output:
!   error         -- logical: Set true if there is an error. False otherwise.
!   attrib_value  -- integer: Value of the attribute.
!-

subroutine hdf5_read_attribute_int_rank0(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)

integer(hid_t) root_id
integer attrib_value, a_val(1)
integer, optional :: dflt_value

logical error, print_error

character(*) attrib_name
character(*), optional :: err_str

call hdf5_read_attribute_int_rank1(root_id, attrib_name, a_val, error, print_error, err_str, dflt_value)
attrib_value = a_val(1)

end subroutine hdf5_read_attribute_int_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_attribute_int_rank1(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)
!
! Routine to read a vector (rank 1) integer attribute array.
! Overloaded by: hdf5_read_attribute_int
!
! Input:
!   root_id         -- integer(hid_t): ID of group or dataset containing the attribute.
!   attrib_name     -- character(*): Name of the attribute.
!   print_error     -- logical: If true, print an error message if there is a problem.
!   err_str       -- character(*), optional: String to use with error message.
!   dflt_value      -- integer, optional: Default value if there is an error. 
!                         If not present, the default is 0.
!
! Output:
!   error           -- logical: Set true if there is an error. False otherwise.
!   attrib_value(:) -- integer: Value of the attribute.
!-

subroutine hdf5_read_attribute_int_rank1(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)

type (hdf5_info_struct) info

integer(hid_t) root_id, a_id
integer attrib_value(:)
integer, optional :: dflt_value
integer h5_err

logical error, print_error

character(*) attrib_name
character(*), optional :: err_str
character(*), parameter :: r_name = 'hdf5_read_attribute_int_rank1'

!

attrib_value = integer_option(0, dflt_value)

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)

if (info%data_class_type == H5T_INTEGER_F) then
  call H5LTget_attribute_int_f(root_id, '.', attrib_name, attrib_value, h5_err)
else
  if (print_error) call out_io (s_error$, r_name, 'ATTRIBUTE IS NOT OF INTEGER TYPE: ' // string_option(attrib_name, err_str))
  return
endif

error = .false.

end subroutine hdf5_read_attribute_int_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_attribute_real_rank0(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)
!
! Routine to read an scaler (rank 0) real attribute value.
! Overloaded by: hdf5_read_attribute_real
!
! Input:
!   root_id       -- integer(hid_t): ID of group or dataset containing the attribute.
!   attrib_name   -- character(*): Name of the attribute.
!   print_error   -- logical: If true, print an error message if there is a problem.
!   err_str       -- character(*), optional: String to use with error message.
!   dflt_value    -- real(rp), optional: Default value if there is an error. 
!                         If not present, the default is 0.
!
! Output:
!   error         -- logical: Set true if there is an error. False otherwise.
!   attrib_value  -- real(rp): Value of the attribute.
!-

subroutine hdf5_read_attribute_real_rank0(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)

integer(hid_t) root_id
integer h5_err

real(rp) attrib_value, val(1)
real(rp), optional :: dflt_value

logical error, print_error

character(*) attrib_name
character(*), optional :: err_str

!

call hdf5_read_attribute_real_rank1(root_id, attrib_name, val, error, print_error, err_str, dflt_value)
attrib_value = val(1)

end subroutine hdf5_read_attribute_real_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_attribute_real_rank1(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)
!
! Routine to read a vector (rank 1) real attribute array
! Overloaded by: hdf5_read_attribute_real
!
! Input:
!   root_id         -- integer(hid_t): ID of group or dataset containing the attribute.
!   attrib_name     -- character(*): Name of the attribute.
!   print_error     -- logical: If true, print an error message if there is a problem.
!   err_str         -- character(*), optional: String to use with error message.
!   dflt_value      -- real(rp), optional: Default value if there is an error. 
!                         If not present, the default is 0.
!
! Output:
!   error           -- logical: Set true if there is an error. False otherwise.
!   attrib_value(:) -- real(rp): Value array of the attribute.
!-

subroutine hdf5_read_attribute_real_rank1(root_id, attrib_name, attrib_value, error, print_error, err_str, dflt_value)

type (hdf5_info_struct) info

integer(hid_t) root_id, a_id
integer h5_err

real(rp) attrib_value(:)
real(rp), optional :: dflt_value

logical error, print_error

character(*) attrib_name
character(*), optional :: err_str
character(*), parameter :: r_name = 'hdf5_read_attribute_real_rank1'

!

attrib_value = real_option(0.0_rp, dflt_value)
error = .true.

info = hdf5_attribute_info(root_id, attrib_name, error, print_error)
if (error) return

if (info%data_class_type == H5T_INTEGER_F .or. info%data_class_type == H5T_FLOAT_F) then
  call H5LTget_attribute_double_f(root_id, '.', attrib_name, attrib_value, h5_err)
else
  if (print_error) call out_io (s_error$, r_name, 'ATTRIBUTE IS NOT OF REAL TYPE: ' // string_option(attrib_name, err_str))
  return
endif

error = .false.

end subroutine hdf5_read_attribute_real_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_real_rank0(root_id, dataset_name, value, error)
!
! Routine to create a dataset with one real value.
! Overloaded by: interface hdf5_write_dataset_real
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value         -- real(rp): Dataset value.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_real_rank0 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(1)
integer h5_err
real(rp) value
real(rp) vector(1)
logical error
character(*) dataset_name

!

error = .true.
v_size = 1
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, [v_size], vector, h5_err);  if (h5_err < 0) return
value = vector(1)
error = .false.

end subroutine hdf5_write_dataset_real_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_real_rank1(root_id, dataset_name, value, error)
!
! Routine to create a dataset with an array of real values.
! Overloaded by: interface hdf5_write_dataset_real
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value(:)      -- real(rp): Dataset value array.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_real_rank1 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(1)
integer h5_err
real(rp) value(:)
logical error
character(*) dataset_name

!

error = .true.
v_size = size(value)
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, v_size, value, h5_err);  if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_dataset_real_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_real_rank2(root_id, dataset_name, value, error)
!
! Routine to create a dataset with a matrix of real values.
! Overloaded by: interface hdf5_write_dataset_real
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value(:,:)    -- real(rp): Dataset value matrix.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_real_rank2 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(2), ds_id
integer h5_err
real(rp) value(:,:)
logical error
character(*) dataset_name

!

error = .true.
v_size = [size(value, 1), size(value, 2)]
call H5LTmake_dataset_double_f(root_id, dataset_name, 2, v_size, value, h5_err);  if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_dataset_real_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_real_rank3(root_id, dataset_name, value, error)
!
! Routine to create a dataset with a 3D array of real values.
! Overloaded by: interface hdf5_write_dataset_real
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value(:,:,:)  -- real(rp): Dataset values
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_real_rank3 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(3)
integer h5_err
real(rp) value(:,:,:)
logical error
character(*) dataset_name

!

error = .true.
v_size = [size(value, 1), size(value, 2), size(value, 3)]
call H5LTmake_dataset_double_f(root_id, dataset_name, 3, v_size, value, h5_err);  if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_dataset_real_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_int_rank0(root_id, dataset_name, value, error)
!
! Routine to create a dataset with one integer value.
! Overloaded by: interface hdf5_write_dataset_int
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value         -- integer: Dataset value.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_int_rank0 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(1)
integer h5_err
integer value
integer vector(1)
logical error
character(*) dataset_name

!

error = .true.
v_size = 1
call H5LTmake_dataset_int_f(root_id, dataset_name, 1, v_size, vector, h5_err);  if (h5_err < 0) return
value = vector(1)
error = .false.

end subroutine hdf5_write_dataset_int_rank0

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_int_rank1(root_id, dataset_name, value, error)
!
! Routine to create a dataset with an array of integer values.
! Overloaded by: interface hdf5_write_dataset_int
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value(:)      -- integer: Dataset value array.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_int_rank1 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(1)
integer h5_err
integer value(:)
logical error
character(*) dataset_name

!

error = .true.
v_size = size(value)
call H5LTmake_dataset_int_f(root_id, dataset_name, 1, v_size, value, h5_err);  if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_dataset_int_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_int_rank2(root_id, dataset_name, value, error)
!
! Routine to create a dataset with a matrix of integer values.
! Overloaded by: interface hdf5_write_dataset_int
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value(:,:)    -- integer: Dataset value matrix.
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_int_rank2 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(2)
integer h5_err
integer value(:,:)
logical error
character(*) dataset_name

!

error = .true.
v_size = [size(value, 1), size(value,2)]
call H5LTmake_dataset_int_f(root_id, dataset_name, 2, v_size, value, h5_err);  if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_dataset_int_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_write_dataset_int_rank3(root_id, dataset_name, value, error)
!
! Routine to create a dataset with a 3D array of integer values.
! Overloaded by: interface hdf5_write_dataset_int
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the dataset is to be put in.
!   dataset_name  -- character(*): Name of the dataset.
!   value(:,:,:)  -- integer: Dataset values
!   error         -- logical Set True if there is an error. False otherwise.
!-

subroutine hdf5_write_dataset_int_rank3 (root_id, dataset_name, value, error)

integer(hid_t) root_id, v_size(3)
integer h5_err
integer value(:,:,:)
logical error
character(*) dataset_name

!

error = .true.
v_size = [size(value, 1), size(value, 2), size(value, 3)]
call H5LTmake_dataset_int_f(root_id, dataset_name, 3, v_size, value, h5_err);  if (h5_err < 0) return
error = .false.

end subroutine hdf5_write_dataset_int_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_read_dataorder(root_id, object_name, axislabels, d_ord)
!
! Routine to read the axisLabels (or old gridDataOrder) attribute for arrays with dimension higher than 1.
!
! Input:
!   root_id       -- integer(hid_t): ID of the group the object is to be put in.
!   object_name   -- character(*): Name of the object containing the data.
!   axislabels(:) -- character(*): Logical axis labels.
!
! Output:
!   d_ord         -- character(*): gridDataOrder value. Set to:
!                     ''    If neither gridDataOrder nor axisLabes attributes exist.
!                     'C'   C/C++ ordering.
!                     'F'   Fortran ordering.
!-

subroutine hdf5_read_dataorder(root_id, object_name, axislabels, d_ord)

type (hdf5_info_struct) info
integer(hid_t) root_id, z_id
integer n
character(*) object_name, d_ord
character(*) :: axislabels(:)
character(16), allocatable :: physical_labels(:)
logical error

! First check for old-style gridDataOrder 

d_ord = ''
info = hdf5_object_info (root_id, object_name, error, .true.);  if (error) return
z_id = hdf5_open_object (root_id, object_name, info, error, .true.);  if (error) return
call hdf5_read_attribute_string(z_id, 'gridDataOrder', d_ord, error, .false.)
call hdf5_close_object(z_id, info)
if (d_ord /= '') return

! Now check for axislabels

call hdf5_read_attribute_string(root_id, 'axisLabels', physical_labels, error, .false.); if (error) return
n = size(axislabels)
if (size(physical_labels) /= n) return

if (all(physical_labels == axislabels)) then
  d_ord = 'C'
  return
endif

if (all(physical_labels == axislabels(n:1:-1))) d_ord = 'F'
return

end subroutine hdf5_read_dataorder

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine hdf5_check_open (root_id)
!
! Routine to print information on what resources are still open.
! This routine is used for debugging.
! 
! Input:
!   root_id   -- integer(hid_t): ID of root node to check.
!-

subroutine hdf5_check_open (root_id)

integer(hid_t) root_id
integer(hid_t), allocatable :: obj_list(:)
integer(size_t) cnt, n_in, n_out
integer i, h5_err, typ
character(100) name

!

call H5Fget_obj_count_f(root_id, H5F_OBJ_ALL_F, cnt, h5_err)
print *, 'Open: ', cnt

allocate(obj_list(cnt))
call H5Fget_obj_ids_f(root_id, H5F_OBJ_ALL_F, cnt, obj_list, h5_err)

do i = 1, cnt
  call H5Iget_type_f(obj_list(i), typ, h5_err)
  n_in = len(name)
  call H5Iget_name_f(obj_list(i), name, n_in, n_out, h5_err)
  print *, '    Open: ', name(1:n_out)
enddo

end subroutine hdf5_check_open

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_double_kind_8_rank0

SUBROUTINE hdf5_read_dataset_real_rank0(loc_id, dset_name, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  REAL(KIND=8), INTENT(INout), TARGET :: buf
  INTEGER :: h5_err 
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  logical error
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_real_rank0'
  !
  f_ptr = C_LOC(buf)
  namelen = LEN(dset_name)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, H5T_NATIVE_DOUBLE, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)
END SUBROUTINE hdf5_read_dataset_real_rank0

!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_double_kind_8_rank1

SUBROUTINE hdf5_read_dataset_real_rank1(loc_id, dset_name, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  REAL(KIND=8), INTENT(INout), TARGET :: buf(:)
  REAL(KIND=8), target :: temp_buf(size(buf))
  INTEGER :: h5_err 
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  logical error
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_real_rank1'
  !
  f_ptr = C_LOC(temp_buf(1))
  namelen = LEN(dset_name)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, H5T_NATIVE_DOUBLE, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)
  buf = temp_buf
END SUBROUTINE hdf5_read_dataset_real_rank1

!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_double_kind_8_rank2

SUBROUTINE hdf5_read_dataset_real_rank2(loc_id, dset_name, axislabels, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  character(*) axislabels(:)
  REAL(KIND=8), INTENT(INout), TARGET :: buf(:,:)
  REAL(KIND=8), target, allocatable :: temp_buf(:,:)
  INTEGER :: h5_err, i
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  logical error
  character(1) d_ord
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_real_rank2'
  !
  call hdf5_read_dataorder(loc_id, dset_name, axislabels, d_ord)
  if (d_ord == 'C') then
    allocate (temp_buf(size(buf,2), size(buf,1)))
  else
    allocate (temp_buf(size(buf,1), size(buf,2)))
  endif

  f_ptr = C_LOC(temp_buf(1,1))
  namelen = LEN(dset_name)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, H5T_NATIVE_DOUBLE, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)

  if (d_ord == 'C') then
    do i = 1, size(buf,1)
      buf(i,:) = temp_buf(:,i)
    enddo
  else
    buf = temp_buf
  endif
END SUBROUTINE hdf5_read_dataset_real_rank2

!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_double_kind_8_rank3

SUBROUTINE hdf5_read_dataset_real_rank3(loc_id, dset_name, axislabels, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  character(*) axislabels(:)
  REAL(KIND=8), INTENT(INout), TARGET :: buf(:, :, :)
  REAL(KIND=8), target, allocatable :: temp_buf(:, :, :)
  INTEGER :: h5_err, i, j
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  logical error
  character(1) d_ord
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_real_rank3'
  !
  call hdf5_read_dataorder(loc_id, dset_name, axislabels, d_ord)
  if (d_ord == 'C') then
    allocate (temp_buf(size(buf, 3), size(buf,2), size(buf,1)))
  else
    allocate (temp_buf(size(buf,1), size(buf,2), size(buf, 3)))
  endif

  f_ptr = C_LOC(temp_buf(1,1,1))
  namelen = LEN(dset_name)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, H5T_NATIVE_DOUBLE, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)
  buf = temp_buf

  if (d_ord == 'C') then
    do i = 1, size(buf,1);  do j = 1, size(buf,2)
      buf(i,j,:) = temp_buf(:,j,i)
    enddo;  enddo
  else
    buf = temp_buf
  endif
END SUBROUTINE hdf5_read_dataset_real_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_int_kind_4_rank0

SUBROUTINE hdf5_read_dataset_int_rank0(loc_id, dset_name, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  INTEGER(KIND=4), INTENT(INout), TARGET :: buf
  INTEGER :: h5_err 
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  INTEGER(hid_t) :: type_id
  logical error
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_int_rank0'
  !
  f_ptr = C_LOC(buf)
  namelen = LEN(dset_name)
  type_id = h5kind_to_type(KIND(buf), H5_INTEGER_KIND)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)
END SUBROUTINE hdf5_read_dataset_int_rank0

!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_int_kind_4_rank1

SUBROUTINE hdf5_read_dataset_int_rank1(loc_id, dset_name, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  INTEGER(KIND=4), INTENT(INout), TARGET :: buf(:)
  INTEGER(KIND=4), target :: temp_buf(size(buf))
  INTEGER :: h5_err 
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  INTEGER(hid_t) :: type_id
  logical error
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_int_rank1'
  !
  f_ptr = C_LOC(temp_buf(1))
  namelen = LEN(dset_name)
  type_id = h5kind_to_type(KIND(buf(1)), H5_INTEGER_KIND)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)
  buf = temp_buf
END SUBROUTINE hdf5_read_dataset_int_rank1

!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_int_kind_4_rank2

SUBROUTINE hdf5_read_dataset_int_rank2(loc_id, dset_name, axislabels, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  character(*) axislabels(:)
  INTEGER(KIND=4), INTENT(INout), TARGET :: buf(:,:)
  INTEGER(KIND=4), target, allocatable :: temp_buf(:,:)
  INTEGER :: h5_err, i
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  INTEGER(hid_t) :: type_id
  logical error
  character(1) d_ord
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_int_rank2'
  !
  call hdf5_read_dataorder(loc_id, dset_name, axislabels, d_ord)
  if (d_ord == 'C') then
    allocate (temp_buf(size(buf,2), size(buf,1)))
  else
    allocate (temp_buf(size(buf,1), size(buf,2)))
  endif

  f_ptr = C_LOC(temp_buf(1,1))
  namelen = LEN(dset_name)
  type_id = h5kind_to_type(KIND(buf(1,1)), H5_INTEGER_KIND)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)
  buf = temp_buf

  if (d_ord == 'C') then
    do i = 1, size(buf,1)
      buf(i,:) = temp_buf(:,i)
    enddo
  else
    buf = temp_buf
  endif
END SUBROUTINE hdf5_read_dataset_int_rank2

!------------------------------------------------------------------------------------------
! Adapted from h5ltread_dataset_int_kind_4_rank3

SUBROUTINE hdf5_read_dataset_int_rank3(loc_id, dset_name, axislabels, buf, error, err_str)
  IMPLICIT NONE
  INTEGER(hid_t)  , INTENT(IN) :: loc_id
  CHARACTER(LEN=*), INTENT(IN) :: dset_name
  character(*) axislabels(:)
  INTEGER(KIND=4), INTENT(INout), TARGET :: buf(:,:,:)
  INTEGER(KIND=4), target, allocatable :: temp_buf(:,:,:)
  INTEGER :: h5_err, i, j
  TYPE(C_PTR) :: f_ptr
  INTEGER(size_t) :: namelen
  INTEGER(hid_t) :: type_id
  logical error
  character(1) d_ord
  character(*), optional :: err_str
  character(*), parameter :: r_name = 'hdf5_read_dataset_int_rank3'
  !
  call hdf5_read_dataorder(loc_id, dset_name, axislabels, d_ord)
  if (d_ord == 'C') then
    allocate (temp_buf(size(buf, 3), size(buf,2), size(buf,1)))
  else
    allocate (temp_buf(size(buf,1), size(buf,2), size(buf, 3)))
  endif

  f_ptr = C_LOC(temp_buf(1,1,1))
  namelen = LEN(dset_name)
  type_id = h5kind_to_type(KIND(buf(1,1,1)), H5_INTEGER_KIND)
  h5_err = h5ltread_dataset_c(loc_id, namelen, dset_name, type_id, f_ptr)
  error = (h5_err < 0)
  if (error .and. present(err_str)) call out_io (s_error$, r_name, 'CANNOT READ: ' // err_str)
  buf = temp_buf

  if (d_ord == 'C') then
    do i = 1, size(buf,1);  do j = 1, size(buf,2)
      buf(i,j,:) = temp_buf(:,j,i)
    enddo;  enddo
  else
    buf = temp_buf
  endif
END SUBROUTINE hdf5_read_dataset_int_rank3

end module
