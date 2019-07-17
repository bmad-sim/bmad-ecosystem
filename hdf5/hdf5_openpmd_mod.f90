module hdf5_openpmd_mod

use hdf5_interface

! Common units

type pmd_unit_struct
 character(8) :: unitSymbol = ''     ! Native units name. EG 'eV'
 real(rp) :: unitSI = 0              ! Conversion to SI
 real(rp) :: unitDimension(7) = 0    ! SI Base Exponents
end type

real(rp), parameter :: dim_1(7)              = [0,0,0,0,0,0,0]
real(rp), parameter :: dim_length(7)         = [1,0,0,0,0,0,0]
real(rp), parameter :: dim_mass(7)           = [0,1,0,0,0,0,0]
real(rp), parameter :: dim_time(7)           = [0,0,1,0,0,0,0]
real(rp), parameter :: dim_current(7)        = [0,0,0,1,0,0,0]
real(rp), parameter :: dim_temperture(7)     = [0,0,0,0,1,0,0]
real(rp), parameter :: dim_mol(7)            = [0,0,0,0,0,1,0]
real(rp), parameter :: dim_luminous(7)       = [0,0,0,0,0,0,1]
real(rp), parameter :: dim_charge(7)         = [0,0,1,1,0,0,0]
real(rp), parameter :: dim_elec_potential(7) = [1,1,-3,-1,0,0,0]
real(rp), parameter :: dim_velocity(7)       = [1,0,-1,0,0,0,0]
real(rp), parameter :: dim_energy(7)         = [2,1,-2,0,0,0,0]
real(rp), parameter :: dim_momentum(7)       = [1,1,-1,0,0,0,0]
real(rp), parameter :: dim_tesla(7)          = [0,1,-2,-1,0,0,0]

type(pmd_unit_struct), parameter :: unit_1          = pmd_unit_struct('', 1.0_rp, dim_1)
type(pmd_unit_struct), parameter :: unit_m          = pmd_unit_struct('m', 1.0_rp, dim_length)
type(pmd_unit_struct), parameter :: unit_kg         = pmd_unit_struct('kg', 1.0_rp, dim_mass)
type(pmd_unit_struct), parameter :: unit_sec        = pmd_unit_struct('sec', 1.0_rp, dim_time)
type(pmd_unit_struct), parameter :: unit_amp        = pmd_unit_struct('Amp', 1.0_rp, dim_current)
type(pmd_unit_struct), parameter :: unit_K          = pmd_unit_struct('K', 1.0_rp, dim_temperture)
type(pmd_unit_struct), parameter :: unit_mol        = pmd_unit_struct('mol', 1.0_rp, dim_mol)
type(pmd_unit_struct), parameter :: unit_cd         = pmd_unit_struct('cd', 1.0_rp, dim_luminous)
type(pmd_unit_struct), parameter :: unit_Coulomb    = pmd_unit_struct('Coulomb', 1.0_rp, dim_charge)
type(pmd_unit_struct), parameter :: unit_charge_num = pmd_unit_struct('charge #', 1.0_rp, dim_charge)
type(pmd_unit_struct), parameter :: unit_V_per_m    = pmd_unit_struct('V/m', 1.0_rp, dim_elec_potential)
type(pmd_unit_struct), parameter :: unit_c_light    = pmd_unit_struct('vel/c', c_light, dim_velocity)
type(pmd_unit_struct), parameter :: unit_m_per_s    = pmd_unit_struct('m/s', 1.0_rp, dim_velocity)
type(pmd_unit_struct), parameter :: unit_eV         = pmd_unit_struct('eV', e_charge, dim_energy)
type(pmd_unit_struct), parameter :: unit_eV_per_c   = pmd_unit_struct('eV/c', e_charge/c_light, dim_momentum)
type(pmd_unit_struct), parameter :: unit_Tesla      = pmd_unit_struct('Tesla', 1.0_rp, dim_tesla)

! 

interface pmd_write_real_to_dataset
  module procedure pmd_write_real_vector_to_dataset
  module procedure pmd_write_real_matrix_to_dataset
  module procedure pmd_write_real_tensor_to_dataset
end interface

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)
!
!-

subroutine pmd_write_real_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)

type (pmd_unit_struct) unit
real(rp) vector(:), v_max, v_min
integer err
integer(HID_T) :: root_id, v_size
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(vector)
v_min = minval(vector)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(vector), error) 
  return
endif

!

v_size = size(vector)
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, [v_size], vector, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_vector_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_matrix_to_dataset(root_id, dataset_name, bmad_name, unit, matrix, error)
!
!-

subroutine pmd_write_real_matrix_to_dataset(root_id, dataset_name, bmad_name, unit, matrix, error)

type (pmd_unit_struct) unit
real(rp) matrix(:,:), v_max, v_min
integer err
integer(HID_T) :: root_id, v_size(3)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(matrix)
v_min = minval(matrix)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(matrix), error) 
  return
endif

!

v_size = size(matrix)
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, v_size, matrix, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_matrix_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_tensor_to_dataset(root_id, dataset_name, bmad_name, unit, tensor, error)
!
!-

subroutine pmd_write_real_tensor_to_dataset(root_id, dataset_name, bmad_name, unit, tensor, error)

type (pmd_unit_struct) unit
real(rp) tensor(:,:,:), v_max, v_min
integer err
integer(HID_T) :: root_id, v_size(3)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(tensor)
v_min = minval(tensor)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(tensor), error) 
  return
endif

!

v_size = size(tensor)
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, v_size, tensor, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_tensor_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)
!
!-

subroutine pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)

type (pmd_unit_struct) unit
integer(hid_t) :: root_id, group_id
integer(size_t) n_dim
real(rp) value
integer v_size(:), err
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_double_f(root_id, dataset_name, 'value', [value], 1_size_t, err)
n_dim = size(v_size)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', v_size, n_dim, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_real_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)
!
!-

subroutine pmd_write_int_vector_to_dataset(root_id, dataset_name, bmad_name, unit, vector, error)

type (pmd_unit_struct) unit
integer vector(:), v_max, v_min
integer err
integer(HID_T) :: root_id, v_size
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(vector)
v_min = minval(vector)

if (v_max == v_min) then
  call pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, [size(vector)], error) 
  return
endif

!

v_size = size(vector)
call H5LTmake_dataset_int_f(root_id, dataset_name, 1, [v_size], vector, err)    

call H5LTset_attribute_int_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_int_vector_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)
!
!-

subroutine pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)

type (pmd_unit_struct) unit
integer(hid_t) :: root_id, group_id
integer(size_t) :: n_dim
integer value
integer err, v_size(:)
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_int_f(root_id, dataset_name, 'value', [value], 1_size_t, err)
n_dim = size(v_size)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', [v_size], n_dim, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_int_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)
!
!-

subroutine pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id
integer err
logical error
character(*) dataset_name, bmad_name

if (bmad_name /= '') call H5LTset_attribute_string_f(root_id, dataset_name, 'localName', bmad_name, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitSI', [unit%unitSI], 1_size_t, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitDimension', [unit%unitDimension], 7_size_t, err)
call H5LTset_attribute_string_f(root_id, dataset_name, 'unitSymbol', unit%unitSymbol, err)

end subroutine pmd_write_units_to_dataset 

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_int_dataset(root_id, name, value, error)
!
!-

subroutine pmd_read_int_dataset(root_id, name, value, error)

type (hdf5_info_struct) info
type (pmd_unit_struct) unit

real(rp) unit_si

integer(HID_T) :: root_id, obj_id
integer h5_err, value(:)

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_int_dataset'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.);  if (error) return
if (abs(unit_si - 1.0_rp) > 1d-6) then
  call out_io (s_error$, r_name, 'CONVERSION TO SI OF A VALUE OTHER THAN 1 DOES NOT MAKE SENSE.', &
                                  'FOR BUNCH PARAMETER: ' // name)
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (any(info%data_dim(2:) /= 0)) then
    call out_io (s_error$, r_name, 'DATA ARRAY IS NOT ONE-DIMENSIONAL! FOR DATA: ' // name)
    return
  endif

  call hdf5_read_dataset_int(root_id, name, value, err)

!

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_int(obj_id, 'value', value(1), error, .true.)
  value = value(1)
endif

!

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_int_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_dataset(root_id, name, conversion_factor, value, error)
!
!-

subroutine pmd_read_real_dataset(root_id, name, conversion_factor, value, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, value(:)
real(rp) unit_si

integer(HID_T) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_real_dataset'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.);  if (error) return

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (any(info%data_dim(2:) /= 0)) then
    call out_io (s_error$, r_name, 'DATA ARRAY IS NOT ONE-DIMENSIONAL! FOR DATA: ' // name)
    return
  endif

  call hdf5_read_dataset_real(root_id, name, value, err)

!

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', value(1), error, .true.)
  value = value(1)
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) value = value * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_tensor_dataset(root_id, name, conversion_factor, value, error)
!
!-

subroutine pmd_read_real_tensor_dataset(root_id, name, conversion_factor, value, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, value(:,:,:)
real(rp) unit_si, val

integer(HID_T) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_real_tensor_dataset'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.);  if (error) return

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (any(info%data_dim == 0)) then
    call out_io (s_error$, r_name, 'DATA ARRAY IS NOT THREE-DIMENSIONAL! FOR DATA: ' // name)
    return
  endif

  call hdf5_read_dataset_real(root_id, name, value, err)

!

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', val, error, .true.)
  value = val
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) value = value * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_tensor_dataset

end module
