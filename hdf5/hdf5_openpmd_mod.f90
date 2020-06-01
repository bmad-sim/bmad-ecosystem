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
real(rp), parameter :: dim_electric_field(7) = [1,1,-3,-1,0,0,0]
real(rp), parameter :: dim_velocity(7)       = [1,0,-1,0,0,0,0]
real(rp), parameter :: dim_energy(7)         = [2,1,-2,0,0,0,0]
real(rp), parameter :: dim_momentum(7)       = [1,1,-1,0,0,0,0]
real(rp), parameter :: dim_tesla(7)          = [0,1,-2,-1,0,0,0]
real(rp), parameter :: dim_hbar(7)           = [2,1,-1,0,0,0,0]

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
type(pmd_unit_struct), parameter :: unit_V_per_m    = pmd_unit_struct('V/m', 1.0_rp, dim_electric_field)
type(pmd_unit_struct), parameter :: unit_c_light    = pmd_unit_struct('vel/c', c_light, dim_velocity)
type(pmd_unit_struct), parameter :: unit_m_per_s    = pmd_unit_struct('m/s', 1.0_rp, dim_velocity)
type(pmd_unit_struct), parameter :: unit_eV         = pmd_unit_struct('eV', e_charge, dim_energy)
type(pmd_unit_struct), parameter :: unit_eV_per_c   = pmd_unit_struct('eV/c', e_charge/c_light, dim_momentum)
type(pmd_unit_struct), parameter :: unit_Tesla      = pmd_unit_struct('Tesla', 1.0_rp, dim_tesla)
type(pmd_unit_struct), parameter :: unit_hbar       = pmd_unit_struct('hbar', e_charge * h_bar_planck, dim_hbar)

! 

interface pmd_write_int_to_dataset
  module procedure pmd_write_int_to_dataset_rank1
  module procedure pmd_write_int_to_dataset_rank2
  module procedure pmd_write_int_to_dataset_rank3
end interface

interface pmd_write_real_to_dataset
  module procedure pmd_write_real_to_dataset_rank1
  module procedure pmd_write_real_to_dataset_rank2
  module procedure pmd_write_real_to_dataset_rank3
end interface

interface pmd_write_complex_to_dataset
  module procedure pmd_write_complex_to_dataset_rank1
  module procedure pmd_write_complex_to_dataset_rank2
  module procedure pmd_write_complex_to_dataset_rank3
end interface


interface pmd_read_int_dataset
  module procedure pmd_read_int_dataset_rank1
  module procedure pmd_read_int_dataset_rank2
  module procedure pmd_read_int_dataset_rank3
end interface

interface pmd_read_real_dataset
  module procedure pmd_read_real_dataset_rank1
  module procedure pmd_read_real_dataset_rank2
  module procedure pmd_read_real_dataset_rank3
end interface

interface pmd_read_complex_dataset
  module procedure pmd_read_complex_dataset_rank1
  module procedure pmd_read_complex_dataset_rank2
  module procedure pmd_read_complex_dataset_rank3
end interface


contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_to_dataset_rank1 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 1.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:)            -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_int_to_dataset_rank1 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
integer array(:), v_max, v_min
integer err
integer(hid_t) :: root_id, v_shape(1)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(array)
v_min = minval(array)

if (v_max == v_min) then
  call pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(array), error) 
  return
endif

!

v_shape = shape(array)
call H5LTmake_dataset_int_f(root_id, dataset_name, 1, v_shape, array, err)    

call H5LTset_attribute_int_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_int_to_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_to_dataset_rank2 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:)          -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_int_to_dataset_rank2 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
integer array(:,:), v_max, v_min
integer err
integer(hid_t) :: root_id, v_shape(2)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(array)
v_min = minval(array)

if (v_max == v_min) then
  call pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(array), error) 
  return
endif

!

v_shape = shape(array)
call H5LTmake_dataset_int_f(root_id, dataset_name, 2, v_shape, array, err)    

call H5LTset_attribute_int_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_int_to_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_to_dataset_rank3 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:,:)        -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_int_to_dataset_rank3 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
integer array(:,:,:), v_max, v_min
integer err
integer(hid_t) :: root_id, v_shape(3)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(array)
v_min = minval(array)

if (v_max == v_min) then
  call pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(array), error) 
  return
endif

!

v_shape = shape(array)
call H5LTmake_dataset_int_f(root_id, dataset_name, 3, v_shape, array, err)    

call H5LTset_attribute_int_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_int_to_dataset_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_shape, error)
!
! Routine to write an openpmd formatted dataset of an array were all the data values are the same.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   value               -- integer: Data value.
!   v_shape(:)          -- integer: Shape of the data array.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_shape, error)

type (pmd_unit_struct) unit
integer(hid_t) :: root_id, group_id
integer(size_t) :: n_dim
integer value
integer err, v_shape(:)
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_int_f(root_id, dataset_name, 'value', [value], 1_size_t, err)
n_dim = size(v_shape)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', v_shape, n_dim, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_int_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_to_dataset_rank1 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 1.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:)            -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_real_to_dataset_rank1 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
real(rp) array(:), v_max, v_min
integer err
integer(hid_t) :: root_id, v_shape(1)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(array)
v_min = minval(array)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(array), error) 
  return
endif

!

v_shape = shape(array)
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, v_shape, array, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_to_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_to_dataset_rank2 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:)          -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_real_to_dataset_rank2 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
real(rp) array(:,:), v_max, v_min
integer err
integer(hid_t) :: root_id, v_shape(2)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(array)
v_min = minval(array)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset (root_id, dataset_name, bmad_name, unit, v_max, shape(array), error) 
  return
endif

!

v_shape = shape(array)
call H5LTmake_dataset_double_f(root_id, dataset_name, 2, v_shape, array, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_to_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_to_dataset_rank3 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:,:)        -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_real_to_dataset_rank3 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
real(rp) array(:,:,:), v_max, v_min
integer err
integer(hid_t) :: root_id, v_shape(3)
character(*) dataset_name, bmad_name
logical error

!

v_max = maxval(array)
v_min = minval(array)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, shape(array), error) 
  return
endif

!

v_shape = shape(array)
call H5LTmake_dataset_double_f(root_id, dataset_name, 3, v_shape, array, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], 1_size_t, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], 1_size_t, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_to_dataset_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_real_to_pseudo_dataset (root_id, dataset_name, bmad_name, unit, value, v_shape, error)
!
! Routine to write an openpmd formatted dataset of an array were all the data values are the same.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   value               -- real(rp): Data value.
!   v_shape(:)          -- integer: Shape of the data array.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_real_to_pseudo_dataset (root_id, dataset_name, bmad_name, unit, value, v_shape, error)

type (pmd_unit_struct) unit
integer(hid_t) :: root_id, group_id
integer(size_t) n_dim
real(rp) value
integer v_shape(:), err
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_double_f(root_id, dataset_name, 'value', [value], 1_size_t, err)
n_dim = size(v_shape)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', v_shape, n_dim, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_real_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_complex_to_dataset_rank1 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 1.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:)            -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_complex_to_dataset_rank1 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
complex(rp) array(:)
integer h5_err
integer(hid_t) :: root_id, z_id
character(*) dataset_name, bmad_name
logical err, error

!

call h5gcreate_f(root_id, dataset_name, z_id, h5_err)
call pmd_write_real_to_dataset (z_id, 'r', 'real', unit, real(array), err)
call pmd_write_real_to_dataset (z_id, 'i', 'imaginary', unit, aimag(array), err)
call h5gclose_f(z_id, h5_err)

end subroutine pmd_write_complex_to_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_complex_to_dataset_rank2 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:)          -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_complex_to_dataset_rank2 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
complex(rp) array(:,:)
integer h5_err
integer(hid_t) :: root_id, z_id
character(*) dataset_name, bmad_name
logical error, err

!

call h5gcreate_f(root_id, dataset_name, z_id, h5_err)
call pmd_write_real_to_dataset (z_id, 'r', 'real', unit, real(array), err)
call pmd_write_real_to_dataset (z_id, 'i', 'imaginary', unit, aimag(array), err)
call h5gclose_f(z_id, h5_err)

end subroutine pmd_write_complex_to_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_complex_to_dataset_rank3 (root_id, dataset_name, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:,:)        -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_complex_to_dataset_rank3 (root_id, dataset_name, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
complex(rp) array(:,:,:)
integer h5_err
integer(hid_t) :: root_id, z_id
character(*) dataset_name, bmad_name
logical error, err

!

call h5gcreate_f(root_id, dataset_name, z_id, h5_err)
call pmd_write_real_to_dataset (z_id, 'r', 'real', unit, real(array), err)
call pmd_write_real_to_dataset (z_id, 'i', 'imaginary', unit, aimag(array), err)
call h5gclose_f(z_id, h5_err)

end subroutine pmd_write_complex_to_dataset_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)
!
! Routine to write the data units and conversion factors to a dataset.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

type (pmd_unit_struct) unit
integer(hid_t) :: root_id
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
! Subroutine pmd_read_int_dataset_rank1 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 1.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:)            -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:)            -- integer: Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_int_dataset_rank1 (root_id, name, conversion_factor, array, error)

type (hdf5_info_struct) info
type (pmd_unit_struct) unit

real(rp) unit_si, conversion_factor

integer(hid_t) :: root_id, obj_id
integer h5_err, array(:), c_val

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_int_dataset_rank1'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.)
if (error) then
  call hdf5_close_object(obj_id, info)
  return
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (info%data_dim(2) /= 0) then
    call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT ONE-DIMENSIONAL! FOR DATA: ' // name)
    call hdf5_close_object(obj_id, info)
    return
  endif

  call hdf5_read_dataset_int(root_id, name, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_int(obj_id, 'value', c_val, error, .true.)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_int_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_int_dataset_rank2 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:,:)          -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:)          -- integer: Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_int_dataset_rank2 (root_id, name, conversion_factor, array, error)

type (hdf5_info_struct) info
type (pmd_unit_struct) unit

real(rp) unit_si, conversion_factor

integer(hid_t) :: root_id, obj_id
integer h5_err, array(:,:), c_val

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_int_dataset_rank2'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.)
if (error) then
  call hdf5_close_object(obj_id, info)
  return
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (info%data_dim(2) == 0 .or. info%data_dim(3) /= 0) then
    call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT TWO-DIMENSIONAL! FOR DATA: ' // name)
    call hdf5_close_object(obj_id, info)
    return
  endif

  call hdf5_read_dataset_int(root_id, name, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_int(obj_id, 'value', c_val, error, .true.)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_int_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_int_dataset_rank3 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:,:,:)        -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:,:)        -- integer: Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_int_dataset_rank3 (root_id, name, conversion_factor, array, error)

type (hdf5_info_struct) info
type (pmd_unit_struct) unit

real(rp) unit_si, conversion_factor

integer(hid_t) :: root_id, obj_id
integer h5_err, array(:,:,:), c_val

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_int_dataset_rank3'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.)
if (error) then
  call hdf5_close_object(obj_id, info)
  return
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (info%data_dim(3) == 0) then
    call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT THREE-DIMENSIONAL! FOR DATA: ' // name)
    call hdf5_close_object(obj_id, info)
    return
  endif

  call hdf5_read_dataset_int(root_id, name, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_int(obj_id, 'value', c_val, error, .true.)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_int_dataset_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_dataset_rank1 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 1.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:)            -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:)            -- real(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_real_dataset_rank1 (root_id, name, conversion_factor, array, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, array(:), c_val
real(rp) unit_si

integer(hid_t) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_real_dataset_rank1'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.)
if (error) then
  call out_io (s_error$, r_name, 'THIS ERROR GENERATED WHILE PARSING: ' // name)
  call hdf5_close_object(obj_id, info)
  return
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (info%data_dim(2) /= 0) then
    call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT ONE-DIMENSIONAL! FOR DATA: ' // name)
    call hdf5_close_object(obj_id, info)
    return
  endif

  call hdf5_read_dataset_real(root_id, name, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', c_val, error, .true.)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_dataset_rank2 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:,:)          -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:)          -- real(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_real_dataset_rank2 (root_id, name, conversion_factor, array, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, array(:,:), c_val
real(rp) unit_si

integer(hid_t) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_real_dataset_rank2'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.)
if (error) then
  call hdf5_close_object(obj_id, info)
  return
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (info%data_dim(2) == 0 .or. info%data_dim(3) /= 0) then
    call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT TWO-DIMENSIONAL! FOR DATA: ' // name)
    call hdf5_close_object(obj_id, info)
    return
  endif

  call hdf5_read_dataset_real(root_id, name, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', c_val, error, .true.)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_dataset_rank3 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:,:,:)        -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:,:)        -- real(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_real_dataset_rank3 (root_id, name, conversion_factor, array, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, array(:,:,:), c_val
real(rp) unit_si

integer(hid_t) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_real_dataset_rank3'

!

info = hdf5_object_info(root_id, name, error, .true.)
obj_id = hdf5_open_object(root_id, name, info, error, .true.)

call hdf5_read_attribute_real(obj_id, 'unitSI', unit_si, error, .true.)
if (error) then
  call hdf5_close_object(obj_id, info)
  return
endif

!

if (info%element_type == H5O_TYPE_DATASET_F) then
  if (info%data_dim(3) == 0) then
    call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT THREE-DIMENSIONAL! FOR DATA: ' // name)
    call hdf5_close_object(obj_id, info)
    return
  endif

  call hdf5_read_dataset_real(root_id, name, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', c_val, error, .true.)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1e-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_complex_dataset_rank1 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 1.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:)            -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:)            -- complex(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_complex_dataset_rank1 (root_id, name, conversion_factor, array, error)

complex(rp) array(:)
real(rp), allocatable :: re(:), im(:)
real(rp) conversion_factor

integer(hid_t) :: root_id, z_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_complex_dataset_rank1'

!

error = .true.

allocate (re(size(array)), im(size(array)))

call h5gopen_f(root_id, name, z_id, h5_err);                     if (h5_err == -1) return
call pmd_read_real_dataset (z_id, 'r', conversion_factor, re, err);      if (err) return
call pmd_read_real_dataset (z_id, 'i', conversion_factor, im, err);      if (err) return
call h5gclose_f(z_id, h5_err)

array = cmplx(re, im, rp)

error = .false.

end subroutine pmd_read_complex_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_complex_dataset_rank2 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:,:)          -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:)          -- complex(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_complex_dataset_rank2 (root_id, name, conversion_factor, array, error)

complex(rp) array(:,:)
real(rp), allocatable :: re(:,:), im(:,:)
real(rp) conversion_factor

integer(hid_t) :: root_id, z_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_complex_dataset_rank2'

!

error = .true.

allocate (re(size(array,1), size(array,2)), im(size(array,1), size(array,2)))

call h5gopen_f(root_id, name, z_id, h5_err);                     if (h5_err == -1) return
call pmd_read_real_dataset (z_id, 'r', conversion_factor, re, err);      if (err) return
call pmd_read_real_dataset (z_id, 'i', conversion_factor, im, err);      if (err) return
call h5gclose_f(z_id, h5_err)

array = cmplx(re, im, rp)

error = .false.

end subroutine pmd_read_complex_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_complex_dataset_rank3 (root_id, name, conversion_factor, array, error)
!
! Routine to read an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   array(:,:,:)        -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:,:)        -- complex(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_complex_dataset_rank3 (root_id, name, conversion_factor, array, error)

complex(rp) array(:,:,:)
real(rp), allocatable :: re(:,:,:), im(:,:,:)
real(rp) conversion_factor

integer(hid_t) :: root_id, z_id
integer h5_err

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_complex_dataset_rank3'

!

error = .true.

allocate (re(size(array,1), size(array,2), size(array,3)), im(size(array,1), size(array,2), size(array,3)))

call h5gopen_f(root_id, name, z_id, h5_err);                     if (h5_err == -1) return
call pmd_read_real_dataset (z_id, 'r', conversion_factor, re, err);      if (err) return
call pmd_read_real_dataset (z_id, 'i', conversion_factor, im, err);      if (err) return
call h5gclose_f(z_id, h5_err)

array = cmplx(re, im, rp)

error = .false.

end subroutine pmd_read_complex_dataset_rank3


end module
