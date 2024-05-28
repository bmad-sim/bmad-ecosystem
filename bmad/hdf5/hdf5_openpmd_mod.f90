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

character(6), parameter :: xyz_axislabels(3) = [character(8):: 'x', 'y', 'z']        
character(6), parameter :: rthetaz_axislabels(3) = [character(8):: 'r', 'theta', 'z']

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
! Subroutine pmd_init_compound_complex(complex_t)
!
! Routine to initialize complex data handling using a compound data type.
! Use complex_t in calls to read and write complex data.
!
! Note: Use pmd_kill_compound_complex to reclaim allocated memory.
!
! Output:
!   complex_t -- integer(hid_t): Complex compound data type identifier.
!-

subroutine pmd_init_compound_complex(complex_t)

integer(hid_t) :: complex_t, re_t
integer hdferr
integer(size_t) offset, re_size

! 

re_t = H5kind_to_type(rp, H5_REAL_KIND)
re_size = storage_size(1_rp) / 8

call H5Tcreate_f (H5T_COMPOUND_F, 2*re_size, complex_t, hdferr)
offset = 0
call H5Tinsert_f (complex_t, "r", offset, re_t, hdferr)
call H5Tinsert_f (complex_t, "i", re_size, re_t, hdferr)

end subroutine pmd_init_compound_complex

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_kill_compound_complex(complex_t)
!
! Routine to do memory cleanup for complex data handling using a compound data type.
! Note: Use pmd_init_compound_complex to first initalize complex_t 
!
! Input:
!   complex_t -- integer(hid_t): Complex compound data type identifier.
!-

subroutine pmd_kill_compound_complex(complex_t)

integer(hid_t) :: complex_t
integer hdferr

! 

call h5tclose_f(complex_t, hdferr)

end subroutine pmd_kill_compound_complex

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
! Subroutine pmd_write_complex_to_dataset_rank1 (root_id, dataset_name, complex_t, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 1.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   complex_t           -- integer(hid_t): Complex type ID obtained from pmd_init_compound_complex.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:)            -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_complex_to_dataset_rank1 (root_id, dataset_name, complex_t, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
complex(rp), target :: array(:), cc(size(array,1))
integer h5_err
integer(hid_t) :: root_id, dspace_id, z_id, complex_t
integer(hsize_t) dims(1)
character(*) dataset_name, bmad_name
logical err, error

! 

!if (all(array == array(1))) then
!  call pmd_write_complex_to_pseudo_dataset(root_id, dataset_name, complex_t, bmad_name, unit, array(1), shape(array), error)
!  return
!endif

! Need to use cc for temp storage since array argument may not be stored in contiguous memory.

dims = size(array)
call H5Screate_simple_f(1, dims, dspace_id, h5_err)  ! Create dataspace
call H5Dcreate_f(root_id, dataset_name, complex_t, dspace_id, z_id, h5_err)
cc = array
call H5dwrite_f(z_id, complex_t, c_loc(cc), h5_err)
call H5Dclose_f(z_id, h5_err)
call H5Sclose_f(dspace_id, h5_err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_complex_to_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_complex_to_dataset_rank2 (root_id, dataset_name, complex_t, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   complex_t           -- integer(hid_t): Complex type ID obtained from pmd_init_compound_complex.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:)          -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_complex_to_dataset_rank2 (root_id, dataset_name, complex_t, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
complex(rp), target :: array(:,:), cc(size(array,1),size(array,2))
integer h5_err
integer(hid_t) :: root_id, dspace_id, z_id, complex_t
integer(hsize_t) dims(2)
character(*) dataset_name, bmad_name
logical error, err

! 

!if (all(array == array(1,1))) then
!  call pmd_write_complex_to_pseudo_dataset(root_id, dataset_name, complex_t, bmad_name, unit, array(1,1), shape(array), error)
!  return
!endif

! Need to use cc for temp storage since array argument may not be stored in contiguous memory.

dims = shape(array)
call H5Screate_simple_f(2, dims, dspace_id, h5_err)  ! Create dataspace
call H5Dcreate_f(root_id, dataset_name, complex_t, dspace_id, z_id, h5_err)
cc = array
call H5dwrite_f(z_id, complex_t, c_loc(cc), h5_err)
call H5Dclose_f(z_id, h5_err)
call H5Sclose_f(dspace_id, h5_err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_complex_to_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_complex_to_dataset_rank3 (root_id, dataset_name, complex_t, bmad_name, unit, array, error)
!
! Routine to write an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   complex_t           -- integer(hid_t): Complex type ID obtained from pmd_init_compound_complex.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   array(:,:,:)        -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_complex_to_dataset_rank3 (root_id, dataset_name, complex_t, bmad_name, unit, array, error)

type (pmd_unit_struct) unit
complex(rp), target :: array(:,:,:), cc(size(array,1),size(array,2),size(array,3))
integer h5_err
integer(hid_t) :: root_id, dspace_id, z_id, complex_t
integer(hsize_t) dims(3)
character(*) dataset_name, bmad_name
logical err, error

! 

!if (all(array == array(1,1,1))) then
!  call pmd_write_complex_to_pseudo_dataset(root_id, dataset_name, complex_t, bmad_name, unit, array(1,1,1), shape(array), error)
!  return
!endif

! Need to use cc for temp storage since array argument may not be stored in contiguous memory.

dims = shape(array)
call H5Screate_simple_f(3, dims, dspace_id, h5_err)  ! Create dataspace
call H5Dcreate_f(root_id, dataset_name, complex_t, dspace_id, z_id, h5_err)
cc = array
call H5dwrite_f(z_id, complex_t, c_loc(cc), h5_err)
call H5Dclose_f(z_id, h5_err)
call H5Sclose_f(dspace_id, h5_err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_complex_to_dataset_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_write_complex_to_pseudo_dataset(root_id, dataset_name, complex_t, bmad_name, unit, value, v_shape, error)
!
! Routine to write an openpmd formatted dataset of an array were all the data values are the same.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   complex_t           -- integer(hid_t): Complex type ID obtained from pmd_init_compound_complex.
!   bmad_name           -- character(*): Equivalent Bmad name.
!   unit                -- pmd_unit_struct: Data units.
!   value               -- complex(rp): Data value.
!   v_shape(:)          -- integer: Shape of the data array.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_write_complex_to_pseudo_dataset(root_id, dataset_name, complex_t, bmad_name, unit, value, v_shape, error)

type (pmd_unit_struct) unit
integer(hid_t) :: root_id, group_id, complex_t
integer(size_t) :: n_dim
complex(rp) value
integer err, v_shape(:)
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call my_H5LTset_attribute_complex(root_id, dataset_name, complex_t, 'value', value, error) 

n_dim = size(v_shape)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', v_shape, n_dim, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_complex_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine my_H5LTset_attribute_complex(root_id, dataset_name, complex_t, attrib_name, value, err) 
!
! Routine to set a complex attribute.
! Modified from H5LTset_attribute_double_f
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   dataset_name        -- character(*): Name of the dataset.
!   complex_t           -- integer(hid_t): Complex type ID obtained from pmd_init_compound_complex.
!   attrib_name         -- character(*): Atrribute name.
!   value               -- complex(rp): Data value.
!
! Output:
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine my_H5LTset_attribute_complex(root_id, dataset_name, complex_t, attrib_name, value, error) 

integer(hid_t) root_id, complex_t, array_size, dataspace_id, attrib_id
integer(size_t) name_len, attrib_len, buf_size
type(c_ptr) f_ptr
character(*) dataset_name, attrib_name
complex(rp), target :: value
integer err
logical error

!

error = .true.

call H5Screate_f(H5S_SCALAR_F, dataspace_id, err)
if (err /= 0) return

call H5Acreate_f(root_id, attrib_name, complex_t, dataspace_id, attrib_id, err, H5P_DEFAULT_F, H5P_DEFAULT_F)
if (err /= 0) return

f_ptr = c_loc(value)
CALL H5Awrite_f(attrib_id, complex_t, f_ptr, err)
if (err /= 0) return

error = .false.

end subroutine my_H5LTset_attribute_complex

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
integer h5_err, n
logical error
character(*) dataset_name, bmad_name

if (bmad_name /= '') call H5LTset_attribute_string_f(root_id, dataset_name, 'localName', bmad_name, h5_err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitSI', [unit%unitSI], 1_size_t, h5_err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitDimension', [unit%unitDimension], 7_size_t, h5_err)
call H5LTset_attribute_string_f(root_id, dataset_name, 'unitSymbol', unit%unitSymbol, h5_err)

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
  call hdf5_read_attribute_int(obj_id, 'value', c_val, error, .true., name)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_int_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_int_dataset_rank2 (root_id, name, conversion_factor, axislabels, array, error)
!
! Routine to read an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   axislabels(:)       -- character(*): Axis labels.
!   array(:,:)          -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:)          -- integer: Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_int_dataset_rank2 (root_id, name, conversion_factor, axislabels, array, error)

type (hdf5_info_struct) info
type (pmd_unit_struct) unit

real(rp) unit_si, conversion_factor

integer(hid_t) :: root_id, obj_id
integer h5_err, array(:,:), c_val

logical error, err

character(*) name
character(*) axislabels(:)
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

  call hdf5_read_dataset_int(root_id, name, axislabels, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_int(obj_id, 'value', c_val, error, .true., name)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_int_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_int_dataset_rank3 (root_id, name, conversion_factor, axislabels, array, error)
!
! Routine to read an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   axislabels(:)       -- character(*): Axis labels.
!   array(:,:,:)        -- integer: Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:,:)        -- integer: Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_int_dataset_rank3 (root_id, name, conversion_factor, axislabels, array, error)

type (hdf5_info_struct) info
type (pmd_unit_struct) unit

real(rp) unit_si, conversion_factor

integer(hid_t) :: root_id, obj_id
integer h5_err, array(:,:,:), c_val

logical error, err

character(*) name
character(*) axislabels(:)
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

  call hdf5_read_dataset_int(root_id, name, axislabels, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_int(obj_id, 'value', c_val, error, .true., name)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

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
  call hdf5_read_attribute_real(obj_id, 'value', c_val, error, .true., name)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_dataset_rank2 (root_id, name, conversion_factor, axislabels, array, error)
!
! Routine to read an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   axislabels(:)       -- character(*): Axis labels.
!   array(:,:)          -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:)          -- real(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_real_dataset_rank2 (root_id, name, conversion_factor, axislabels, array, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, array(:,:), c_val
real(rp) unit_si

integer(hid_t) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*) axislabels(:)
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

  call hdf5_read_dataset_real(root_id, name, axislabels, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', c_val, error, .true., name)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_real_dataset_rank3 (root_id, name, conversion_factor, axislabels, array, error)
!
! Routine to read an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   axislabels(:)       -- character(*): Axis labels.
!   array(:,:,:)        -- real(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:,:)        -- real(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_real_dataset_rank3 (root_id, name, conversion_factor, axislabels, array, error)

type (hdf5_info_struct) info

real(rp) conversion_factor, array(:,:,:), c_val
real(rp) unit_si

integer(hid_t) :: root_id, obj_id
integer h5_err

logical error, err

character(*) name
character(*), optional :: axislabels(:)
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

  call hdf5_read_dataset_real(root_id, name, axislabels, array, err)

else  ! Must be a "constant record component" as defined by the openPMD standard
  call hdf5_read_attribute_real(obj_id, 'value', c_val, error, .true., name)
  array = c_val
endif

!

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)
call hdf5_close_object(obj_id, info)

end subroutine pmd_read_real_dataset_rank3

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_complex_dataset_rank1 (root_id, name, complex_t, conversion_factor, array, error)
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

subroutine pmd_read_complex_dataset_rank1 (root_id, name, complex_t, conversion_factor, array, error)

type (hdf5_info_struct) info

complex(rp) array(:)
complex(rp), target :: cc(size(array,1))
real(rp), allocatable :: re(:), im(:)
real(rp) conversion_factor, unit_si

integer(hid_t) :: root_id, z_id, complex_t
integer h5_err
type(c_ptr) :: f_ptr

logical error, err

character(*) name
character(*), parameter :: r_name = 'pmd_read_complex_dataset_rank1'

! Non-compound data means old format

error = .true.

info = hdf5_object_info(root_id, name, error, .true.)
z_id = hdf5_open_object(root_id, name, info, error, .true.)

if (info%data_class_type /= H5T_COMPOUND_F) then
  allocate (re(size(array,1)), im(size(array,1)))

  call h5gopen_f(root_id, name, z_id, h5_err);                     if (h5_err == -1) return
  call pmd_read_real_dataset (z_id, 'r', conversion_factor, re, err);      if (err) return
  call pmd_read_real_dataset (z_id, 'i', conversion_factor, im, err);      if (err) return
  call h5gclose_f(z_id, h5_err)

  array = cmplx(re, im, rp)
  error = .false.
  return
endif

! Need to use cc for temp storage since array argument may not be stored in contiguous memory.

if (info%data_dim(1) /= size(array)) then
  call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT OF THE CORRECT SIZE! FOR DATA: ' // name)
  return
endif

f_ptr = c_loc(cc)
call H5Dread_f(z_id, complex_t, f_ptr, h5_err)
array = cc

call hdf5_read_attribute_real(z_id, 'unitSI', unit_si, error, .true.)
if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

call hdf5_close_object(z_id, info)

end subroutine pmd_read_complex_dataset_rank1

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_complex_dataset_rank2 (root_id, name, complex_t, conversion_factor, axislabels, array, error)
!
! Routine to read an openpmd formatted dataset of rank 2.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   axislabels(:)       -- character(*): Axis labels.
!   array(:,:)          -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:)          -- complex(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_complex_dataset_rank2 (root_id, name, complex_t, conversion_factor, axislabels, array, error)

type (hdf5_info_struct) info

complex(rp) array(:,:)
complex(rp), target, allocatable :: cc(:,:)
real(rp), allocatable :: re(:,:), im(:,:)
real(rp) conversion_factor, unit_si

integer(hid_t) :: root_id, z_id, complex_t
integer h5_err, i
type(c_ptr) :: f_ptr

logical error, err

character(*) name
character(*) axislabels(:)
character(1) d_ord
character(*), parameter :: r_name = 'pmd_read_complex_dataset_rank2'

! Non-compound data means old format

error = .true.

info = hdf5_object_info(root_id, name, error, .true.)

if (info%data_class_type /= H5T_COMPOUND_F) then
  allocate (re(size(array,1), size(array,2)), im(size(array,1), size(array,2)))

  call h5gopen_f(root_id, name, z_id, h5_err); if (h5_err == -1) return
  call pmd_read_real_dataset (z_id, 'r', conversion_factor, axislabels, re, err); if (err) return
  call pmd_read_real_dataset (z_id, 'i', conversion_factor, axislabels, im, err); if (err) return
  call h5gclose_f(z_id, h5_err)

  array = cmplx(re, im, rp)
  error = .false.
  return
endif

! Need to use cc for temp storage since array argument may not be stored in contiguous memory.

call hdf5_read_dataorder(root_id, name, axislabels, d_ord)
if (d_ord == 'C') then
  allocate (cc(size(array,2), size(array,1)))
else
  allocate (cc(size(array,1), size(array,2)))
endif

if (any(info%data_dim(1:2) /= shape(cc))) then
  call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT OF THE CORRECT SIZE! FOR DATA: ' // name)
  return
endif

z_id = hdf5_open_object(root_id, name, info, error, .true.)
f_ptr = c_loc(cc)
call H5Dread_f(z_id, complex_t, f_ptr, h5_err)
call hdf5_read_attribute_real(z_id, 'unitSI', unit_si, error, .true.)
call hdf5_close_object(z_id, info)

if (d_ord == 'C') then
  do i = 1, size(array,1)
    array(i,:) = cc(:,i)
  enddo
else
  array = cc
endif

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

end subroutine pmd_read_complex_dataset_rank2

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine pmd_read_complex_dataset_rank3 (root_id, name, complex_t, conversion_factor, axislabels, array, error)
!
! Routine to read an openpmd formatted dataset of rank 3.
!
! Input:
!   root_id             -- integer(hid_t): Root group containing the dataset.
!   name                -- character(*): Name of the dataset.
!   conversion_factor   -- real(rp): Conversion factor from SI units to Bmad units.
!   axislabels(:)       -- character(*): Axis labels.
!   array(:,:,:)        -- complex(rp): Array to hold the data. Must be of the correct size.
!
! Output:
!   array(:,:,:)        -- complex(rp): Array with data.
!   error               -- logical: Set true if there is an error. False otherwise.
!-

subroutine pmd_read_complex_dataset_rank3 (root_id, name, complex_t, conversion_factor, axislabels, array, error)

type (hdf5_info_struct) info

complex(rp) array(:,:,:)
complex(rp), target, allocatable :: cc(:,:,:)
real(rp), allocatable :: re(:,:,:), im(:,:,:)
real(rp) conversion_factor, unit_si

integer(hid_t) :: root_id, z_id, complex_t
integer h5_err, i, j
type(c_ptr) :: f_ptr

logical error, err

character(*) name
character(*) axislabels(:)
character(1) d_ord
character(*), parameter :: r_name = 'pmd_read_complex_dataset_rank3'

! Non-compound data means old format

error = .true.

info = hdf5_object_info(root_id, name, error, .true.)

if (info%data_class_type /= H5T_COMPOUND_F) then
  allocate (re(size(array,1), size(array,2), size(array,3)), im(size(array,1), size(array,2), size(array,3)))

  call h5gopen_f(root_id, name, z_id, h5_err);                     if (h5_err == -1) return
  call pmd_read_real_dataset (z_id, 'r', conversion_factor, axislabels, re, err);      if (err) return
  call pmd_read_real_dataset (z_id, 'i', conversion_factor, axislabels, im, err);      if (err) return
  call h5gclose_f(z_id, h5_err)

  array = cmplx(re, im, rp)
  error = .false.
  return
endif

! Need to use cc for temp storage since array argument may not be stored in contiguous memory.

call hdf5_read_dataorder(root_id, name, axislabels, d_ord)
if (d_ord == 'C') then
  allocate (cc(size(array,3), size(array,2), size(array,1)))
else
  allocate (cc(size(array,1), size(array,2), size(array,3)))
endif

if (any(info%data_dim(1:3) /= shape(cc))) then
  call out_io (s_error$, r_name, 'STORED DATA ARRAY IS NOT OF THE CORRECT SIZE! FOR DATA: ' // name)
  return
endif

z_id = hdf5_open_object(root_id, name, info, error, .true.)
f_ptr = c_loc(cc)
call H5Dread_f(z_id, complex_t, f_ptr, h5_err)
call hdf5_read_attribute_real(z_id, 'unitSI', unit_si, error, .true.)
call hdf5_close_object(z_id, info)

if (d_ord == 'C') then
  do i = 1, size(array,1);  do j = 1, size(array,2)
    array(i,j,:) = cc(:,j,i)
  enddo;  enddo
else
  array = cc
endif

if (abs(unit_si - conversion_factor) > 1d-15 * conversion_factor) array = array * (conversion_factor / unit_si)

end subroutine pmd_read_complex_dataset_rank3

end module
