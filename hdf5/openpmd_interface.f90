module openpmd_interface

use h5lt
use hdf5
use iso_fortran_env
use sim_utils

implicit none

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

type(pmd_unit_struct), parameter :: unit_1        = pmd_unit_struct('', 1.0_rp, dim_1)
type(pmd_unit_struct), parameter :: unit_m        = pmd_unit_struct('m', 1.0_rp, dim_length)
type(pmd_unit_struct), parameter :: unit_kg       = pmd_unit_struct('kg', 1.0_rp, dim_mass)
type(pmd_unit_struct), parameter :: unit_sec      = pmd_unit_struct('sec', 1.0_rp, dim_time)
type(pmd_unit_struct), parameter :: unit_amp      = pmd_unit_struct('Amp', 1.0_rp, dim_current)
type(pmd_unit_struct), parameter :: unit_K        = pmd_unit_struct('K', 1.0_rp, dim_temperture)
type(pmd_unit_struct), parameter :: unit_mol      = pmd_unit_struct('mol', 1.0_rp, dim_mol)
type(pmd_unit_struct), parameter :: unit_cd       = pmd_unit_struct('cd', 1.0_rp, dim_luminous)
type(pmd_unit_struct), parameter :: unit_Coulomb  = pmd_unit_struct('Coulomb', 1.0_rp, dim_charge)
type(pmd_unit_struct), parameter :: unit_V_per_m  = pmd_unit_struct('V/m', 1.0_rp, dim_elec_potential)
type(pmd_unit_struct), parameter :: unit_c_light  = pmd_unit_struct('vel/c', c_light, dim_velocity)
type(pmd_unit_struct), parameter :: unit_m_per_s  = pmd_unit_struct('m/s', 1.0_rp, dim_velocity)
type(pmd_unit_struct), parameter :: unit_eV       = pmd_unit_struct('eV', e_charge, dim_energy)
type(pmd_unit_struct), parameter :: unit_eV_per_c = pmd_unit_struct('eV/c', e_charge/c_light, dim_momentum)

! Header information

type pmd_header_struct
  character(:), allocatable :: openPMD
  character(:), allocatable :: openPMDextension
  character(:), allocatable :: basePath
  character(:), allocatable :: particlesPath
  character(:), allocatable :: author          != 'anonymous'
  character(:), allocatable :: software        != 'Bmad'
  character(:), allocatable :: softwareVersion ! = 'Revision XXX'
  character(:), allocatable :: date            ! = ''
  character(:), allocatable :: latticeFile
  character(:), allocatable :: latticeName
end type

! Misc

integer(hsize_t), parameter :: h5_size_1 = 1
integer(hsize_t), parameter :: h5_size_7 = 7

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_string_attrib(id, name, string)
integer(HID_T) :: id
character(*) :: name, string
integer error
!
call H5LTset_attribute_string_f(id, '.', name, trim(string), error)
end subroutine pmd_write_string_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_integer_attrib(id, name, i)
integer(HID_T) :: id
character(*) :: name
integer :: i
integer error
!
call H5LTset_attribute_int_f(id, '.', name, [i], h5_size_1, error)
end subroutine pmd_write_integer_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_real_attrib(id, name, x)
integer(HID_T) :: id
character(*) :: name
real(rp) :: x
integer error
call H5LTset_attribute_double_f(id, '.', name, [x], h5_size_1, error)
end subroutine pmd_write_real_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

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
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, size(vector), error) 
  return
endif

!

v_size = size(vector)
call H5LTmake_dataset_double_f(root_id, dataset_name, 1, [v_size], vector, err)    

call H5LTset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], h5_size_1, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], h5_size_1, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_real_vector_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_real_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id, group_id
real(rp) value
integer v_size, err
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_double_f(root_id, dataset_name, 'value', [value], h5_size_1, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', [v_size], h5_size_1, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_real_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

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
  call pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, v_max, size(vector), error) 
  return
endif

!

v_size = size(vector)
call H5LTmake_dataset_int_f(root_id, dataset_name, 1, [v_size], vector, err)    

call H5LTset_attribute_int_f(root_id, dataset_name, 'minValue', [v_min], h5_size_1, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'maxValue', [v_max], h5_size_1, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

end subroutine pmd_write_int_vector_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_int_to_pseudo_dataset(root_id, dataset_name, bmad_name, unit, value, v_size, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id, group_id
integer value
integer err, v_size
character(*) dataset_name, bmad_name
logical error

!

call h5gcreate_f(root_id, dataset_name, group_id, err)

call H5LTset_attribute_int_f(root_id, dataset_name, 'value', [value], h5_size_1, err)
call H5LTset_attribute_int_f(root_id, dataset_name, 'shape', [v_size], h5_size_1, err)

call pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

call h5gclose_f(group_id, err)

end subroutine pmd_write_int_to_pseudo_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_units_to_dataset (root_id, dataset_name, bmad_name, unit, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id
integer err
logical error
character(*) dataset_name, bmad_name

if (bmad_name /= '') call H5LTset_attribute_string_f(root_id, dataset_name, 'localName', bmad_name, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitSI', [unit%unitSI], h5_size_1, err)
call H5LTset_attribute_double_f(root_id, dataset_name, 'unitDimension', [unit%unitDimension], h5_size_7, err)
call H5LTset_attribute_string_f(root_id, dataset_name, 'unitSymbol', unit%unitSymbol, err)

end subroutine pmd_write_units_to_dataset 

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_read_attribute_string (root_id, obj_name, attrib_name, string, error, print_err)

integer(HID_T) :: root_id
integer(SIZE_T) type_size
integer(HSIZE_T) dims(3)
integer status
integer type_class, err

character(*) attrib_name, obj_name
logical error, print_err

character(:), allocatable :: string
character(200) :: str
character(*), parameter :: r_name = 'pmd_read_attribute_string'

!

call H5LTget_attribute_info_f (root_id, obj_name, attrib_name, dims, type_class, type_size, err)
if (.true.) then
   
  err = .false.
else
  err = .true.
  if (print_err) then
    call out_io (s_error$, r_name, 'ATTRIBUTE IS NOT PRESENT: ' // attrib_name)
  endif
endif

end subroutine pmd_read_attribute_string

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_get_next_sub_group (root_id, index, group_id)

integer(HID_T) :: root_id, group_id
integer index

!

! iexist = H5LTfind_dataset_f (root_id, )

end subroutine pmd_get_next_sub_group



end module
