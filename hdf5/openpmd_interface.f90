module openpmd_interface

use h5lt
use hdf5
use iso_fortran_env
use physical_constants

implicit none

! Common units

type pmd_unit_struct
 character(8) :: unitSymbol = ''     ! Native units name. EG 'eV'
 real(rp) :: unitSI = 0              ! Conversion to SI
 real(rp) :: unitDimension(7) = 0    ! SI Base Exponents
end type

type(pmd_unit_struct), parameter :: unit_1        = pmd_unit_struct('1', 1.0_rp, [0,0,0,0,0,0,0])
type(pmd_unit_struct), parameter :: unit_m        = pmd_unit_struct('m', 1.0_rp, [1,0,0,0,0,0,0])
type(pmd_unit_struct), parameter :: unit_kg       = pmd_unit_struct('kg', 1.0_rp, [0,1,0,0,0,0,0])
type(pmd_unit_struct), parameter :: unit_s        = pmd_unit_struct('s', 1.0_rp, [0,0,1,0,0,0,0])
type(pmd_unit_struct), parameter :: unit_A        = pmd_unit_struct('A', 1.0_rp, [0,0,0,1,0,0,0])
type(pmd_unit_struct), parameter :: unit_K        = pmd_unit_struct('K', 1.0_rp, [0,0,0,0,1,0,0])
type(pmd_unit_struct), parameter :: unit_mol      = pmd_unit_struct('mol', 1.0_rp, [0,0,0,0,0,1,0])
type(pmd_unit_struct), parameter :: unit_cd       = pmd_unit_struct('cd', 1.0_rp, [0,0,0,0,0,0,1])
type(pmd_unit_struct), parameter :: unit_Coulomb  = pmd_unit_struct('Coulomb', 1.0_rp, [0,0,1,1,0,0,0])
type(pmd_unit_struct), parameter :: unit_V_per_m  = pmd_unit_struct('V/m', 1.0_rp, [1,1,-3,-1,0,0,0])
type(pmd_unit_struct), parameter :: unit_clight   = pmd_unit_struct('clight', c_light, [1,0,-1,0,0,0,0])
type(pmd_unit_struct), parameter :: unit_m_per_s  = pmd_unit_struct('m/s', 1.0_rp, [1,0,-1,0,0,0,0])
type(pmd_unit_struct), parameter :: unit_eV       = pmd_unit_struct('eV', e_charge, [2,1,-2,0,0,0,0])
type(pmd_unit_struct), parameter :: unit_eV_per_c = pmd_unit_struct('eV/C', e_charge/c_light, [1,1,-1,0,0,0,0]) ! 5.3442859e-28

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

  ! BeamPhysics optional
  character(:), allocatable :: latticeFile
  character(:), allocatable :: latticeName
end type

! Misc

integer(hsize_t), parameter :: pmd_size_1 = 1
integer(hsize_t), parameter :: pmd_size_7 = 7

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_string_attrib(id, name, string)
integer(HID_T) :: id
character(*) :: name, string
integer error
!
call h5ltset_attribute_string_f(id, '.', name, trim(string), error)
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
call h5ltset_attribute_int_f(id, '.', name, [i], pmd_size_1, error)
end subroutine pmd_write_integer_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_real_attrib(id, name, x)
integer(HID_T) :: id
character(*) :: name
real(rp) :: x
integer error
call h5ltset_attribute_double_f(id, '.', name, [x], pmd_size_1, error)
end subroutine pmd_write_real_attrib

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_real_vector_to_dataset(root_id, dataset_name, unit, vector, error)

type (pmd_unit_struct) unit
real(rp) vector(:), v_max, v_min
integer error
integer(HID_T) :: root_id, v_size
character(*) dataset_name

!

v_max = maxval(vector)
v_min = minval(vector)

if (v_max == v_min) then
  call pmd_write_real_to_pseudo_dataset(root_id, dataset_name, unit, v_max, size(vector), error) 
  return
endif

!

v_size = size(vector)
call h5ltmake_dataset_double_f(root_id, dataset_name, 1, [v_size], vector, error)    

call h5ltset_attribute_double_f(root_id, dataset_name, 'minValue', [v_min], pmd_size_1, error)
call h5ltset_attribute_double_f(root_id, dataset_name, 'maxValue', [v_max], pmd_size_1, error)

call h5ltset_attribute_double_f(root_id, dataset_name, 'unitSI', [unit%unitSI], pmd_size_1, error)
call h5ltset_attribute_double_f(root_id, dataset_name, 'unitDimension', unit%unitDimension, pmd_size_7, error)
call h5ltset_attribute_string_f(root_id, dataset_name, 'unitSymbol', unit%unitSymbol, error)

end subroutine pmd_write_real_vector_to_dataset

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine pmd_write_real_to_pseudo_dataset(root_id, dataset_name, unit, value, v_size, error)

type (pmd_unit_struct) unit
integer(HID_T) :: root_id, group_id
real(rp) value
integer error, v_size
character(*) dataset_name

!

call h5gcreate_f(root_id, dataset_name, group_id, error)

call h5ltset_attribute_double_f(root_id, dataset_name, 'value', [value], pmd_size_1, error)
call h5ltset_attribute_int_f(root_id, dataset_name, 'shape', [v_size], pmd_size_1, error)

call h5ltset_attribute_double_f(root_id, dataset_name, 'unitSI', [unit%unitSI], pmd_size_1, error)
call h5ltset_attribute_double_f(root_id, dataset_name, 'unitDimension', unit%unitDimension, pmd_size_7, error)
call h5ltset_attribute_string_f(root_id, dataset_name, 'unitSymbol', unit%unitSymbol, error)

call h5gclose_f(group_id, error)

end subroutine pmd_write_real_to_pseudo_dataset

end module
