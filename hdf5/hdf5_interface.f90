module hdf5_interface

use h5lt
use hdf5
use iso_fortran_env
use physical_constants

implicit none

! Common units

type unit
 character(8) :: name = '' ! Unit name
 real(rp) :: factor = 0    ! Conversion to SI
 real(rp) :: exp(7) = 0    ! SI Base Exponents
end type

type(unit), parameter :: unit_1        = unit('1', 1.0_rp, [0,0,0,0,0,0,0])
type(unit), parameter :: unit_m        = unit('m', 1.0_rp, [1,0,0,0,0,0,0])
type(unit), parameter :: unit_kg       = unit('kg', 1.0_rp, [0,1,0,0,0,0,0])
type(unit), parameter :: unit_s        = unit('s', 1.0_rp, [0,0,1,0,0,0,0])
type(unit), parameter :: unit_A        = unit('A', 1.0_rp, [0,0,0,1,0,0,0])
type(unit), parameter :: unit_K        = unit('K', 1.0_rp, [0,0,0,0,1,0,0])
type(unit), parameter :: unit_mol      = unit('mol', 1.0_rp, [0,0,0,0,0,1,0])
type(unit), parameter :: unit_cd       = unit('cd', 1.0_rp, [0,0,0,0,0,0,1])
type(unit), parameter :: unit_Coulomb  = unit('Coulomb', 1.0_rp, [0,0,1,1,0,0,0])
type(unit), parameter :: unit_V_per_m  = unit('V/m', 1.0_rp, [1,1,-3,-1,0,0,0])
type(unit), parameter :: unit_clight   = unit('clight', c_light, [1,0,-1,0,0,0,0])
type(unit), parameter :: unit_m_per_s  = unit('m/s', 1.0_rp, [1,0,-1,0,0,0,0])
type(unit), parameter :: unit_eV       = unit('eV', e_charge, [2,1,-2,0,0,0,0])
type(unit), parameter :: unit_eV_per_c = unit('eV/C', e_charge/c_light, [1,1,-1,0,0,0,0]) ! 5.3442859e-28

contains

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------

subroutine hdf5_write_string_attrib(id, name, string)
integer(HID_T) :: id
character(*) :: name, string
integer error
call h5ltset_attribute_string_f(id, '.', name, trim(string), error)
end subroutine

subroutine hdf5_write_integer_attrib(id, name, i)
integer(HID_T) :: id
character(*) :: name
integer :: i
integer error
call h5ltset_attribute_int_f(id, '.', name, [i], 1_hsize_t, error)
end subroutine

subroutine hdf5_write_real_attrib(id, name, x)
integer(HID_T) :: id
character(*) :: name
real(rp) :: x
integer error
call h5ltset_attribute_double_f(id, '.', name, [x], 1_hsize_t, error)
end subroutine

end module
