module xraylib

use iso_c_binding

type crystal_struct
  real volume
end type


contains

function crystal_get_crystal(material) result (cryst)
character(*) material
type (crystal_struct), pointer :: cryst
end function

function crystal_f_h_structure_factor (cryst, E_kev, i, j, k, debye, angle) result (f_h)
type (crystal_struct) cryst
real(c_float) E_kev, debye, angle
integer i, j, k
complex(8) f0_tot
end function

function atomicnumbertosymbol(n) result (sym)
integer n
character(16) sym
end function

function atomicweight(n) result (weight)
integer n
real(c_float) weight
end function

function atomicdensity(n) result (density)
integer n
real(c_float) density
end function

subroutine atomic_factors (n, E_kev, q, debye, f0, fp, fpp) 
integer n
real(c_float) E_kev, q, debye, f0, fp, fpp
end subroutine


end module
