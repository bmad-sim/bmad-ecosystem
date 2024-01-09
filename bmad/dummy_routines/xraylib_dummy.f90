!+
! Module xraylib
!
! Dummy module for the xraylib interface.
!
! This can be used when the real xraylib is not easily available and X-rays are not being tracked.
! To use: Compile this file along with any program files.
!
! Note: This module is NOT compiled into the Bmad library.
! Note: This module is OBSOLETE since xraylib is now part of the Bmad Distribution.
!-

module xraylib

use iso_c_binding

implicit none

type crystal_struct
  real volume
end type

type compounddatanist
  integer nElements, elements(1)
  real massFractions(1), density
end type

TYPE, BIND(C) :: xrlComplex_C
        REAL (C_DOUBLE) :: re
        REAL (C_DOUBLE) :: im
ENDTYPE

real, parameter :: r_e = 0

!------------------------------------------------------------------

contains

function crystal_getcrystal(material) result (cryst)
character(*) material
type (crystal_struct), pointer :: cryst
nullify(cryst)
end function

function crystal_f_h_structurefactor (cryst, E_kev, i, j, k, debye, angle) result (f_h)
type (crystal_struct) cryst
real(c_double) E_kev, debye, angle
integer i, j, k
complex(8) f0_tot
complex(c_double) f_h
f_h = cmplx(0.0, 0.0)
end function

function atomicnumbertosymbol(n) result (sym)
integer n
character(16) sym
sym = ''
end function

function atomicweight(n) result (weight)
integer n
real(c_double) weight
weight = 0
end function

function elementdensity(n) result (density)
integer n
real(c_double) density
density = 0
end function

function atomic_factors (n, E_kev, q, debye, f0, fp, fpp) result (int_err)
integer n, int_err
real(c_double) E_kev, q, debye, f0, fp, fpp
end subroutine

function elementdensity(n) result (density)
integer n
real(c_double) density
density = 0
end function

function GetCompoundDataNISTByIndex(n) result (compound)
integer n
type (compoundDataNIST), pointer :: compound
nullify (compound)
end function

function Crystal_dSpacing(cryst, hkl1, hkl2, hkl3) result (spacing)
type (crystal_struct) cryst
integer hkl1, hkl2, hkl3
real spacing
spacing = 0
end function

subroutine FreeCompoundDataNIST(compound)
type (compoundDataNIST), pointer :: compound
nullify (compound)
end subroutine

end module
