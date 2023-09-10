!+
! Subroutine track_a_foil (orbit, ele, param, mat6, make_matrix)
!
! Bmad_standard tracking through an foil element.
!
! Input:
!   orbit       -- Coord_struct: Starting position.
!   ele         -- ele_struct: foil element.
!   param       -- lat_param_struct: Lattice parameters.
!   make_matrix -- logical, optional: Propagate the transfer matrix? Default is False.
!
! Output:
!   orbit      -- coord_struct: End position.
!   mat6(6,6)  -- real(rp), optional: Transfer matrix through the element.
!-

subroutine track_a_foil (orbit, ele, param, mat6, make_matrix)

use bmad_interface, except_dummy => track_a_foil
use xraylib, dummy => r_e
use random_mod

implicit none

type (coord_struct) :: orbit, orb0
type (ele_struct), target :: ele
type (lat_param_struct) :: param

real(rp), optional :: mat6(6,6)
real(rp) x0, xx0, sigma, z, rnd(2), density, p
real(rp), parameter :: S2 = 13.6, epsilon = 0.088

integer material, atomic_num

logical, optional :: make_matrix

character(*), parameter :: r_name = 'track_a_foil'

!

material = species_id(ele%component_name)
x0 = x0_radiation_length(material)   ! kg/m^2
if (x0 == real_garbage$) then
  call out_io(s_error$, r_name, 'CANNOT HANDLE NON-ATOMIC MATERIAL_TYPE: ' // species_name(material))
  orbit%state = lost$
  return
endif

atomic_num = atomic_number(material)
density = ElementDensity(atomic_num) * 1e3_rp  ! Convert to kg/m^3

!

z = atomic_number(orbit%species)
xx0 = ele%value(thickness$) / (x0 / density)
p = (1.0_rp + orbit%vec(6)) * orbit%p0c
! Factor of 1e6 is due to original formula using MeV/c for momentum
sigma = 1.0e6_rp * S2 * z * sqrt(xx0) / (p * orbit%beta) * (1.0_rp + epsilon * log10(xx0*z*z/orbit%beta**2))

!

call ran_gauss(rnd)
orbit%vec(2) = orbit%vec(2) + rnd(1) * sigma
orbit%vec(4) = orbit%vec(4) + rnd(2) * sigma

!

if (logic_option(.false., make_matrix)) then
  call mat_make_unit(mat6)
endif

end subroutine track_a_foil
