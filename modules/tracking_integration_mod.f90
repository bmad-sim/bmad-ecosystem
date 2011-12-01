!+
! Module for helper routines for tracking that is done step-by-step.
! Types of tracking that can use these routines:
!   runge_kutta
!   boris
!   symplectic integration
!-

module tracking_integration_mod

use bmad_struct
use bmad_interface

contains

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine lcavity_reference_energy_correction (ele, param, orbit)
!
! For elements where the reference energy is changing the reference energy in the body is 
! taken by convention to be the reference energy at the exit end.
! Elements where the reference energy can change:
!   lcavity
!   custom
!
! This routine should be called at the start of any tracking integration.
!
! Input:
!   ele       -- Ele_struct: Element being tracked through.
!   param     -- lat_param_struct:
!   orbit     -- Coord_struct: Coordinates to correct.
!
! Output:
!   orbit     -- Coord_struct: Coordinates to correct.
!-

subroutine lcavity_reference_energy_correction (ele, param, orbit)

implicit none

type (ele_struct) :: ele
type (lat_param_struct), intent(in) :: param
type (coord_struct) :: orbit

real(rp) p0, p1, e_start
real(rp), save :: vec6_start
character(40), parameter :: r_name = 'lcavity_reference_energy_correction'

!

if (ele_type_has_constant_reference_energy(ele%key)) return

orbit%vec(2:4:2) = orbit%vec(2:4:2) * ele%value(p0c_start$) / ele%value(p0c$)
orbit%vec(5) = orbit%vec(5) * (ele%value(p0c$) / ele%value(e_tot$)) / &
                          (ele%value(p0c_start$) / ele%value(e_tot_start$))
orbit%vec(6) = (1 + orbit%vec(6)) * ele%value(p0c_start$) / ele%value(p0c$) - 1

end subroutine lcavity_reference_energy_correction


end module
