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
use spin_mod

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

if (ele_has_constant_reference_energy(ele)) return

orbit%vec(2:4:2) = orbit%vec(2:4:2) * ele%value(p0c_start$) / ele%value(p0c$)
orbit%vec(5) = orbit%vec(5) * (ele%value(p0c$) / ele%value(e_tot$)) / &
                          (ele%value(p0c_start$) / ele%value(e_tot_start$))
orbit%vec(6) = (1 + orbit%vec(6)) * ele%value(p0c_start$) / ele%value(p0c$) - 1

end subroutine lcavity_reference_energy_correction

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!+
! Subroutine spin_track_a_step (ele, param, field, s_here, ds, orbit) 
!
! Routine to track the spin for one step of length ds.
!
! Input:
!   ele       -- Ele_struct: Element being tracked through.
!   param     -- lat_param_struct:
!   field     -- Em_field_struct: E & B fields.
!   s_here    -- Real(rp): Longitudinal position with respect to the element entrance end.
!   ds        -- Real(rp): Step size.
!   orbit     -- Coord_struct: Particle coordinates.
!     %spin     -- Spin coordinates
!
! Output:
!   orbit     -- Coord_struct: Particle coordinates
!     %spin     -- Spin coordinates
!-

subroutine spin_track_a_step (ele, param, field, s_here, ds, orbit) 

implicit none

type (ele_struct) :: ele
type (lat_param_struct), intent(in) :: param
type (coord_struct) :: orbit
type (em_field_struct) field

real(rp) s_here, ds
real(rp) :: Omega(3)
complex(rp) :: dspin_dz(2), quaternion(2,2)

! this uses a modified Omega' = -Omega/v_z

Omega = spin_omega_at (field, orbit, ele, param, s_here)
quaternion = (i_imaginary/2.0_rp)* (pauli(1)%sigma*Omega(1) + pauli(2)%sigma*Omega(2) + pauli(3)%sigma*Omega(3))

!   quaternion = normalized_quaternion (quaternion)

dspin_dz = matmul(quaternion, orbit%spin)
orbit%spin = orbit%spin + dspin_dz * ds

end subroutine spin_track_a_step

end module
