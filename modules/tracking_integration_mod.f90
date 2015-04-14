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
! Subroutine spin_track_a_step (ele, field, ds, orbit) 
!
! Routine to track the spin for one step of length ds.
!
! Input:
!   ele       -- Ele_struct: Element being tracked through.
!   field     -- em_field_struct: Field at particle.
!   orbit     -- Coord_struct: Particle coordinates.
!     %spin     -- Spin coordinates
!
! Output:
!   orbit     -- Coord_struct: Particle coordinates
!     %spin     -- Spin coordinates
!-

subroutine spin_track_a_step (ele, field, ds, orbit) 

implicit none

type (ele_struct) :: ele
type (coord_struct) :: orbit
type (em_field_struct) field

real(rp) ds
real(rp) :: Omega(3)
complex(rp) :: dspin_dz(2), quaternion(2,2)

! this uses a modified Omega' = Omega/v_z

Omega = spin_omega (field, orbit, ele)
quaternion = -(i_imaginary/2.0_rp)* (pauli(1)%sigma*Omega(1) + pauli(2)%sigma*Omega(2) + pauli(3)%sigma*Omega(3))

!   quaternion = normalized_quaternion (quaternion)

dspin_dz = matmul(quaternion, orbit%spin)
orbit%spin = orbit%spin + dspin_dz * ds

end subroutine spin_track_a_step

end module
