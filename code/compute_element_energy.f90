!+
! Subroutine compute_element_energy (ring)
!
! Subroutine to compute the energy of the reference particle for each 
! element in a ring structure.
! This routine is for linac lattices.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: Input ring.
!     %ele_(0)%value(energy$) -- Energy at the start.
!
! Output:
!   ring -- Ring_struct
!     %ele_(:)%value(energy$) -- Energy at the end of the element.
!-

#include "CESR_platform.inc"

subroutine compute_element_energy (ring)

  use bmad_struct

  type (ring_struct) ring

  real(rp) energy

  integer i

! compute element energy if a linac_lattice

  if (ring%param%lattice_type /= linac_lattice$) return

  energy = ring%ele_(0)%value(energy$)

  do i = 1, ring%n_ele_ring
    if (ring%ele_(i)%key == lcavity$) then
      ring%ele_(i)%value(energy_start$) = energy
      energy = energy + ring%ele_(i)%value(gradient$) * &
          ring%ele_(i)%value(l$) * cos(twopi*ring%ele_(i)%value(phi0$))
    endif
    ring%ele_(i)%value(energy$) = energy
  enddo

end subroutine
