!+
! Subroutine compute_element_energy (ring)
!
! Subroutine to compute the energy of the reference particle for each 
! element in a ring structure.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: Input ring.
!     %ele_(0)%value(beam_energy$) -- Energy at the start.
!
! Output:
!   ring -- Ring_struct
!     %ele_(:)%value(beam_energy$) -- Energy at the end of the element.
!-

#include "CESR_platform.inc"

subroutine compute_element_energy (ring)

  use bmad_struct

  type (ring_struct) ring

  real(rp) beam_energy, charge

  integer i

! compute element energy if a linac_lattice

  if (ring%param%lattice_type == linac_lattice$) then

    beam_energy = ring%ele_(0)%value(beam_energy$)
    ring%param%beam_energy = beam_energy

    charge = ring%param%charge  ! save charge
    ring%param%charge = 0

    do i = 1, ring%n_ele_ring
      if (ring%ele_(i)%key == lcavity$ .and. ring%ele_(i)%is_on) then
        ring%ele_(i)%value(energy_start$) = beam_energy
        beam_energy = beam_energy + ring%ele_(i)%value(gradient$) * &
            ring%ele_(i)%value(l$) * cos(twopi*ring%ele_(i)%value(phi0$))
      endif
      ring%ele_(i)%value(beam_energy$) = beam_energy
    enddo

    ring%param%charge = charge

    return

  endif

! Otherwise everyone gets the same beam_energy

  ring%ele_(0:ring%n_ele_max)%value(beam_energy$) = ring%param%beam_energy


end subroutine
