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
  use bmad_utils_mod

  implicit none

  type (ring_struct) ring

  real(rp) beam_energy, p0c

  integer i

! Init energy

  beam_energy = ring%ele_(0)%value(beam_energy$)
  call energy_to_kinetic (beam_energy, ring%param%particle, p0c = p0c)

! compute element energy if a linac_lattice

  if (ring%param%lattice_type == linac_lattice$) then

    do i = 1, ring%n_ele_use

      if (ring%ele_(i)%key == lcavity$ .and. ring%ele_(i)%is_on) then
        ring%ele_(i)%value(energy_start$) = beam_energy
        beam_energy = beam_energy + ring%ele_(i)%value(gradient$) * &
            ring%ele_(i)%value(l$) * cos(twopi*ring%ele_(i)%value(phi0$)) 
        if (bmad_com%sr_wakes_on) beam_energy = beam_energy - &
                            ring%ele_(i)%value(e_loss$) * ring%param%charge
        call energy_to_kinetic (beam_energy, ring%param%particle, p0c = p0c)
      endif

      ring%ele_(i)%value(beam_energy$) = beam_energy
      ring%ele_(i)%value(p0c$) = p0c

    enddo

! Otherwise everyone gets the same beam_energy

  else

    ring%ele_(:)%value(beam_energy$) = beam_energy
    ring%ele_(:)%value(p0c$)         = p0c

  endif

end subroutine
