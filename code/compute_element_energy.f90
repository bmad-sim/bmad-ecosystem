!+
! Subroutine compute_element_energy (ring)
!
! Subroutine to compute the energy of the reference particle for each 
! element in a ring structure.
! This routine is for linac lattices.
! The energy at the start of the "ring" is set to ring%param%energy.
!
! Modules needed:
!   use bmad
!
! Input:
!   ring -- Ring_struct: Input ring.
!
! Output:
!   ring -- Ring_struct
!     %ele_(0)%value(energy$) -- Set to ring%param%energy.
!     %ele_(:)%value(energy$) -- Energy at the end of the element.
!-

subroutine compute_element_energy (ring)

  use bmad_struct

  type (ring_struct) ring

  real(rdef) energy

  integer i

! compute element energy if a linac_lattice

  if (ring%param%lattice_type /= linac_lattice$) return

  energy = ring%param%beam_energy
  ring%ele_(0)%value(energy$) = energy

  do i = 1, ring%n_ele_ring
    if (ring%ele_(i)%key == linac_rf_cavity$) then
      ring%ele_(i)%value(energy_start$) = energy
      energy = energy + ring%ele_(i)%value(gradiant$) * &
          ring%ele_(i)%value(l$) * cos(twopi*ring%ele_(i)%value(phase_0$))
    endif
    ring%ele_(i)%value(energy$) = energy
  enddo

end subroutine
