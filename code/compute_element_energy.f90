!+
! Subroutine compute_element_energy (lattice)
!
! Subroutine to compute the energy and momentum of the reference particle for 
! each element in a ring structure.
!
! Modules needed:
!   use bmad
!
! Input:
!   lattice -- Ring_struct: Input lattice.
!     %ele_(0)%value(beam_energy$) -- Energy at the start.
!
! Output:
!   lattice -- Ring_struct
!     %ele_(:)%value(beam_energy$) -- Energy at the end of the element.
!     %ele_(:)%value(p0c$)         -- Momentum at the end of the element.
!-

#include "CESR_platform.inc"

subroutine compute_element_energy (lattice)

  use bmad_struct
  use bmad_utils_mod

  implicit none

  type (ring_struct) lattice
  type (ele_struct), pointer :: ele, lord
  real(rp) beam_energy, p0c, phase

  integer i, j, k

! Init energy

  beam_energy = lattice%ele_(0)%value(beam_energy$)
  call energy_to_kinetic (beam_energy, lattice%param%particle, p0c = p0c)
  lattice%ele_(0)%value(p0c$) = p0c

! propagate the energy through the lattice

  do i = 1, lattice%n_ele_use
    ele => lattice%ele_(i)
    if (ele%key == lcavity$ .and. ele%is_on) then
      ele%value(energy_start$) = beam_energy
      phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
      beam_energy = beam_energy + ele%value(gradient$) * &
                                                  ele%value(l$) * cos(phase)
      if (bmad_com%sr_wakes_on) beam_energy = beam_energy - &
                            ele%value(e_loss$) * lattice%param%charge
      call energy_to_kinetic (beam_energy, lattice%param%particle, p0c = p0c)
        

    elseif (ele%key == custom$) then
      beam_energy = beam_energy + ele%value(delta_e$)
      call energy_to_kinetic (beam_energy, lattice%param%particle, p0c = p0c)
    endif

    ele%value(beam_energy$) = beam_energy
    ele%value(p0c$) = p0c

! If the ele is a super slave put the energy values in the lords also
! If the ele is also an lcavity then we want then energy at the end of
! the cavity

    if (ele%control_type == super_slave$) then
      do j = ele%ic1_lord, ele%ic2_lord
        k = lattice%ic_(j)
        lord => lattice%ele_(lattice%control_(k)%ix_lord)
        if (ele%key == lcavity$) then
          if (k == lord%ix1_slave) then  ! beginning of cavity
            lord%value(energy_start$) = ele%value(energy_start$)
            cycle
          elseif (k /= lord%ix2_slave) then
            cycle
          endif
        endif
        lord%value(beam_energy$) = beam_energy
        lord%value(p0c$) = p0c
      enddo
    endif

  enddo

end subroutine
