!+
! Subroutine compute_reference_energy (lattice, compute)
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
!   bmad_com -- Bmad_com_struct: Bmad global common block.
!     %compute_ref_energy -- Logical: If set False then do not recompute the
!                 reference energy.
!   compute -- Logical, optional: If present then overrides the setting of
!                bmad_com%compute_ref_energy
!
! Output:
!   lattice -- Ring_struct
!     %ele_(:)%value(beam_energy$) -- Reference energy at the end of the element.
!     %ele_(:)%value(p0c$)         -- Reference momentum at the end of the element.
!-

#include "CESR_platform.inc"

subroutine compute_reference_energy (lattice, compute)

  use bmad_struct
  use bmad_utils_mod

  implicit none

  type (ring_struct) lattice
  type (ele_struct), pointer :: ele, lord
  real(rp) beam_energy, p0c, phase

  integer i, j, k, ix
  logical, optional :: compute

! Init energy

  if (.not. logic_option(bmad_com%compute_ref_energy, compute)) return

  beam_energy = lattice%ele_(0)%value(beam_energy$)
  call energy_to_kinetic (beam_energy, lattice%param%particle, p0c = p0c)
  lattice%ele_(0)%value(p0c$) = p0c

! propagate the energy through the lattice

  do i = 1, lattice%n_ele_use
    ele => lattice%ele_(i)
    if (ele%key == lcavity$) then
      ele%value(energy_start$) = beam_energy
      ele%value(p0c_start$) = p0c
      if (ele%is_on) then
        phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
        beam_energy = beam_energy + ele%value(gradient$) * &
                                                  ele%value(l$) * cos(phase)
        if (bmad_com%sr_wakes_on) beam_energy = beam_energy - &
                   ele%value(e_loss$) * lattice%param%n_part * e_charge
        call energy_to_kinetic (beam_energy, lattice%param%particle, p0c = p0c)
      endif

    elseif (ele%key == custom$) then
      beam_energy = beam_energy + ele%value(gradient$) * ele%value(l$)
      call energy_to_kinetic (beam_energy, lattice%param%particle, p0c = p0c)
    endif

    ele%value(beam_energy$) = beam_energy
    ele%value(p0c$) = p0c

  enddo

! Put energy in the lord elements. 

  do i = lattice%n_ele_use+1, lattice%n_ele_max

    ele => lattice%ele_(i)
    do
      ix = ele%ix2_slave
      j = lattice%control_(ix)%ix_slave
      ele => lattice%ele_(j)
      if (j <= lattice%n_ele_use) exit
    enddo
    lattice%ele_(i)%value(p0c$) = ele%value(p0c$)
    lattice%ele_(i)%value(beam_energy$) = ele%value(beam_energy$)

    ele => lattice%ele_(i)
    if (ele%key == lcavity$ .or. ele%key == custom$) then
      ix = ele%ix1_slave
      j = lattice%control_(ix)%ix_slave
      ele%value(energy_start$) = lattice%ele_(j)%value(energy_start$)
      ele%value(p0c_start$) = lattice%ele_(j)%value(p0c_start$)
    endif

  enddo

end subroutine
