!+
! Subroutine compute_reference_energy (lattice, compute)
!
! Subroutine to compute the energy and momentum of the reference particle for 
! each element in a lat structure.
!
! Modules needed:
!   use bmad
!
! Input:
!   lattice -- lat_struct: Input lattice.
!     %ele(0)%value(E_TOT$) -- Energy at the start.
!   bmad_com -- bmad_common_struct: Bmad global common block.
!     %compute_ref_energy -- Logical: If set False then do not recompute the
!                 reference energy.
!   compute -- Logical, optional: If present then overrides the setting of
!                bmad_com%compute_ref_energy
!
! Output:
!   lattice -- lat_struct
!     %ele(:)%value(E_TOT$) -- Reference energy at the end of the element.
!     %ele(:)%value(p0c$)         -- Reference momentum at the end of the element.
!-

#include "CESR_platform.inc"

subroutine compute_reference_energy (lattice, compute)

  use bmad_struct
  use bmad_utils_mod

  implicit none

  type (lat_struct) lattice
  type (ele_struct), pointer :: ele, lord, slave
  real(rp) E_TOT, p0c, phase

  integer i, j, k, ix
  logical, optional :: compute

! Init energy

  if (.not. logic_option(bmad_com%compute_ref_energy, compute)) return

  E_TOT = lattice%ele(0)%value(E_TOT$)
  call convert_total_energy_to (E_TOT, lattice%param%particle, pc = p0c)
  lattice%ele(0)%value(p0c$) = p0c

! propagate the energy through the lattice

  do i = 1, lattice%n_ele_track
    ele => lattice%ele(i)

    select case (ele%key)
    case (lcavity$) 
      ele%value(E_TOT_START$) = E_TOT
      ele%value(p0c_start$) = p0c

      phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
      E_TOT = E_TOT + ele%value(gradient$) * &
                                                ele%value(l$) * cos(phase)
      E_TOT = E_TOT - &
                        ele%value(e_loss$) * lattice%param%n_part * e_charge
      call convert_total_energy_to (E_TOT, &
                                           lattice%param%particle, pc = p0c)

    case (custom$) 
      E_TOT = E_TOT + ele%value(gradient$) * ele%value(l$)
      call convert_total_energy_to (E_TOT, lattice%param%particle, pc = p0c)
    end select

    ele%value(E_TOT$) = E_TOT
    ele%value(p0c$) = p0c

  enddo

! Put energy in the lord elements. 

  do i = lattice%n_ele_track+1, lattice%n_ele_max

    lord => lattice%ele(i)
    if (lord%ix2_slave < 1) cycle  ! lord has no slaves

    slave => lord
    do
      ix = slave%ix2_slave
      j = lattice%control(ix)%ix_slave
      slave => lattice%ele(j)
      if (j <= lattice%n_ele_track) exit
    enddo

    lord%value(p0c$) = slave%value(p0c$)
    lord%value(E_TOT$) = slave%value(E_TOT$)

    if (lord%key == lcavity$ .or. lord%key == custom$) then
      ix = lord%ix1_slave
      j = lattice%control(ix)%ix_slave
      lord%value(E_TOT_START$) = lattice%ele(j)%value(E_TOT_START$)
      lord%value(p0c_start$) = lattice%ele(j)%value(p0c_start$)
    endif

  enddo

end subroutine
