!+
! Subroutine compute_reference_energy (lat, compute)
!
! Subroutine to compute the energy and momentum of the reference particle for 
! each element in a lat structure.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat     -- lat_struct: Input lattice.
!     %ele(0)%value(E_tot$) -- Energy at the start.
!   compute -- Logical, optional: If present then overrides the setting of
!                bmad_com%compute_ref_energy
!
!   bmad_com -- bmad_common_struct: Bmad global common block.
!     %compute_ref_energy -- Logical: If set False then do not recompute the
!                 reference energy.
!
! Output:
!   lat -- lat_struct
!     %ele(:)%value(E_tot$) -- Reference energy at the end of the element.
!     %ele(:)%value(p0c$)         -- Reference momentum at the end of the element.
!-

#include "CESR_platform.inc"

subroutine compute_reference_energy (lat, compute)

use bmad_struct
use bmad_utils_mod

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: ele, lord, slave
real(rp) E_tot, p0c, phase

integer i, j, k, ix, ixs
logical, optional :: compute

! Init energy

if (.not. logic_option(bmad_com%compute_ref_energy, compute)) return

E_tot = lat%ele(0)%value(E_tot$)
call convert_total_energy_to (E_tot, lat%param%particle, pc = p0c)
lat%ele(0)%value(p0c$) = p0c

! propagate the energy through the tracking part of the lattice

do i = 1, lat%n_ele_track
  ele => lat%ele(i)

  select case (ele%key)
  case (lcavity$) 
    ele%value(E_tot_START$) = E_tot
    ele%value(p0c_start$) = p0c

    phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
    E_tot = E_tot + ele%value(gradient$) * ele%value(l$) * cos(phase)
    E_tot = E_tot - ele%value(e_loss$) * lat%param%n_part * e_charge
    call convert_total_energy_to (E_tot, lat%param%particle, pc = p0c)

  case (custom$) 
    E_tot = E_tot + ele%value(gradient$) * ele%value(l$)
    call convert_total_energy_to (E_tot, lat%param%particle, pc = p0c)
  end select

  ele%value(E_tot$) = E_tot
  ele%value(p0c$) = p0c

enddo

! Put the appropriate energy values in the lord elements...

do i = lat%n_ele_track+1, lat%n_ele_max

  lord => lat%ele(i)
  if (lord%ix2_slave < 1) cycle  ! lord has no slaves

  ! Multipass lords have their own reference energy.
  ! If n_multipass_ref /= 0 then get this reference energy from the n^th slave.

  if (lord%control_type == multipass_lord$) then
    ix = nint(lord%value(n_multipass_ref$))
    if (ix /= 0) then
      j = lord%ix1_slave + ix - 1
      ixs = lat%control(j)%ix_slave
      lord%value(e_tot_ref_geometry$) = lat%ele(ixs)%value(e_tot$)
      lord%value(p0c_ref_geometry$)   = lat%ele(ixs)%value(p0c$)

    elseif (lord%value(e_tot_ref_geometry$) /= 0) then
      call convert_total_energy_to (lord%value(e_tot_ref_geometry$), &
                            lat%param%particle, pc = lord%value(p0c_ref_geometry$))
    endif

    lord%value(e_tot$) = lord%value(e_tot_ref_geometry$)
    lord%value(p0c$) = lord%value(p0c_ref_geometry$)

    cycle
  endif

  ! Now for everything but multipass_lord elements...
  ! The lord inherits the energy from the last slave.
  ! First find this slave.

  slave => lord
  do
    ix = slave%ix2_slave
    j = lat%control(ix)%ix_slave
    slave => lat%ele(j)
    if (j <= lat%n_ele_track) exit
  enddo

  ! Now transfer the information to the lord.

  lord%value(p0c$) = slave%value(p0c$)
  lord%value(E_tot$) = slave%value(E_tot$)

  ! Transfer the starting energy if needed.

  if (lord%key == lcavity$ .or. lord%key == custom$) then
    ix = lord%ix1_slave
    j = lat%control(ix)%ix_slave
    lord%value(E_tot_start$) = lat%ele(j)%value(E_tot_start$)
    lord%value(p0c_start$) = lat%ele(j)%value(p0c_start$)
  endif

enddo

end subroutine
