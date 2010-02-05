!+
! Subroutine compute_reference_energy (lat, compute)
!
! Subroutine to compute the energy, momentum and time of the reference particle for 
! each element in a lat structure.
!
! Modules needed:
!   use bmad
!
! Input:
!   lat     -- lat_struct: Input lattice.
!     %ele(0)%value(E_tot$) -- Energy at the start of the lattice.
!   compute -- Logical, optional: If present then overrides the setting of
!                bmad_com%compute_ref_energy
!
!   bmad_com -- bmad_common_struct: Bmad global common block.
!     %compute_ref_energy -- Logical: If set False then do not recompute the
!                 reference energy.
!
! Output:
!   lat -- lat_struct
!     %ele(:)%value(E_tot$) -- Reference energy at the exit end.
!     %ele(:)%value(p0c$)   -- Reference momentum at the exit end.
!     %ele(:)%ref_time      -- Reference time from the beginning at the exit end.
!-

subroutine compute_reference_energy (lat, compute)

use bmad_struct
use lat_ele_loc_mod

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: ele, lord, slave
real(rp) E_tot, p0c, phase

integer i, k, ix, ixs
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
    ele%value(E_tot_start$) = E_tot
    ele%value(p0c_start$) = p0c

    phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
    E_tot = E_tot + ele%value(gradient$) * ele%value(l$) * cos(phase)
    E_tot = E_tot - ele%value(e_loss$) * lat%param%n_part * e_charge
    call convert_total_energy_to (E_tot, lat%param%particle, pc = p0c)

  case (custom$) 
    E_tot = E_tot + ele%value(gradient$) * ele%value(l$)
    call convert_total_energy_to (E_tot, lat%param%particle, pc = p0c)

  case (hybrid$)
    ele%value(E_tot_start$) = E_tot
    ele%value(p0c_start$) = p0c
    E_tot = E_tot + ele%value(delta_e$)
    call convert_total_energy_to (E_tot, lat%param%particle, pc = p0c)

  end select

  ele%value(E_tot$) = E_tot
  ele%value(p0c$) = p0c

  if (ele%key == lcavity$ .and. lat%ele(i-1)%value(E_tot$) /= E_tot) then
    ele%ref_time = lat%ele(i-1)%ref_time + ele%value(l$) * &
              (p0c - lat%ele(i-1)%value(p0c$)) / ((E_tot - lat%ele(i-1)%value(E_tot$)) * c_light)
  elseif (ele%key == hybrid$) then
    ele%ref_time = lat%ele(i-1)%ref_time + ele%value(delta_ref_time$)
  else
    ele%ref_time = lat%ele(i-1)%ref_time + ele%value(l$) * p0c / (E_tot * c_light)
  endif

enddo

! Put the appropriate energy values in the lord elements...

do i = lat%n_ele_track+1, lat%n_ele_max

  lord => lat%ele(i)
  if (lord%ix2_slave < 1) cycle  ! lord has no slaves

  ! Multipass lords have their own reference energy if n_ref_pass /= 0.

  if (lord%lord_status == multipass_lord$) then
    ix = nint(lord%value(n_ref_pass$))
    if (ix /= 0) then  
      slave => pointer_to_slave(lat, lord, ix)
      lord%value(e_tot$) = slave%value(e_tot$)
      lord%value(p0c$)   = slave%value(p0c$)
    elseif (lord%value(e_tot$) == 0 .and. lord%value(p0c$) /= 0) then
      call convert_pc_to (lord%value(p0c$), lat%param%particle, e_tot = lord%value(e_tot$))
    elseif (lord%value(p0c$) == 0 .and. lord%value(e_tot$) /= 0) then
      call convert_total_energy_to (lord%value(e_tot$), lat%param%particle, pc = lord%value(p0c$))
    endif
    cycle
  endif

  ! Now for everything but multipass_lord elements...
  ! The lord inherits the energy from the last slave.
  ! First find this slave.

  slave => lord
  do
    if (slave%n_slave == 0) exit
    slave => pointer_to_slave (lat, slave, slave%n_slave)
  enddo

  ! Now transfer the information to the lord.

  lord%value(p0c$) = slave%value(p0c$)
  lord%value(E_tot$) = slave%value(E_tot$)

  ! Transfer the starting energy if needed.

  if (lord%key == lcavity$ .or. lord%key == custom$) then
    slave => pointer_to_slave (lat, lord, 1)
    lord%value(E_tot_start$) = slave%value(E_tot_start$)
    lord%value(p0c_start$)   = slave%value(p0c_start$)
  endif

enddo

end subroutine
