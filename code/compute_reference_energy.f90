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

use lat_ele_loc_mod
use rf_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave
type (branch_struct), pointer :: branch
type (coord_struct) start_orb, end_orb

real(rp) E_tot, p0c, phase

integer i, k, ib, ix, ixs
logical, optional :: compute
logical no_mode_adjust_done

! Init energy

if (.not. logic_option(bmad_com%compute_ref_energy, compute)) return

! propagate the energy through the tracking part of the lattice

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (ib /= 0) branch%ele(0)%value(E_tot$) = &
        lat%branch(branch%ix_from_branch)%ele(branch%ix_from_ele)%value(E_tot$)
  E_tot = branch%ele(0)%value(E_tot$)
  call convert_total_energy_to (E_tot, branch%param%particle, pc = p0c)
  branch%ele(0)%value(p0c$) = p0c

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)
    no_mode_adjust_done = .true.

    ! Calculate the energy

    select case (ele%key)
    case (lcavity$) 
      ele%value(E_tot_start$) = E_tot
      ele%value(p0c_start$) = p0c

      ! Different tracking methods can give different results so need to track if not bmad_standard tracking.
      if (ele%tracking_method == bmad_standard$) then
        phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
        E_tot = E_tot + ele%value(gradient$) * ele%value(l$) * cos(phase)
        E_tot = E_tot - e_loss_sr_wake (ele%value(e_loss$), branch%param)
        call convert_total_energy_to (E_tot, branch%param%particle, pc = p0c)
      else
        ! A zero e_tot can mess up tracking so put in a temp value if needed.
        if (ele%value(e_tot$) == 0) then ! Can happen on first pass through this routine
          ele%value(E_tot$) = E_tot      ! Temp value. Does not affect phase & amp adjustment.
          ele%value(p0c$) = p0c
          ele%ref_time = branch%ele(i-1)%ref_time
        endif
        if (associated(ele%rf%field)) call rf_accel_mode_adjust_phase_and_amp (ele, lat%param)
        no_mode_adjust_done = .false.
        call track1 (start_orb, ele, branch%param, end_orb)
        p0c = (1 + end_orb%vec(6)) * ele%value(p0c$)
        call convert_pc_to (p0c, branch%param%particle, e_tot = e_tot)
      endif

    case (custom$, hybrid$)
      ele%value(E_tot_start$) = E_tot
      ele%value(p0c_start$) = p0c
      E_tot = E_tot + ele%value(delta_e$)
      call convert_total_energy_to (E_tot, branch%param%particle, pc = p0c)

    case (crystal$, mirror$)
      ele%value(ref_wavelength$) = c_light * h_planck / E_tot

    case (patch$) 
      ele%value(E_tot_start$) = E_tot
      ele%value(p0c_start$) = p0c

      if (ele%is_on .and. ele%value(e_tot_offset$) /= 0) then
        e_tot = e_tot + ele%value(e_tot_offset$)
        call convert_total_energy_to (e_tot, branch%param%particle, pc = p0c)
      endif

    end select

    ele%value(E_tot$) = E_tot
    ele%value(p0c$) = p0c

    ! If there is an RF field map then scale the field to match.

    if (associated(ele%rf%field) .and. no_mode_adjust_done) call rf_accel_mode_adjust_phase_and_amp (ele, lat%param)

    ! Calculate the ref_time

    if (ele%key == lcavity$ .and. branch%ele(i-1)%value(E_tot$) /= E_tot) then
      if (ele%tracking_method == bmad_standard$) then
        ele%ref_time = branch%ele(i-1)%ref_time + ele%value(l$) * &        ! lcavity with non-zero acceleration formula
                  (p0c - branch%ele(i-1)%value(p0c$)) / ((E_tot - branch%ele(i-1)%value(E_tot$)) * c_light)
      else
        ele%ref_time = ele%ref_time - end_orb%vec(5) * E_tot / (p0c * c_light)
      endif
    elseif (ele%key == hybrid$) then
      ele%ref_time = branch%ele(i-1)%ref_time + ele%value(delta_ref_time$)
    else
      ele%ref_time = branch%ele(i-1)%ref_time + ele%value(l$) * E_tot / (p0c * c_light)
    endif

    ele%value(delta_ref_time$) = ele%ref_time - branch%ele(i-1)%ref_time

    ! %old_value is changed in tandem so changes in delta_ref_time do not trigger unnecessary bookkeeping.
    ele%old_value(delta_ref_time$) = ele%value(delta_ref_time$) 

  enddo
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
