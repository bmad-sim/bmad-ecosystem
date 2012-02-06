!+
! Subroutine lat_compute_reference_energy (lat, err_flag)
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
!
! Output:
!   lat      -- lat_struct
!     %ele(:)%value(E_tot$) -- Reference energy at the exit end.
!     %ele(:)%value(p0c$)   -- Reference momentum at the exit end.
!     %ele(:)%ref_time      -- Reference time from the beginning at the exit end.
!   err_flag -- Logical, optional: Set true if there is an error. False otherwise.
!-

subroutine lat_compute_reference_energy (lat, err_flag)

use lat_ele_loc_mod
use bookkeeper_mod
use multipass_mod
use rf_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, slave, branch_ele, ele0
type (branch_struct), pointer :: branch

integer i, k, ib, ix, ixs, ibb, ix_slave, ixl, ix_pass, n_links

logical did_set, stale, err
logical, optional :: err_flag

character(24), parameter :: r_name = 'lat_compute_reference_energy'

! propagate the energy through the tracking part of the lattice

if (present(err_flag)) err_flag = .true.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  if (bmad_com%auto_bookkeeper) then
    stale = .true.
  else
    if (branch%param%status%ref_energy /= stale$) cycle
    stale = .false.
  endif

  branch%param%status%ref_energy = ok$

  ! Init energy at beginning of branch if needed.

  ele0 => branch%ele(0)

  if (branch%ix_from_branch >= 0 .and. (stale .or. ele0%status%ref_energy == stale$)) then

    branch_ele => pointer_to_ele (lat, branch%ix_from_ele, branch%ix_from_branch)
    branch%param%particle = nint(branch_ele%value(particle$))
    branch%param%lattice_type = nint(branch_ele%value(lattice_type$))

    did_set = .false.

    if (branch_ele%value(E_tot_start$) == 0) then
      ele0%value(E_tot$) = branch_ele%value(E_tot$)
      call convert_total_energy_to (ele0%value(E_tot$), branch%param%particle, pc = ele0%value(p0c$))
    else
      ele0%value(E_tot$) = branch_ele%value(E_tot_start$)
      did_set = .true.
    endif

    if (branch_ele%value(p0c_start$) == 0) then
      ele0%value(p0c$) = branch_ele%value(p0c$)
     call convert_pc_to (ele0%value(p0c$), branch%param%particle, e_tot = ele0%value(e_tot$))
    else
      ele0%value(p0c$) = branch_ele%value(p0c_start$)
      did_set = .true.
    endif

    if (.not. did_set .and. mass_of(branch%param%particle) /= &
                            mass_of(lat%branch(branch_ele%ix_branch)%param%particle)) then
      call out_io (s_fatal$, r_name, &
        'E_TOT_START OR P0C_START MUST BE SET IN A BRANCHING ELEMENT IF THE PARTICLE IN ', &
        'THE "FROM" BRANCH IS DIFFERENT FROM THE PARTICLE IN THE "TO" BRANCH.', &
        'PROBLEM OCCURS WITH BRANCH ELEMENT: ' // branch_ele%name) 
      if (bmad_status%exit_on_error) call err_exit
      return
    endif

    stale = .true.
    ele0%status%ref_energy = ok$

  endif

  ! Loop over all elements in the branch

  do i = 1, branch%n_ele_track
    ele => branch%ele(i)

    if (.not. stale .and. ele%status%ref_energy /= stale$) cycle

    stale = .true.
    ele%status%ref_energy = ok$

    if (ele%key == branch$ .or. ele%key == photon_branch$) then
      ibb = nint(ele%value(ix_branch_to$))
      call set_ele_status_stale (lat%branch(ibb)%ele(0), lat%branch(ibb)%param, ref_energy_group$)
    endif

    ! If this element is the first super_slave or multipass_slave of a lord with varying reference energy,
    ! then, if needed, make sure that the lord has its RF phase and amplitude properly adjusted.

    ele0 => branch%ele(i-1)

    if ((ele%slave_status == super_slave$ .or. ele%slave_status == multipass_slave$) .and. &
                                            .not. ele_has_constant_reference_energy(ele)) then
      do ixl = 1, ele%n_lord
        lord => pointer_to_lord (ele, ixl, ix_slave = ix_slave)
        if (ele_has_constant_reference_energy(lord) .or. ix_slave /= 1) cycle
        if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)
        if (lord%tracking_method == bmad_standard$) cycle
        if (lord%lord_status == multipass_lord$) then
          call multipass_chain(ele, lat, ix_pass, n_links)
          if (ix_pass /= 1) cycle
        endif
        ! This adjusts the RF phase and amplitude
        call compute_ele_reference_energy (lord, branch%param, &
                                              ele0%value(e_tot$), ele0%value(p0c$), ele0%ref_time, err)
        if (err) return
        call control_bookkeeper (lat, lord)
      enddo
    endif

    ! Calculate the energy at the end of the present element.

    call compute_ele_reference_energy (ele, branch%param, ele0%value(e_tot$), ele0%value(p0c$), ele0%ref_time, err)
    if (err) return

    call set_ele_status_stale (ele, branch%param, attribute_group$)
    call set_lords_status_stale (ele, lat, ref_energy_group$)

  enddo

enddo

! Put the appropriate energy values in the lord elements...

lat%param%status%ref_energy = ok$

do i = lat%n_ele_track+1, lat%n_ele_max

  lord => lat%ele(i)

  if (.not. bmad_com%auto_bookkeeper .and. lord%status%ref_energy /= stale$) cycle

  call set_ele_status_stale (lord, lat%param, attribute_group$)
  lord%status%ref_energy = ok$

  ! Multipass lords have their own reference energy if n_ref_pass /= 0.

  if (lord%lord_status == multipass_lord$) then
    ix = nint(lord%value(n_ref_pass$))
    if (ix /= 0) then  
      slave => pointer_to_slave(lord, ix)
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
    slave => pointer_to_slave(slave, slave%n_slave)
  enddo

  ! Now transfer the information to the lord.

  lord%value(p0c$) = slave%value(p0c$)
  lord%value(E_tot$) = slave%value(E_tot$)

  ! Transfer the starting energy.

  !!if (lord%key == lcavity$ .or. lord%key == custom$) then
    slave => pointer_to_slave(lord, 1)
    lord%value(E_tot_start$) = slave%value(E_tot_start$)
    lord%value(p0c_start$)   = slave%value(p0c_start$)
  !!endif

  ! Autophase rfcavity lords.

  if (lord%key == rfcavity$) then
    slave => pointer_to_slave(lord, 1)
    call rf_auto_scale_phase_and_amp (lord, lat%branch(slave%ix_branch)%param, err)
    if (err) return
  endif

enddo

if (present(err_flag)) err_flag = .false.

end subroutine lat_compute_reference_energy

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine compute_ele_reference_energy (ele, param, e_tot_start, p0c_start, ref_time_start, err_flag)
!
! Routine to compute the reference energy and reference time at the end of an element 
! given the reference enegy and reference time at the start of the element.
!
! Input:
!   ele            -- Ele_struct: Lattice element
!   param          -- lat_Param_struct: Lattice parameters.
!   e_tot_start    -- Real(rp): Entrance end energy.
!   p0c_start      -- Real(rp): Entrance end momentum
!   ref_time_start -- Real(rp): Entrance end reference time
!   err_flag       -- Logical: Set true if there is an error. False otherwise.
!
! Output:
!   ele         -- Ele_struct: Lattice element with reference energy and time.
!-

subroutine compute_ele_reference_energy (ele, param, e_tot_start, p0c_start, ref_time_start, err_flag)

use lat_ele_loc_mod
use rf_mod

implicit none

type (ele_struct) ele
type (ele_struct), save :: ele2
type (lat_param_struct) :: param
type (coord_struct) start_orb, end_orb

real(rp) E_tot_start, p0c_start, ref_time_start, e_tot, p0c, phase
integer key
logical err_flag, err

character(32), parameter :: r_name = 'compute_ele_reference_energy'

! Treat an accelerating em_field element like an lcavity

err_flag = .true.

key = ele%key
if (ele%key == em_field$ .and. .not. ele_has_constant_reference_energy(ele)) key = lcavity$

ele%value(E_tot_start$) = E_tot_start
ele%value(p0c_start$) = p0c_start

select case (key)
case (lcavity$) 

  ! We can only use the formula dE = voltage * cos(phase) with bmad_standard$ tracking since with other 
  ! tracking there is no guarantee that dE varies as cos(phase). Additionally, for multipass elements 
  ! with tracking /= bmad_standard$, dE at phase = 0 may not even be equal to the voltage if this
  ! is not the reference pass. 

  if (ele%tracking_method == bmad_standard$) then
    phase = twopi * (ele%value(phi0$) + ele%value(dphi0$)) 
    E_tot = E_tot_start + ele%value(gradient$) * ele%value(l$) * cos(phase)
    call convert_total_energy_to (E_tot, param%particle, pc = p0c, err_flag = err)
    if (err) return
    ele%value(E_tot$) = E_tot
    ele%value(p0c$) = p0c

    if (E_tot_start == E_tot) then
      ele%ref_time = ref_time_start + ele%value(l$) * E_tot / (p0c * c_light)
    else
      ele%ref_time = ref_time_start + ele%value(l$) * &  ! lcavity with non-zero acceleration formula
                (p0c - p0c_start) / ((E_tot - E_tot_start) * c_light)
    endif

  else
    ! A zero e_tot (Can happen on first pass through this routine) can mess up tracking so put 
    ! in a temp value if needed. This does not affect the phase & amp adjustment.
    if (ele%value(e_tot$) == 0) then ! 
      ele%value(E_tot$) = ele%value(e_tot_start$)      
      ele%value(p0c$) = ele%value(p0c_start$)
      ele%ref_time = ref_time_start
    endif

    if (ele%slave_status /= super_slave$ .and. ele%slave_status /= multipass_slave$) then
      call rf_auto_scale_phase_and_amp (ele, param, err)
      if (err) return
    endif

    ! For reference energy tracking need to turn off any element offsets and kicks.
    ! If a super_slave, only want to track through the accelerating element.

    ele2 = ele
    if (ele2%slave_status == super_slave$) then
      
    else
      call zero_ele_offsets (ele2)
      call zero_ele_kicks (ele2)
    endif

    ele2%value(phi0_err$) = 0
    ele2%value(gradient_err$) = 0
    call track1 (start_orb, ele2, param, end_orb)
    if (end_orb%status == dead$) then
      call out_io (s_error$, r_name, 'PARTICLE LOST IN TRACKING LCAVITY: ' // ele%name, &
                                     'CANNOT COMPUTE REFERENCE ENERGY')
      return
    endif
    E_tot = ele%value(E_tot$)
    p0c = ele%value(p0c$)
    ele%ref_time = ref_time_start + end_orb%t
    ele%value(p0c$) = p0c * (1 + end_orb%vec(6))
    call convert_pc_to (ele%value(p0c$), param%particle, E_tot = ele%value(E_tot$), err_flag = err)
    if (err) return
  endif

case (custom$, hybrid$)
  ele%value(E_tot$) = E_tot_start + ele%value(delta_e$)
  call convert_total_energy_to (ele%value(E_tot$), param%particle, pc = ele%value(p0c$), err_flag = err)
  if (err) return

  ele%ref_time = ref_time_start + ele%value(delta_ref_time$)

case (e_gun$)
  ele%value(E_tot$) = E_tot_start + ele%value(voltage$)
  call convert_total_energy_to (ele%value(E_tot$), param%particle, pc = ele%value(p0c$), err_flag = err)
  if (err) return

  ele2 = ele
  call zero_ele_offsets (ele2)
  call zero_ele_kicks (ele2)

  start_orb%status = inside$ !to avoid entrance kick in time tracking
  call track1 (start_orb, ele2, param, end_orb)
  E_tot = ele%value(E_tot$)
  p0c = ele%value(p0c$)
  ele%ref_time = ref_time_start + ele%value(delta_ref_time$) - end_orb%vec(5) * E_tot / (p0c * c_light)

case (crystal$, mirror$, multilayer_mirror$)
  ele%value(ref_wavelength$) = c_light * h_planck / E_tot_start
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start
  ele%ref_time = ref_time_start

case (patch$) 
  if (ele%is_on .and. ele%value(e_tot_offset$) /= 0) then
    ele%value(E_tot$) = e_tot_start + ele%value(e_tot_offset$)
    call convert_total_energy_to (ele%value(E_tot$), param%particle, pc = ele%value(p0c$), err_flag = err)
    if (err) return
  else
    ele%value(E_tot$) = E_tot_start
    ele%value(p0c$) = p0c_start
  endif

  ele%ref_time = ref_time_start

case default
  if (ele%key == em_field$) then
    ele%value(E_tot_start$) = E_tot_start
    ele%value(p0c_start$) = p0c_start
  endif
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start
  ele%ref_time = ref_time_start + ele%value(l$) * E_tot_start / (p0c_start * c_light)

  if (ele%key == rfcavity$ .and. ele%slave_status /= super_slave$ .and. ele%slave_status /= multipass_slave$) then
    call rf_auto_scale_phase_and_amp (ele, param, err)
    if (err) return
  endif

end select


! %old_value is changed in tandem so changes in delta_ref_time do not trigger unnecessary bookkeeping.

ele%value(delta_ref_time$) = ele%ref_time - ref_time_start
ele%old_value(delta_ref_time$) = ele%value(delta_ref_time$) 

err_flag = .false.

end subroutine compute_ele_reference_energy
