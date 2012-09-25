!+
! Subroutine lat_compute_ref_energy_and_time (lat, err_flag)
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

subroutine lat_compute_ref_energy_and_time (lat, err_flag)

use lat_ele_loc_mod
use bookkeeper_mod
use multipass_mod
use rf_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, lord2, slave, branch_ele, ele0, gun_ele, ele_init
type (branch_struct), pointer :: branch
type (coord_struct) start_orb, end_orb

real(rp) pc

integer j, k, ie, ib, ix, ixs, ibb, ix_slave, ixl, ix_pass, n_links
integer ix_super_end, ix_e_gun

logical did_set, stale, err
logical, optional :: err_flag

character(24), parameter :: r_name = 'lat_compute_ref_energy_and_time'

! propagate the energy through the tracking part of the lattice

if (present(err_flag)) err_flag = .true.

do ib = 0, ubound(lat%branch, 1)

  branch => lat%branch(ib)
  ele_init => branch%ele(0)

  if (bmad_com%auto_bookkeeper) then
    stale = .true.
  else
    if (branch%param%bookkeeping_state%ref_energy /= stale$) cycle
    stale = .false.
  endif

  branch%param%bookkeeping_state%ref_energy = ok$

  ! Init energy at beginning of branch if needed.

  if (stale .or. ele_init%bookkeeping_state%ref_energy == stale$) then
    if (branch%ix_from_branch >= 0) then

      branch_ele => pointer_to_ele (lat, branch%ix_from_ele, branch%ix_from_branch)
      branch%param%particle = nint(branch_ele%value(particle$))
      branch%param%lattice_type = nint(branch_ele%value(lattice_type$))

      did_set = .false.

      if (branch_ele%value(E_tot_start$) == 0) then
        ele_init%value(E_tot$) = branch_ele%value(E_tot$)
        call convert_total_energy_to (ele_init%value(E_tot$), branch%param%particle, pc = ele_init%value(p0c$))
      else
        ele_init%value(E_tot$) = branch_ele%value(E_tot_start$)
        did_set = .true.
      endif

      if (branch_ele%value(p0c_start$) == 0) then
        ele_init%value(p0c$) = branch_ele%value(p0c$)
       call convert_pc_to (ele_init%value(p0c$), branch%param%particle, e_tot = ele_init%value(e_tot$))
      else
        ele_init%value(p0c$) = branch_ele%value(p0c_start$)
        did_set = .true.
      endif

      if (.not. did_set .and. mass_of(branch%param%particle) /= &
                              mass_of(lat%branch(branch_ele%ix_branch)%param%particle)) then
        call out_io (s_fatal$, r_name, &
          'E_TOT_START OR P0C_START MUST BE SET IN A BRANCHING ELEMENT IF THE PARTICLE IN ', &
          'THE "FROM" BRANCH IS DIFFERENT FROM THE PARTICLE IN THE "TO" BRANCH.', &
          'PROBLEM OCCURS WITH BRANCH ELEMENT: ' // branch_ele%name) 
        if (global_com%exit_on_error) call err_exit
        return
      endif

      stale = .true.
      ele_init%bookkeeping_state%ref_energy = ok$
      ele_init%time_ref_orb_out = 0
      ele_init%value(delta_ref_time$) = 0
      ele_init%value(ref_time_start$) = ele_init%ref_time

    endif
  endif

  if (ele_init%value(E_tot$) == 0) then
    ele_init%value(E_tot$) = ele_init%value(E_tot_start$)
    ele_init%value(p0c$) = ele_init%value(p0c_start$)
  elseif (ele_init%value(E_tot_start$) == 0) then
    ele_init%value(E_tot_start$) = ele_init%value(E_tot$)
    ele_init%value(p0c_start$) = ele_init%value(p0c$)
  endif

  if (ele_init%bookkeeping_state%ref_energy == stale$) ele_init%bookkeeping_state%ref_energy = ok$

  ! Look for an e_gun and if found then the starting energy must be computed accordingly.
  ! Remember that there may be markers or null_eles before an e_gun in the lattice but nothing else.

  ix_e_gun = 0
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key == marker$) cycle
    if (ele%key == null_ele$) cycle
    ! Test first non-marker, non-null_ele element...
    if (ele%slave_status == super_slave$) then
      do j = 1, ele%n_lord
        lord => pointer_to_lord(ele, j)
        if (lord%key == e_gun$) then
          gun_ele => lord
          ix_e_gun = ie
          exit
        endif
      enddo
    elseif (ele%key == e_gun$) then
      gun_ele => ele
      ix_e_gun = ie
    endif
    exit
  enddo

  if (ix_e_gun /= 0) then ! Have found an e_gun...

    if (lat%rf_auto_scale_amp) then
      ele_init%value(e_tot$) = ele_init%value(e_tot_start$) + gun_ele%value(voltage$)
      call convert_total_energy_to (ele_init%value(e_tot$), branch%param%particle, pc = ele_init%value(p0c$))
 
    else
      gun_ele%value(e_tot$) = 2 * mass_of(branch%param%particle)   ! Dummy numbers so can do tracking
      call convert_total_energy_to (gun_ele%value(e_tot$), branch%param%particle, pc = gun_ele%value(p0c$))
      gun_ele%value(e_tot_start$) = gun_ele%value(e_tot$)
      gun_ele%value(p0c_start$) = gun_ele%value(p0c$)
      call init_coord (start_orb, ele = gun_ele, at_exit_end = .false., particle = branch%param%particle)
      start_orb%vec(6) = (ele_init%value(p0c_start$) - gun_ele%value(p0c_start$)) / gun_ele%value(p0c_start$)
      call track1 (start_orb, gun_ele, branch%param, end_orb, ignore_radiation = .true.)
      if (.not. particle_is_moving_forward(end_orb, branch%param%particle)) then
        call out_io (s_fatal$, r_name, 'PARTICLE LOST IN TRACKING E_GUN: ' // gun_ele%name, &
                                       'CANNOT COMPUTE REFERENCE TIME & ENERGY.')
        if (global_com%exit_on_error) call err_exit
        return
      endif
      ele_init%value(p0c$) = (1 + end_orb%vec(6)) * gun_ele%value(p0c$)
      call convert_pc_to (ele_init%value(p0c$), branch%param%particle, e_tot = ele_init%value(e_tot$))
    endif

  endif

  !----------------------------
  ! Loop over all elements in the branch

  ix_super_end = 0  ! End index of current super_lord_region

  do ie = 1, branch%n_ele_track

    ele0 => branch%ele(ie-1)
    ele => branch%ele(ie)

    if (.not. stale .and. ele%bookkeeping_state%ref_energy /= stale$) cycle

    stale = .true.

    if (ele%key == branch$ .or. ele%key == photon_branch$) then
      ibb = nint(ele%value(ix_branch_to$))
      call set_ele_status_stale (lat%branch(ibb)%ele(0), ref_energy_group$)
    endif

    ! If we are in the middle of a super_lord region then the "zero" orbit is just the continuation
    ! of the "zero" orbit of the previous element. This is important in wigglers. If the fact
    ! that the "zero" orbit is not truely zero is not taken into account, splitting 
    ! wigglers would result in z-position shifts when tracking particles.

    if (ix_super_end < ie) then       ! If not in super_lord region...
      ele%time_ref_orb_in = 0   ! Want zero orbit except if this is an e_gun then must set pz.
      if (ie == ix_e_gun) then
        ele%time_ref_orb_in(6) = (ele_init%value(p0c_start$) - ele_init%value(p0c$)) / ele_init%value(p0c$)
      endif 

    else                              ! In super_lord region
      ele%time_ref_orb_in = ele0%time_ref_orb_out
    endif

    ! Find where the current super lord region ends.

    if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
      do j = 1, ele%n_lord
        lord => pointer_to_lord(ele, j)
        lord2 => pointer_to_slave(lord, lord%n_slave) ! last element of lord
        ix_super_end = max(ix_super_end, lord2%ix_ele)
      enddo
    endif

    ! If this element is the first super_slave or multipass_slave of a lord that needs to be auto phased,
    ! make sure that the lord has its RF phase and amplitude properly adjusted.

    if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$ .or. ele%slave_status == multipass_slave$) then 
      do ixl = 1, ele%n_lord
        lord => pointer_to_lord (ele, ixl, ix_slave = ix_slave)
        if (ix_slave /= 1) cycle
        if (lord%lord_status /= super_lord$ .and. lord%lord_status /= multipass_lord$) cycle
        if (lord%slave_status == multipass_slave$) lord => pointer_to_lord(lord, 1)
        if (lord%lord_status == multipass_lord$) then
          call multipass_chain(ele, ix_pass, n_links)
          if (ix_pass /= 1) cycle
        endif
        ! This adjusts the RF phase and amplitude.
        ! Note: Any multipass lord element where the reference energy is not constant must have n_ref_pass = 1.
        if (lord%lord_status == multipass_lord$ .and. lord%value(n_ref_pass$) == 0) then
          if (lord%value(p0c$) /= 0) then
            call convert_pc_to(lord%value(p0c$), branch%param%particle, lord%value(e_tot$))
          elseif (lord%value(e_tot$) /= 0) then
            call convert_total_energy_to(lord%value(e_tot$), branch%param%particle, pc = lord%value(p0c$))
          else
            ! Marker element, for example, may not set either p0c or e_tot.
            cycle
          endif
          call ele_compute_ref_energy_and_time (lord, branch%param, &
                                              lord%value(e_tot$), lord%value(p0c$), 0.0_rp, err)
        else
          call ele_compute_ref_energy_and_time (lord, branch%param, &
                                              ele0%value(e_tot$), ele0%value(p0c$), ele0%ref_time, err)
        endif
        if (err) return
        call control_bookkeeper (lat, lord)
      enddo
    endif

    ! Calculate the energy and reference time at the end of the present element.

    call ele_compute_ref_energy_and_time (ele, branch%param, ele0%value(e_tot$), ele0%value(p0c$), ele0%ref_time, err)
    if (err) return

    ele%bookkeeping_state%ref_energy = ok$
    call set_ele_status_stale (ele, attribute_group$)
    call set_lords_status_stale (ele, lat, ref_energy_group$)

  enddo

enddo ! Branch loop

! Put the appropriate energy values in the lord elements...

lat%param%bookkeeping_state%ref_energy = ok$

do ie = lat%n_ele_track+1, lat%n_ele_max

  lord => lat%ele(ie)

  if (.not. bmad_com%auto_bookkeeper .and. lord%bookkeeping_state%ref_energy /= stale$) cycle
  if (lord%n_slave == 0) cycle   ! Can happen with null_ele$ elements for example.

  call set_ele_status_stale (lord, attribute_group$)
  lord%bookkeeping_state%ref_energy = ok$

  ! Multipass lords have their enegy computed above.

  if (lord%lord_status == multipass_lord$) cycle

  ! Now for everything but multipass_lord elements...
  ! The lord inherits the energy from the last slave.
  ! First find this slave.

  slave => lord
  do
    if (slave%n_slave == 0) exit
    slave => pointer_to_slave(slave, slave%n_slave)
  enddo

  ! Now transfer the information to the lord.

  lord%value(p0c$)   = slave%value(p0c$)
  lord%value(E_tot$) = slave%value(E_tot$)
  lord%ref_time      = slave%ref_time

  ! Transfer the starting energy.

  if (lord%lord_status == super_lord$) then
    slave => pointer_to_slave(lord, 1)
    lord%value(E_tot_start$) = slave%value(E_tot_start$)
    lord%value(p0c_start$)   = slave%value(p0c_start$)
  endif

enddo ! Branch loop

lat%lord_state%ref_energy = ok$
if (present(err_flag)) err_flag = .false.

end subroutine lat_compute_ref_energy_and_time

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine ele_compute_ref_energy_and_time (ele, param, e_tot_start, p0c_start, ref_time_start, err_flag)
!
! Routine to compute the reference energy and reference time at the end of an element 
! given the reference enegy and reference time at the start of the element.
!
! Input:
!   ele            -- Ele_struct: Lattice element
!     %time_ref_orb_in  -- Starting orbit for ref time calc.
!   param          -- lat_Param_struct: Lattice parameters.
!   e_tot_start    -- Real(rp): Entrance end energy.
!   p0c_start      -- Real(rp): Entrance end momentum
!   ref_time_start -- Real(rp): Entrance end reference time
!   err_flag       -- Logical: Set true if there is an error. False otherwise.
!
! Output:
!   ele         -- Ele_struct: Lattice element with reference energy and time.
!     %time_ref_orb_out  -- Ending orbit for ref time calc.
!-

subroutine ele_compute_ref_energy_and_time (ele, param, e_tot_start, p0c_start, ref_time_start, err_flag)

use lat_ele_loc_mod
use rf_mod

implicit none

type (ele_struct) ele
type (lat_param_struct) :: param
type (coord_struct) orb_start, orb_end

real(rp) E_tot_start, p0c_start, ref_time_start, e_tot, p0c, phase
real(rp) old_delta_ref_time, old_p0c
integer key
logical err_flag, err, changed

character(32), parameter :: r_name = 'ele_compute_ref_energy_and_time'

!

err_flag = .true.
old_delta_ref_time = ele%value(delta_ref_time$)
old_p0c = ele%value(p0c$)

ele%value(E_tot_start$)    = E_tot_start
ele%value(p0c_start$)      = p0c_start
ele%value(ref_time_start$) = ref_time_start

ele%time_ref_orb_out = ele%time_ref_orb_in  ! This should be true when we don't have to track.

key = ele%key
if (key == em_field$ .and. ele%sub_key == nonconst_ref_energy$) key = lcavity$

select case (key)

case (lcavity$)

  ! A zero e_tot (Can happen on first pass through this routine) can mess up tracking so put 
  ! in a temp value if needed. This does not affect the phase & amp adjustment.
  if (ele%value(e_tot$) == 0) then ! 
    ele%value(E_tot$) = ele%value(e_tot_start$)      
    ele%value(p0c$) = ele%value(p0c_start$)
    ele%ref_time = ref_time_start
  endif

  if (ele%slave_status /= super_slave$ .and. ele%slave_status /= slice_slave$ .and. ele%slave_status /= multipass_slave$) then
    call rf_auto_scale_phase_and_amp (ele, param, err)
    if (err) return
  endif

  ! Track

  call track_this_ele (.false.)

  ele%value(p0c$) = ele%value(p0c$) * (1 + orb_end%vec(6))
  call convert_pc_to (ele%value(p0c$), param%particle, E_tot = ele%value(E_tot$), err_flag = err)
  if (err) return

  call calc_time_ref_orb_out

case (custom$, hybrid$)
  ele%value(E_tot$) = E_tot_start + ele%value(delta_e$)
  call convert_total_energy_to (ele%value(E_tot$), param%particle, pc = ele%value(p0c$), err_flag = err)
  if (err) return

  ele%ref_time = ref_time_start + ele%value(delta_ref_time$)

case (e_gun$)
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start

  call track_this_ele (.true.)
  call calc_time_ref_orb_out

case (crystal$, mirror$, multilayer_mirror$)
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start

  ele%value(ref_wavelength$) = c_light * h_planck / E_tot_start
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

  ele%ref_time = ref_time_start + ele%value(t_offset$)

case default
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start

  if (ele%key == rfcavity$ .and. ele%slave_status /= super_slave$ .and. &
                        ele%slave_status /= slice_slave$ .and. ele%slave_status /= multipass_slave$) then
    call rf_auto_scale_phase_and_amp (ele, param, err)
    if (err) return
  endif

  if (ele_has_constant_ds_dt_ref(ele)) then
    if (ele%value(l$) == 0) then     ! Must avoid problem with zero length markers and p0c = 0.
      ele%ref_time = ref_time_start
    else
      ele%ref_time = ref_time_start + ele%value(l$) * E_tot_start / (p0c_start * c_light)
    endif

  else
    call track_this_ele (.false.)
    call calc_time_ref_orb_out
  endif

end select

! If delta_ref_time has shifted then any taylor map must be updated.

ele%value(delta_ref_time$) = ele%ref_time - ref_time_start

if (abs(ele%value(delta_ref_time$) - old_delta_ref_time) > bmad_com%significant_length / c_light) then
  if (associated (ele%taylor(1)%term) .and. ele%key /= taylor$) call kill_taylor (ele%taylor)
  ele%bookkeeping_state%mat6 = stale$
endif

! %old_value(delta_ref_time$) is changed in tandem so changes in delta_ref_time do not trigger unnecessary bookkeeping.
! However changes in the reference energy should trigger bookkeeping

ele%old_value(delta_ref_time$) = ele%value(delta_ref_time$) 
ele%old_value(p0c$) = old_p0c

err_flag = .false.

!---------------------------------------------------------------------------------
contains

subroutine track_this_ele (is_inside)

logical is_inside, auto_bookkeeper_saved

! Set auto_bookkeeper to prevent track1 calling attribute_bookkeeper which will overwrite
! ele%old_value

auto_bookkeeper_saved = bmad_com%auto_bookkeeper
bmad_com%auto_bookkeeper = .false.

call zero_errors_in_ele (ele, changed)
call init_coord (orb_start, ele%time_ref_orb_in, ele, .false., param%particle)
if (is_inside) orb_start%location = inside$ !to avoid entrance kick in time tracking
call track1 (orb_start, ele, param, orb_end, ignore_radiation = .true.)
if (.not. particle_is_moving_forward(orb_end, param%particle)) then
  call out_io (s_fatal$, r_name, 'PARTICLE LOST IN TRACKING: ' // ele%name, &
                                 'CANNOT COMPUTE REFERENCE TIME & ENERGY.')
  if (global_com%exit_on_error) call err_exit
  return
endif
call restore_errors_in_ele (ele)

ele%ref_time = ref_time_start + (orb_end%t - orb_start%t)

bmad_com%auto_bookkeeper = auto_bookkeeper_saved

end subroutine track_this_ele

!---------------------------------------------------------------------------------
! contains

recursive subroutine zero_errors_in_ele (ele, changed)

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer i
logical changed, has_changed

! For reference energy tracking need to turn off any element offsets and kicks and zero any errors.
! If the element is a super_slave then the errors must be zeroed in the super_lord elements also.

ele%old_value = ele%value
ele%bmad_logic = ele%is_on
ele%is_on = .true.
ele%ix_value = ele%tracking_method

has_changed = .false.

if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status /= super_lord$) cycle
    call zero_errors_in_ele (lord, changed)
    if (changed) has_changed = .true.
  enddo
endif

if (ele_has_offset(ele)) then
  call zero_ele_offsets (ele)
  has_changed = .true.
endif

if (ele_has_nonzero_kick(ele)) then
  call zero_ele_kicks (ele)
  has_changed = .true.
endif

select case (ele%key)
case (lcavity$)
  if (ele%value(phi0_err$) /= 0) then
    ele%value(phi0_err$) = 0
    has_changed = .true.
  endif
  if (ele%value(gradient_err$) /= 0) then
    ele%value(gradient_err$) = 0
    has_changed = .true.
  endif
end select

! For speed, use symp_lie_bmad tracking if the taylor map does not exist or if the taylor
! map is not valid for the element with kicks and offsets removed.

changed = has_changed
if (ele%tracking_method == taylor$) then
  if (.not. associated (ele%taylor(1)%term) .or. (changed .and. ele%map_with_offsets)) &
                                                       ele%tracking_method = symp_lie_bmad$
endif

end subroutine zero_errors_in_ele

!---------------------------------------------------------------------------------
! contains

recursive subroutine restore_errors_in_ele (ele)

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer i

! 

ele%value = ele%old_value
ele%is_on = ele%bmad_logic
ele%tracking_method = ele%ix_value

if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status /= super_lord$) cycle
    call restore_errors_in_ele (lord)
  enddo
endif

end subroutine restore_errors_in_ele

!---------------------------------------------------------------------------------
! contains

subroutine calc_time_ref_orb_out ()

! The tracking did not have the correct delta_ref_time end exit end ref energy so need to correct for this.
! Notice that here the delta_ref_time value (but not ele%value(p0c$)) the value used in tracking and not the 
! corrected value computed later.

ele%time_ref_orb_out = orb_end%vec
ele%time_ref_orb_out(2) = ele%time_ref_orb_out(2) / (1 + orb_end%vec(6))
ele%time_ref_orb_out(4) = ele%time_ref_orb_out(4) / (1 + orb_end%vec(6))
ele%time_ref_orb_out(5) = ele%time_ref_orb_out(5) + &
            (orb_end%t - orb_start%t - ele%value(delta_ref_time$)) * orb_end%beta * c_light
ele%time_ref_orb_out(6) = ((1 + orb_end%vec(6)) * orb_end%p0c - ele%value(p0c$)) / ele%value(p0c$)

end subroutine

end subroutine ele_compute_ref_energy_and_time
