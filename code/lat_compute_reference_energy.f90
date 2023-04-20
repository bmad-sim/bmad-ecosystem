!+
! Subroutine lat_compute_ref_energy_and_time (lat, err_flag)
!
! Subroutine to compute the energy, momentum and time of the reference particle for 
! each element in a lat structure.
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
!   err_flag -- Logical: Set true if there is an error. False otherwise.
!-

subroutine lat_compute_ref_energy_and_time (lat, err_flag)

use autoscale_mod, dummy2 => lat_compute_ref_energy_and_time

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele, lord, lord2, slave, fork_ele, ele0, gun_ele, begin_ele, ele2
type (branch_struct), pointer :: branch
type (coord_struct) start_orb, end_orb

real(rp) pc, abs_tol(3)
real(rp), parameter :: zero6(6) = 0

integer j, k, n, ie, ib, ix, ixs, ibb, ix_slave, ixl, ix_pass, n_links
integer ix_super_end, ix_e_gun, it

logical stale, err, lord_compute
logical :: err_flag

character(*), parameter :: r_name = 'lat_compute_ref_energy_and_time'

! propagate the energy through the tracking part of the lattice

err_flag = .true.

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)
  branch%param%unstable_factor = 0   ! Set to positive number if there is a problem.
  begin_ele => branch%ele(0)
  nullify(fork_ele)

  if (branch%ix_from_branch > -1 .and. branch%ix_to_ele == 0) then
    fork_ele => lat%branch(branch%ix_from_branch)%ele(branch%ix_from_ele)
    if (fork_ele%key == photon_fork$) branch%param%particle = photon$
  endif

  if (begin_ele%ref_species /= branch%param%particle) then
    begin_ele%ref_species = branch%param%particle
    begin_ele%bookkeeping_state%ref_energy = stale$
    branch%param%bookkeeping_state%ref_energy = stale$
  endif

  if (branch%param%bookkeeping_state%ref_energy /= stale$) cycle
  stale = (begin_ele%bookkeeping_state%ref_energy == stale$)

  ! Init energy at beginning of branch if needed.
  ! If there is a fork into this branch at element 0, the reference energy is may be inherited from the branch forked from.

  if (stale) then
    if (begin_ele%value(inherit_from_fork$) == real_garbage$) then  ! Happens first time this routine is called from bmad_parser. 
      begin_ele%value(inherit_from_fork$) = false$
      if (associated(fork_ele)) then
        if (begin_ele%ref_species == fork_ele%ref_species .or. begin_ele%ref_species == not_set$) &
                                                                           begin_ele%value(inherit_from_fork$) = true$
      endif
    endif

    ! If inherit from fork.
    if (is_true(begin_ele%value(inherit_from_fork$))) then
      begin_ele%value(E_tot_start$) = fork_ele%value(E_tot$)
      begin_ele%value(p0c_start$) = fork_ele%value(p0c$)
      begin_ele%value(E_tot$) = fork_ele%value(E_tot$)
      begin_ele%value(p0c$) = fork_ele%value(p0c$)
      begin_ele%ref_species = fork_ele%ref_species
      branch%param%particle = begin_ele%ref_species
      call init_coord (begin_ele%time_ref_orb_in, zero6, begin_ele, upstream_end$)
      call init_coord (begin_ele%time_ref_orb_out, zero6, begin_ele, downstream_end$)

    ! Else no inherit from fork.
    else
      if (begin_ele%ref_species == not_set$) then
        if (associated(fork_ele)) then
          begin_ele%ref_species = fork_ele%ref_species
          branch%param%particle = fork_ele%ref_species
        elseif (ib == 0) then
          begin_ele%ref_species = positron$
          branch%param%particle = begin_ele%ref_species
        else
          begin_ele%ref_species = lat%param%particle
          branch%param%particle = lat%param%particle
        endif
      endif

      if (begin_ele%value(p0c$) >= 0) then
        call convert_pc_to (begin_ele%value(p0c$), begin_ele%ref_species, e_tot = begin_ele%value(e_tot$))
      elseif (begin_ele%value(e_tot$) >= mass_of(begin_ele%ref_species)) then
        call convert_total_energy_to (begin_ele%value(e_tot$), begin_ele%ref_species, pc = begin_ele%value(p0c$))
      elseif (branch%ix_from_branch < 0 .and. branch%ix_branch > 0 .and. &
                              lat%param%particle == begin_ele%ref_species) then ! A root branch.
        begin_ele%value(e_tot$) = lat%ele(0)%value(e_tot$)
        begin_ele%value(p0c$)   = lat%ele(0)%value(p0c$)
      else
        if (begin_ele%value(e_tot$) < 0 .and. begin_ele%value(p0c$) < 0) then
          if (begin_ele%ref_species == photon$) then
            call out_io (s_fatal$, r_name, 'REFERENCE ENERGY IS NOT SET IN BRANCH: ' // branch%name)
          else
            call out_io (s_fatal$, r_name, 'REFERENCE ENERGY IS NOT SET IN BRANCH: ' // branch%name)
          endif
        else
          call out_io (s_fatal$, r_name, 'REFERENCE ENERGY IS SET BELOW MC^2 IN BRANCH ' // branch%name)
        endif
        if (global_com%exit_on_error) call err_exit
        return
      endif

      ! If E_tot or p0c is modified, set_flags_for_changed_attribute will set E_tot_start = E_tot and 
      ! p0c_start = E_tot. So only need to consider the case when starting E_tot/p0c is initialized.

      if (begin_ele%value(e_tot_start$) <= 0) then
        begin_ele%value(e_tot_start$) = begin_ele%value(e_tot$)
        begin_ele%value(p0c_start$) = begin_ele%value(p0c$)
      endif

      begin_ele%value(delta_ref_time$) = 0
      begin_ele%value(ref_time_start$) = begin_ele%ref_time
    endif

    begin_ele%bookkeeping_state%ref_energy = ok$

    !

    if (branch%param%geometry == 0) then   ! Not yet set
      if (branch%param%particle == photon$) then
        branch%param%geometry = open$
      elseif (any(branch%ele(:)%key == lcavity$)) then
        if (branch%ix_branch == 0) then
          call out_io (s_warn$, r_name, 'NOTE: THIS LATTICE HAS AN LCAVITY.', 'SETTING THE GEOMETRY TO OPEN.')
        else
          call out_io (s_warn$, r_name, 'NOTE: BRANCH ' // trim(branch%name) // ' HAS AN LCAVITY.', 'SETTING THE GEOMETRY TO OPEN.')
        endif
        branch%param%geometry = open$
      else
        branch%param%geometry = closed$
      endif
    endif

  endif    ! if (stale)

  !---------------
  ! Look for an e_gun.
  ! Remember that there may be markers or null_eles before an e_gun in the lattice but nothing else.

  ix_e_gun = 0
  do ie = 1, branch%n_ele_track
    ele => branch%ele(ie)
    if (ele%key == marker$ .or. ele%key == detector$) cycle
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

  ! If there is an e_gun then compute the energy at the exit. 
  ! This must be done by tracking since with autoscale off or if the field is not DC, the 
  ! voltage is not a reliable number.

  if (ix_e_gun /= 0) then ! Have found an e_gun...
    do j = 1, ix_e_gun  ! Also mark marker elements before gun
      branch%ele(j)%value(e_tot_ref_init$) = begin_ele%value(e_tot_start$)
      branch%ele(j)%value(p0c_ref_init$) = begin_ele%value(p0c_start$)
    enddo
    gun_ele%value(e_tot_ref_init$) = begin_ele%value(e_tot_start$) ! In case gun is a super_lord.
    gun_ele%value(p0c_ref_init$) = begin_ele%value(p0c_start$)
    gun_ele%ref_species = begin_ele%ref_species

    if (gun_ele%value(e_tot_ref_init$) + gun_ele%value(voltage$) < mass_of(default_tracking_species(branch%param)) .and. &
        (is_true(gun_ele%value(autoscale_amplitude$)) .or. gun_ele%tracking_method == bmad_standard$)) then
      call out_io (s_fatal$, r_name, '(INITIAL ENERGY) + (E_GUN VOLTAGE) MUST BE NON-NEGATIVE! ' // gun_ele%name, &
                                     'CANNOT COMPUTE REFERENCE TIME & ENERGY.')
      if (global_com%exit_on_error) call err_exit
      return
    endif

    ! p0c_start and p0c, need to be set for tracking and they need to be nonzero.
    ! Since p0c_ref_init and the voltage may both be zero, just use a dummy number in this case.
    if (gun_ele%value(p0c$) == 0) then
      gun_ele%value(p0c$) = 1d3 + abs(gun_ele%value(voltage$))
    endif
    call convert_pc_to (gun_ele%value(p0c$), branch%param%particle, E_tot= gun_ele%value(E_tot$))
    gun_ele%value(e_tot_start$) = gun_ele%value(e_tot$)
    gun_ele%value(p0c_start$)   = gun_ele%value(p0c$)

    do it = 1, 3
      call autoscale_phase_and_amp (gun_ele, branch%param, err, call_bookkeeper = .false.)
      if (err) then
        call out_io (s_fatal$, r_name, 'AUTOSCALE FAILED FOR ELEMENT: ' // ele%name)
        if (global_com%exit_on_error) call err_exit
        return
      endif

      call init_coord (start_orb, zero6, gun_ele, upstream_end$)
      call track1 (start_orb, gun_ele, branch%param, end_orb, ignore_radiation = .true.)
      if (.not. particle_is_moving_forward(end_orb)) then
        call out_io (s_fatal$, r_name, 'PARTICLE LOST IN TRACKING E_GUN: ' // gun_ele%name, &
                                       'CANNOT COMPUTE REFERENCE TIME & ENERGY.')
        if (global_com%exit_on_error) call err_exit
        return
      endif

      ! e_gun exit energy gets put into begin_ele exit energy
      pc = (1.0_rp + end_orb%vec(6)) * gun_ele%value(p0c$)
      if (begin_ele%value(p0c$) /= pc) stale = .true.
      begin_ele%value(p0c$) = pc
      call convert_pc_to (begin_ele%value(p0c$), branch%param%particle, e_tot = begin_ele%value(e_tot$))

      ! Now propagate this energy to the e_gun, and any markers in between.
      start_orb%vec(6) = start_orb%p0c * (1.0_rp + start_orb%vec(6)) / begin_ele%value(p0c$) - 1.0_rp
      start_orb%p0c = begin_ele%value(p0c$)

      gun_ele%value(p0c_start$)      = begin_ele%value(p0c$)
      gun_ele%value(e_tot_start$)    = begin_ele%value(e_tot$)
      gun_ele%value(p0c$)            = begin_ele%value(p0c$)
      gun_ele%value(e_tot$)          = begin_ele%value(e_tot$)
      gun_ele%time_ref_orb_in        = start_orb
      gun_ele%value(delta_ref_time$) = end_orb%t - start_orb%t

      do ie = 1, ix_e_gun
        ele => branch%ele(ie)
        ele%value(p0c_start$)     = begin_ele%value(p0c$)
        ele%value(e_tot_start$)   = begin_ele%value(e_tot$)
        ele%value(p0c$)           = begin_ele%value(p0c$)
        ele%value(e_tot$)         = begin_ele%value(e_tot$)
        ele%time_ref_orb_in       = start_orb
      enddo

      if (abs(end_orb%vec(6)) <= bmad_com%rel_tol_tracking) exit
    enddo
  endif

  ! Since Bmad is S-based it cannot handle zero reference energy. 
  ! To avoid roundoff problems set a lower limit of 1d-6 eV.

  if (begin_ele%value(p0c$) < 1d-6) then
    if (ix_e_gun == 0) then
      call out_io (s_fatal$, r_name, 'INITIAL REFERENCE MOMENTUM LESS THAN 1E-6 eV WHICH IS TOO LOW IN BRANCH \i0\ ' // branch%name, &
                                     i_array = [ib])
    else
      call out_io (s_fatal$, r_name, 'INITIAL REFERENCE MOMENTUM LESS THAN 1E-6 eV WHICH IS TOO LOW.', &
                                     'PROBLEM IN PART IS ZERO OR LOW VOLTAGE ON THE E-GUN.')
    endif
    if (global_com%exit_on_error) call err_exit
    return
  endif

  !----------------------------
  ! Loop over all tracking elements in the branch

  ix_super_end = 0  ! End index of current super_lord_region

  do ie = 1, branch%n_ele_track

    ele0 => branch%ele(ie-1)
    ele => branch%ele(ie)

    if (.not. stale .and. ele%bookkeeping_state%ref_energy /= stale$) cycle
    stale = .true.

    !

    if (ele%key == patch$ .or. ele%key == floor_shift$) then
      if (ele0%key == patch$) then
        ele%value(upstream_coord_dir$) = ele0%value(downstream_coord_dir$)
      else
        ele%value(upstream_coord_dir$) = ele0%orientation
      endif

      if (cos(ele%value(x_pitch$)) * cos(ele%value(y_pitch$)) > 0) then
        ele%value(downstream_coord_dir$) = ele%value(upstream_coord_dir$)
      else
        ele%value(downstream_coord_dir$) = -ele%value(upstream_coord_dir$)
      endif
    endif

    ! If we are in the middle of a super_lord region then the "zero" orbit is just the continuation
    ! of the "zero" orbit of the previous element. This is important in wigglers. If the fact
    ! that the "zero" orbit is not truely zero is not taken into account, splitting 
    ! wigglers would result in z-position shifts when tracking particles.

    if (ix_super_end < ie) then       ! If not in super_lord region...
      call init_coord (ele%time_ref_orb_in, zero6, ele0, downstream_end$) ! Note: wrt ele0.
      ele%time_ref_orb_in%location = upstream_end$
    else                              ! In super_lord region
      ele%time_ref_orb_in = ele0%time_ref_orb_out
      ele%time_ref_orb_in%location = upstream_end$
    endif

    ! Find where the current super lord region ends.

    if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
      do j = 1, ele%n_lord
        lord => pointer_to_lord(ele, j)
        if (ele%slave_status == super_slave$ .and. lord%lord_status /= super_lord$) cycle
        lord2 => pointer_to_slave(lord, lord%n_slave) ! last element of lord
        ix_super_end = max(ix_super_end, lord2%ix_ele)
      enddo
    endif

    ! If this element is the first super_slave or multipass_slave of a lord that needs to be auto phased,
    ! make sure that the lord has its RF phase and amplitude properly adjusted.

    if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$ .or. ele%slave_status == multipass_slave$) then 
      do ixl = 1, ele%n_lord
        lord => pointer_to_lord (ele, ixl, ix_slave_back = ix_slave)
        if (ix_slave /= 1) cycle
        if (lord%lord_status /= super_lord$ .and. lord%lord_status /= multipass_lord$) cycle

        nullify(lord2)
        if (lord%slave_status == multipass_slave$) then
          lord2 => lord   ! lord2 is a combination super_lord/multipass_slave element 
          lord => pointer_to_lord(lord, 1)
        endif

        lord_compute = .true.
        if (lord%lord_status == multipass_lord$) then
          call multipass_chain(ele, ix_pass, n_links)
          if (ix_pass /= 1) lord_compute = .false.
        endif

        ! This adjusts the RF phase and amplitude.
        ! Note: Any multipass lord element where the reference energy is not constant must have multipass_ref_energy = first_pass.
        if (lord_compute) then
          lord%time_ref_orb_in = ele%time_ref_orb_in
          if (lord%lord_status == multipass_lord$ .and. nint(lord%value(multipass_ref_energy$)) == user_set$) then
            select case (ele%key)
            case (lcavity$, em_field$, custom$)
              ! These elements store energy values in e_tot_start and p0c_start rather than e_tot and p0c.
              lord%value(p0c$) = lord%value(p0c_start$)
              lord%value(e_tot$) = lord%value(e_tot_start$)
            end select

            if (lord%value(p0c$) > 0) then
              call convert_pc_to(lord%value(p0c$), branch%param%particle, lord%value(e_tot$))
            elseif (lord%value(e_tot$) > 0) then
              call convert_total_energy_to(lord%value(e_tot$), branch%param%particle, pc = lord%value(p0c$))
            else
              ! Marker element, for example, may not set either p0c or e_tot.
              cycle
            endif

            lord%ref_time = 0
            call ele_compute_ref_energy_and_time (lord, lord, branch%param, err)
          else
            call ele_compute_ref_energy_and_time (ele0, lord, branch%param, err)
          endif
          if (err) return
        endif

        ! Compute energy/time for a combo super_lord/multipass_slave element
        if (associated(lord2)) then
          lord2%time_ref_orb_in = ele%time_ref_orb_in
          call ele_compute_ref_energy_and_time (ele0, lord2, branch%param, err)
          if (err) return
        endif
      enddo
    endif

    ! Calculate the energy and reference time at the end of the present element.

    call ele_compute_ref_energy_and_time (ele0, ele, branch%param, err)
    if (err) then
      call out_io (s_error$, r_name, 'CANNOT COMPUTE REFERENCE ENERGY FOR: ' // ele%name, &
                                      'REFERENCE ENERGY NOT COMPUTED FOR REST OF LATTICE')
      return
    endif

    ele%bookkeeping_state%ref_energy = ok$

  enddo

  ! Need to do this at the end of the loop since if it is at the beginning the set may be overwritten.

  branch%param%bookkeeping_state%ref_energy = ok$

enddo ! Branch loop

! Correct particle_start info

!call init_coord (lat%particle_start, lat%particle_start, lat%ele(0), downstream_end$, &
!                                                            E_photon = lat%particle_start%p0c, shift_vec6 = .false.)

! Put the appropriate energy values in the lord elements...

do ie = lat%n_ele_track+1, lat%n_ele_max

  lord => lat%ele(ie)

  if (lord%bookkeeping_state%ref_energy /= stale$) cycle
  lord%bookkeeping_state%ref_energy = ok$

  if (lord%n_slave == 0) cycle   ! Can happen with null_ele$ elements for example.

  ! Multipass lords have their enegy computed above.

  if (lord%lord_status == multipass_lord$) cycle
  if (lord%lord_status == ramper_lord$) cycle

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
  lord%ref_species   = slave%ref_species

  ! Transfer the starting energy.

  if (lord%lord_status == super_lord$) then
    slave => pointer_to_slave(lord, 1)
    lord%value(E_tot_start$) = slave%value(E_tot_start$)
    lord%value(p0c_start$)   = slave%value(p0c_start$)
  endif

  abs_tol(1:2) = 1d-3 + bmad_com%rel_tol_adaptive_tracking * lord%value(p0c_start$)
  abs_tol(3)   = 1d-3 + bmad_com%rel_tol_adaptive_tracking * lord%value(p0c$)

  if (ele_value_has_changed(lord, [p0c$, e_tot$, delta_ref_time$], abs_tol, .true.)) then
    call set_ele_status_stale (lord, attribute_group$)
  endif

enddo ! Branch loop

lat%lord_state%ref_energy = ok$
err_flag = .false.

end subroutine lat_compute_ref_energy_and_time

!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!------------------------------------------------------------------------------------------
!+
! Subroutine ele_compute_ref_energy_and_time (ele0, ele, param, err_flag)
!
! Routine to compute the reference energy and reference time at the end of an element 
! given the reference enegy and reference time at the start of the element.
!
! Input:
!   ele0           -- ele_struct: Previous element in lattice with starting energy and time values.
!   ele            -- Ele_struct: Lattice element
!     %time_ref_orb_in  -- Starting orbit for ref time calc.
!   param          -- lat_Param_struct: Lattice parameters.
!   err_flag       -- Logical: Set true if there is an error. False otherwise.
!
! Output:
!   ele         -- Ele_struct: Lattice element with reference energy and time.
!     %time_ref_orb_out  -- Ending orbit for ref time calc.
!-

subroutine ele_compute_ref_energy_and_time (ele0, ele, param, err_flag)

use bmad_interface, dummy => ele_compute_ref_energy_and_time
use autoscale_mod, only: autoscale_phase_and_amp
use radiation_mod, only: track1_radiation

implicit none

type (ele_struct), target :: ele, ele0
type (ele_struct), pointer :: lord
type (branch_struct), pointer :: branch
type (floor_position_struct) old_floor
type (lat_param_struct) :: param
type (coord_struct) orb_start, orb_end, orb1, orb2
type (bmad_common_struct) bmad_com_saved

real(rp) E_tot_start, p0c_start, ref_time_start, e_tot, p0c, phase, velocity, abs_tol(3)
real(rp) value_saved(num_ele_attrib$), beta0, ele_ref_time

integer i, key
logical err_flag, err, changed, saved_is_on, energy_stale, do_track

character(*), parameter :: r_name = 'ele_compute_ref_energy_and_time'

! bmad_com%auto_bookkeeper is set False to prevent track1 calling attribute_bookkeeper 
! which will overwrite ele%old_value.

err_flag = .true.

bmad_com_saved = bmad_com
bmad_com%radiation_fluctuations_on = .false.
bmad_com%radiation_damping_on = .false.
bmad_com%radiation_zero_average = .false.
bmad_com%auto_bookkeeper = .false.

E_tot_start    = ele0%value(E_tot$)
p0c_start      = ele0%value(p0c$)
ref_time_start = ele0%ref_time
value_saved = ele%value

ele%value(E_tot_start$)    = E_tot_start
ele%value(p0c_start$)      = p0c_start
ele%value(ref_time_start$) = ref_time_start

ele%time_ref_orb_out = ele%time_ref_orb_in  ! This should be true when we don't have to track.
ele%time_ref_orb_out%location = downstream_end$

key = ele%key
if (key == em_field$ .and. is_false(ele%value(constant_ref_energy$))) key = lcavity$

if (key == converter$) then
  ele%ref_species = ele%converter%species_out
else
  ele%ref_species = ele0%ref_species
endif

!

select case (key)

case (converter$)
  ele%ref_time = ref_time_start + ele%value(l$) * E_tot_start / (p0c_start * c_light)

  if (ele%value(p0c$) == 0 .and. ele%value(E_tot$) == 0)  then
    call out_io (s_abort$, r_name, 'NEITHER P0C NOR E_TOT IS SET FOR CONVERTER ' // ele%name)
    goto 9000
  elseif (ele%value(p0c$) == 0) then
    call convert_total_energy_to (ele%value(E_tot$), ele%converter%species_out, pc = ele%value(p0c$))
  else
    call convert_pc_to (ele%value(p0c$), ele%converter%species_out, E_tot = ele%value(E_tot$))
  endif

case (lcavity$)

  ! A zero e_tot (Can happen on first pass through this routine) can mess up tracking so put 
  ! in a temp value if needed. This does not affect the phase & amp adjustment.
  if (ele%value(e_tot$) == 0) then
    ele%value(E_tot$) = ele%value(e_tot_start$)      
    ele%value(p0c$) = ele%value(p0c_start$)
    ele%ref_time = ref_time_start
    call attribute_bookkeeper(ele, .true.)
  endif

  if (ele%slave_status /= super_slave$ .and. ele%slave_status /= slice_slave$ .and. ele%slave_status /= multipass_slave$) then
    call autoscale_phase_and_amp (ele, param, err, call_bookkeeper = .false.)
    if (err) goto 9000
  endif

  ! If ele is super_slave and the super_lord and ele are essentially the same element, do not need to track.

  do_track = .true.
  if (ele%slave_status == super_slave$) then
    do i = 1, ele%n_lord
      lord => pointer_to_lord(ele, i)
      if (lord%key /= lcavity$) cycle  ! May be a pipe$
      if (lord%n_slave /= 1) cycle
      do_track = .false.
      ele%value(p0c$) = lord%value(p0c$)
      ele%value(E_tot$) = lord%value(E_tot$)
      ele%value(delta_ref_time$) = lord%value(delta_ref_time$)
      ele%ref_time = lord%ref_time
      call calc_time_ref_orb_out(lord%time_ref_orb_out)
    enddo
  endif

  ! Track. With runge_kutta (esp fixed step with only a few steps), a shift in the end energy can cause 
  ! small changes in the tracking. So if there has been a shift in the end energy, track again.
  

  if (do_track) then
    do i = 1, 5
      call track_this_ele (orb_start, orb_end, ref_time_start, .false., err); if (err) goto 9000
      call calc_time_ref_orb_out(orb_end)
      ele%value(p0c$) = orb_end%p0c * (1 + orb_end%vec(6))
      call convert_pc_to (ele%value(p0c$), param%particle, E_tot = ele%value(E_tot$), err_flag = err)
      if (err) goto 9000
      if (abs(orb_end%vec(6)) < small_rel_change$) exit
    enddo
  endif

  ele%value(delta_ref_time$) = ele%ref_time - ref_time_start
  ele%time_ref_orb_out%p0c = ele%value(p0c$)  ! To prevent small roundoff errors

case (custom$, hybrid$)
  ele%value(E_tot$) = E_tot_start + ele%value(delta_e_ref$)   ! Delta_ref_time is an independent attrib here.
  call convert_total_energy_to (ele%value(E_tot$), param%particle, pc = ele%value(p0c$), err_flag = err)
  if (err) goto 9000

  ele%ref_time = ref_time_start + ele%value(delta_ref_time$)

case (taylor$)
  ele%value(E_tot$) = E_tot_start + ele%value(delta_e_ref$)
  call convert_total_energy_to (ele%value(E_tot$), param%particle, pc = ele%value(p0c$), err_flag = err)
  if (err) goto 9000

  ! Delta_ref_time is an independent attrib if L = 0.
  if (ele%value(l$) /= 0) then
    ele%value(delta_ref_time$) = ele%value(l$) * E_tot_start / (p0c_start * c_light)
  endif
  ele%ref_time = ref_time_start + ele%value(delta_ref_time$)

case (e_gun$)
  ! Note: Due to the coupling between an e_gun and the begin_ele, autoscaling is
  ! done in lat_compute_ref_energy_and_time.
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start

  branch => pointer_to_branch(ele)
  if (associated(branch)) then
    if (ele%s_start == branch%ele(0)%s) then
      orb_start%vec(6) = branch%ele(0)%value(p0c_start$) / branch%ele(0)%value(p0c$) - 1
    endif
  endif

  call track_this_ele (orb_start, orb_end, ref_time_start, .true., err); if (err) goto 9000
  call calc_time_ref_orb_out(orb_end)
  ele%value(delta_ref_time$) = ele%ref_time - ref_time_start

case (crystal$, mirror$, multilayer_mirror$, diffraction_plate$, photon_init$, mask$)
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start

  ele%value(ref_wavelength$) = c_light * h_planck / E_tot_start
  ele%value(delta_ref_time$) = 0
  ele%ref_time = ref_time_start

case (patch$)
  ! Ignore flexible patch error messages that can be generated.
  if (is_true(ele%value(flexible$))) then
    call ele_geometry(ele0%floor, ele, ignore_patch_err = .true.)
  endif

  if (ele%is_on .and. ele%value(e_tot_offset$) /= 0) then
    if (ele%orientation == 1) then
      ele%value(E_tot$) = e_tot_start + ele%value(e_tot_offset$)
    else
      ele%value(E_tot$) = e_tot_start - ele%value(e_tot_offset$)
    endif
    call convert_total_energy_to (ele%value(E_tot$), param%particle, pc = ele%value(p0c$), err_flag = err)
    if (err) goto 9000

  elseif (ele%value(E_tot_set$) /= 0) then
    ele%value(E_tot$) = ele%value(E_tot_set$)
    call convert_total_energy_to (ele%value(E_tot_set$), param%particle, pc = ele%value(p0c$))

  elseif (ele%value(p0c_set$) /= 0) then
    ele%value(p0c$) = ele%value(p0c_set$)
    call convert_pc_to (ele%value(p0c$), param%particle, E_tot = ele%value(E_tot$))

  else
    ele%value(E_tot$) = E_tot_start
    ele%value(p0c$) = p0c_start
  endif

  velocity = c_light * ele%value(p0c$) / ele%value(E_tot$)
  ele%value(delta_ref_time$) = ele%value(t_offset$) + ele%value(l$) / velocity
  ele%ref_time = ref_time_start + ele%value(delta_ref_time$)

case (marker$, fork$, photon_fork$)
  ele%value(E_tot$) = E_tot_start
  ele%value(p0c$) = p0c_start
  ele%value(delta_ref_time$) = 0
  ele%ref_time = ref_time_start

case default
  if (significant_difference(ele%value(p0c$), p0c_start, rel_tol = small_rel_change$)) then
    ele%value(E_tot$) = E_tot_start
    ele%value(p0c$) = p0c_start
    ! Need to call attribute_bookkeeper since num_steps is not set until the energy is set.
    call attribute_bookkeeper(ele, .true.) 
  endif
  
  if (ele%key == rfcavity$ .and. ele%slave_status /= super_slave$ .and. &
                        ele%slave_status /= slice_slave$ .and. ele%slave_status /= multipass_slave$) then
    call autoscale_phase_and_amp (ele, param, err, call_bookkeeper = .false.)
    if (err) goto 9000
  endif

  if (ele_has_constant_ds_dt_ref(ele)) then
    if (ele%value(l$) == 0) then     ! Must avoid problem with zero length markers and p0c = 0.
      ele%value(delta_ref_time$) = 0
    else
      ele%value(delta_ref_time$) = ele%value(l$) * E_tot_start / (p0c_start * c_light)
    endif
    ele%ref_time = ref_time_start + ele%value(delta_ref_time$)

  else
    call track_this_ele (orb_start, orb_end, ref_time_start, .false., err); if (err) goto 9000
    call calc_time_ref_orb_out(orb_end)
    ele%value(delta_ref_time$) = ele%ref_time - ref_time_start
  endif
end select

! If delta_ref_time has shifted then any taylor map must be updated.

ele%old_value = value_saved

if (abs(ele%value(delta_ref_time$) - ele%old_value(delta_ref_time$)) * c_light > bmad_com%significant_length) then
  if (associated (ele%taylor(1)%term) .and. ele%key /= taylor$) then
    call kill_taylor (ele%taylor)
    call kill_taylor (ele%spin_taylor)
  endif
  ele%bookkeeping_state%mat6 = stale$
endif

! Changes in the reference energy should trigger bookkeeping. 
! Example: A lattice with bends with field_master = True, flexible patch, and absolute time tracking has
! the reference energy dependent upon the geometry and the geometry depends upon the ref energy.

abs_tol(1:2) = 1d-3 + bmad_com%rel_tol_adaptive_tracking * ele%value(p0c$)
abs_tol(3) = bmad_com%significant_length/c_light

energy_stale = ele_value_has_changed(ele, [p0c$, e_tot$, delta_ref_time$], abs_tol, .false.)
if (energy_stale.or. ele%bookkeeping_state%control /= ok$ .or. ele%bookkeeping_state%floor_position /= ok$) then
  ! Transfer ref energy to super_lord before bookkeeping done. This is important for bends.
  if (ele%slave_status == super_slave$ .and. (ele%key == sbend$ .or. ele%key == rf_bend$)) then
    do i = 1, ele%n_lord
      lord => pointer_to_lord(ele, i)
      lord%value(e_tot$) = ele%value(e_tot$); lord%value(e_tot_start$) = ele%value(e_tot_start$)
      lord%value(p0c$) = ele%value(p0c$); lord%value(p0c_start$) = ele%value(p0c_start$)
      call control_bookkeeper (ele%branch%lat, lord)
    enddo
  endif

  call control_bookkeeper (ele%branch%lat, ele)
  old_floor = ele%floor
  call ele_geometry (ele0%floor, ele, ignore_patch_err = .true.)
  ele%bookkeeping_state%floor_position = ele0%bookkeeping_state%floor_position

  if (ele%ix_ele > 0) then ! If element is associated with a lattice
    if (.not. old_floor == ele%floor) call set_lords_status_stale (ele, floor_position_group$) ! Need to update girders
    call set_lords_status_stale (ele, ref_energy_group$)!, control_bookkeeping = .true.)
  endif
endif

!

err_flag = .false.

9000 continue
bmad_com = bmad_com_saved

!---------------------------------------------------------------------------------
contains

subroutine track_this_ele (orb_start, orb_end, ref_time_start, is_inside, error)

type (coord_struct) orb_start, orb_end
type (ele_struct), pointer :: lord
real(rp) ref_time_start
integer ie
logical error, is_inside, totalpath_saved

! With use_totalpath = False (the standard setting) and with ptc tracking, orb_end%t is calculated
! using ele%ref_time. But this is what we want to calculate here.

error = .true.
call set_ptc_base_state('TOTALPATH', .true., totalpath_saved)

call zero_errors_in_ele (ele, changed)
call init_coord (orb_start, ele%time_ref_orb_in, ele, upstream_end$, shift_vec6 = .false.)
if (is_inside) orb_start%location = inside$ ! To avoid entrance kick in time RK tracking
call track1 (orb_start, ele, param, orb_end, ignore_radiation = .true.)
if (.not. particle_is_moving_forward(orb_end)) then
  call out_io (s_fatal$, r_name, 'PARTICLE LOST IN TRACKING: ' // ele%name, &
                                 'CANNOT COMPUTE REFERENCE TIME & ENERGY.')
  if (global_com%exit_on_error) call err_exit
  return
endif
call restore_errors_in_ele (ele)

! Note: delta_ref_time is used in runge_kutta tracking.
! Note: There is a potential problem with RF elements with RK tracking and slicing and where 
! the particle is non-relativistic. To avoid z-shifts with slicing, use the lord delta_ref_time.

if (ele%tracking_method /= bmad_standard$ .and. &
          (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$)) then
  do ie = 1, ele%n_lord
    lord => pointer_to_lord(ele, ie)
    if (lord%key /= pipe$) exit
  enddo

  if (lord%value(l$) == 0) then
    ele%value(delta_ref_time$) = 0
  else
    ele%value(delta_ref_time$) = lord%value(delta_ref_time$) * ele%value(l$) / lord%value(l$)
  endif
else
  ele%value(delta_ref_time$) = orb_end%t - orb_start%t
endif

ele%ref_time = ref_time_start + ele%value(delta_ref_time$)

call set_ptc_base_state('TOTALPATH', totalpath_saved)
error = .false.

end subroutine track_this_ele

!---------------------------------------------------------------------------------
! contains

recursive subroutine zero_errors_in_ele (ele, changed)

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer i
logical changed, has_changed

! For reference energy tracking need to turn off any element offsets and kicks and zero any errors.
! If the element is a super_slave or multipass_slave then the errors must be zeroed in the lord elements also.

ele%old_value = ele%value
saved_is_on = ele%is_on
ele%is_on = .true.
ele%iyy = ele%tracking_method

has_changed = .false.

if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$ .or. ele%slave_status == multipass_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status == girder_lord$ .or. lord%lord_status == overlay_lord$ .or. lord%lord_status == group_lord$) cycle
    call zero_errors_in_ele (lord, changed)
    if (changed) has_changed = .true.
  enddo
endif

if (ele_has_nonzero_offset(ele)) then
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

  if (ele%value(gradient_err$) /= 0 .or. ele%value(voltage_err$) /= 0) then
    ele%value(gradient_err$) = 0
    ele%value(voltage_err$)  = 0
    ele%value(gradient_tot$) = ele%value(gradient$)
    ele%value(voltage_tot$)  = ele%value(voltage$)
    has_changed = .true.
  endif

case (e_gun$)
  if (ele%value(phi0_err$) /= 0) then
    ele%value(phi0_err$) = 0
    has_changed = .true.
  endif

  if (ele%value(gradient_err$) /= 0 .or. ele%value(voltage_err$) /= 0) then
    ele%value(gradient_err$) = 0
    ele%value(voltage_err$)  = 0
    ele%value(gradient_tot$) = ele%value(gradient$)
    ele%value(voltage_tot$)  = ele%value(voltage$)
    has_changed = .true.
  endif

case (rfcavity$)
  if (bmad_com%rf_phase_below_transition_ref) then
    if (ele%value(phi0$) /= 0.5_rp) then
      ele%value(phi0$) = 0.5_rp
      has_changed = .true.
    endif    
  else
    if (ele%value(phi0$) /= 0) then
      ele%value(phi0$) = 0
      has_changed = .true.
    endif    
  endif
end select

! For speed, use symp_lie_ptc tracking if the taylor map does not exist or if the taylor
! map is not valid for the element with kicks and offsets removed.

changed = has_changed
if (ele%tracking_method == taylor$ .and. .not. associated (ele%taylor(1)%term)) ele%tracking_method = symp_lie_ptc$

end subroutine zero_errors_in_ele

!---------------------------------------------------------------------------------
! contains

recursive subroutine restore_errors_in_ele (ele)

type (ele_struct) ele
type (ele_struct), pointer :: lord
integer i

!

ele%value = ele%old_value
ele%is_on = saved_is_on
ele%tracking_method = ele%iyy

if (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$  .or. &
                                                                 ele%slave_status == multipass_slave$) then
  do i = 1, ele%n_lord
    lord => pointer_to_lord(ele, i)
    if (lord%lord_status == girder_lord$ .or. lord%lord_status == overlay_lord$ .or. lord%lord_status == group_lord$) cycle
    call restore_errors_in_ele (lord)
  enddo
endif

end subroutine restore_errors_in_ele

!---------------------------------------------------------------------------------
! contains

subroutine calc_time_ref_orb_out (orb_end)

type (coord_struct) orb_end
type (ele_struct), pointer :: lord

!

ele%time_ref_orb_out = orb_end
ele%time_ref_orb_out%location = downstream_end$

! The tracking did not have the correct delta_ref_time end exit end ref energy so need to correct for this.
! Notice that here the delta_ref_time value (but not ele%value(p0c$)) the value used in tracking and not the 
! corrected value computed later.

! Note: zero length slice of e_gun can have orb_end%vec(6) = -1
if (orb_end%vec(6) /= -1) then
  ele%time_ref_orb_out%vec(2) = ele%time_ref_orb_out%vec(2) / (1 + orb_end%vec(6))
  ele%time_ref_orb_out%vec(4) = ele%time_ref_orb_out%vec(4) / (1 + orb_end%vec(6))
endif

if (ele%key == lcavity$ .and. ele%tracking_method == bmad_standard$) then
  ele%time_ref_orb_out%vec(5) = 0
elseif (ele%slave_status == super_slave$ .or. ele%slave_status == slice_slave$) then
  lord => pointer_to_lord(ele, 1)
  beta0 = lord%value(l$) / (c_light * lord%value(delta_ref_time$))
  ele_ref_time = ele%value(ref_time_start$) + ele%value(l$) / (c_light * beta0)
  ele%time_ref_orb_out%vec(5) = c_light * orb_end%beta * (ele_ref_time - orb_end%t)
else
  ele%time_ref_orb_out%vec(5) = c_light * orb_end%beta * (ele%ref_time - orb_end%t)
endif

!ele%time_ref_orb_out%vec(5) = ele%time_ref_orb_out%vec(5) + &
!            (orb_end%t - orb_start%t - ele%value(delta_ref_time$)) * orb_end%beta * c_light
ele%time_ref_orb_out%vec(6) = ((1 + orb_end%vec(6)) * orb_end%p0c - ele%value(p0c$)) / ele%value(p0c$)

end subroutine

end subroutine ele_compute_ref_energy_and_time
