program bbu_program

use bbu_track_mod

implicit none

type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in, lat0
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele, ele2

integer i, ix, j, n_hom, n, n_ele, ix_pass, i_lr

integer istep
integer :: nstep = 100

integer irep

real(rp) deldr,dr, charge_try
real(rp) hom_power_gain, charge0, charge1, charge_old, charge_threshold
real(rp) trtb, currth, growth_rate, growth_rate0, growth_rate1, growth_rate_old
real(rp) time

logical lost
logical, allocatable :: keep_ele(:)

namelist / bbu_params / bbu_param, beam_init 

! Read in parameters

bbu_param%lat_file_name = 'erl.lat'  ! Bmad lattice file name
bbu_param%simulation_turns_max = 20 ! 
bbu_param%bunch_freq = 1.3e9         ! Freq in Hz.
bbu_param%init_particle_offset = 1e-8  ! Initial particle offset for particles born 
                                       !  in the first turn period.
bbu_param%limit_factor = 2           ! Init_hom_amp * limit_factor = simulation unstable limit
bbu_param%hybridize = .true.         ! Combine non-hom elements to speed up simulation?
bbu_param%keep_overlays_and_groups = .false. ! keep when hybridizing?
bbu_param%keep_all_lcavities       = .false. ! keep when hybridizing?
bbu_param%current = 20e-3            ! Starting current (amps)
bbu_param%rel_tol = 1e-2             ! Final threshold current accuracy.
bbu_param%write_hom_info = .true.  
bbu_param%drscan = .true.        ! If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3.
bbu_param%nstep = 100
bbu_param%begdr = 5.234
bbu_param%enddr = 6.135
bbu_param%use_interpolated_threshold = .true.
bbu_param%nrep = 1     ! Number of times to repeat threshold calculation
bbu_param%ran_seed = 0
bbu_param%ran_gauss_sigma_cut = -1
bbu_param%stable_orbit_anal = .false.

beam_init%n_particle = 1

open (1, file = 'bbu.init', status = 'old')
read (1, nml = bbu_params)
close (1)

if (bbu_param%stable_orbit_anal) bbu_param%nstep = 1

! Define distance between bunches

beam_init%dt_bunch = 1 / bbu_param%bunch_freq
call ran_seed_put (bbu_param%ran_seed)
print *, 'Random number seed:', bbu_param%ran_seed

if (bbu_param%ran_gauss_sigma_cut > 0) then
  call ran_gauss_converter (set_sigma_cut = bbu_param%ran_gauss_sigma_cut)
  print *, 'ran_gauss sigma cut: ', bbu_param%ran_gauss_sigma_cut 
endif

! Init

print *, 'Lattice file: ', trim(bbu_param%lat_file_name)
call bmad_parser (bbu_param%lat_file_name, lat_in)
call twiss_propagate_all (lat_in)
call run_timer ('START')

! Initialize DR scan

nstep = 1
dr=0
if (bbu_param%drscan) then
! If DRSCAN requested, set number of calculation repetitions to 1
  bbu_param%nrep = 1
! Determine element index of variable-length element
  do i = 1, lat_in%n_ele_max
    if (lat_in%ele(i)%name .eq. bbu_param%elname)then
      bbu_param%elindex = i
      write(6,'(a,i6,a,a40)')&
            ' Element index ',i,' found for DR scan element ',bbu_param%elname
! Open DRSCAN output file for threshold current calculation comparison
      open (50, file = 'drscan.out', status = 'unknown') 
      nstep = bbu_param%nstep
      exit  ! Exit loop over elements
    else
      cycle ! Continue loop over elements
    endif
    print *,' No DR scan element found, disabling scan'
    bbu_param%drscan = .false.
  enddo
endif

! Open file for storing output of repeated threshold calculations
if (bbu_param%nrep.gt.1)then
  write(6,'(a,i10,a)')&
        ' Opening output file for',bbu_param%nrep,' repetitions'
  open (55, file = 'rep.out', status = 'unknown')
endif

! Open file for stable orbit analysis data
if (bbu_param%stable_orbit_anal.eq..true.)then
  write(6,'(a,i10,a)')&
        ' Opening output file for stable orbit analysis'
  open (56, file = 'stable_orbit.out', status = 'unknown')
! Write number of repetitions and nr wake elements
!  write(56,'(2i10)'),bbu_param%nrep,size(bbu_beam%stage)
endif

!-------------------------------------------------------
! DRSCAN loop. NSTEP=1 if no DRSCAN.

do istep = 1, nstep

  if (bbu_param%drscan) then
    ! Change length of taylor element
    deldr = 0.0
    if(nstep > 1) deldr = (bbu_param%enddr - bbu_param%begdr)/(nstep-1)
    dr = bbu_param%begdr + (istep-1) * deldr
    lat_in%ele(bbu_param%elindex)%value(l$) = dr * c_light / bbu_param%bunch_freq
    write(6,'(a,2f8.3)')' DRSCAN analysis step: dr, scan element length = ', &
                 dr, lat_in%ele(bbu_param%elindex)%value(l$)
    call lattice_bookkeeper(lat_in)
  endif

  if (bbu_param%hybridize) then
    print *, 'Note: Hybridizing lattice...'
    allocate (keep_ele(lat_in%n_ele_max))
    keep_ele = .false.
    do i = 1, lat_in%n_ele_max
      if (lat_in%ele(i)%key /= lcavity$) cycle
      if (.not. bbu_param%keep_all_lcavities) then
        if (.not. associated (lat_in%ele(i)%wake)) cycle
        if (size(lat_in%ele(i)%wake%lr) == 0) cycle
      endif
      call update_hybrid_list (lat_in, i, keep_ele, bbu_param%keep_overlays_and_groups)
    enddo
    call make_hybrid_lat (lat_in, keep_ele, .true., lat)
    deallocate (keep_ele)
  else
    lat = lat_in
  endif

  lat0 = lat

  ! Print some information and get the analytic approximation result for the threshold current

  print '(2a)', ' Lattice File: ', trim(lat%input_file_name)
  print '(2a)', ' Lattice Name: ', trim(lat%lattice)
  print '(a, i7)', ' Nr Tracking Elements: ', lat%n_ele_track
  print '(a, es12.2)', ' Beam Energy: ', lat%ele(0)%value(e_tot$)

  if (bbu_param%write_hom_info) then
   call write_homs(lat, bbu_param%bunch_freq, trtb, currth)
  ! Update starting current according to analytic approximation
   if (currth.gt.0.)bbu_param%current = currth
  else
   trtb=dr
  endif

  call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)

  print *, 'Number of lr wake elements in tracking lattice:', size(bbu_beam%stage)

  n_ele = 0
  do i = 1, size(bbu_beam%stage)
    j = bbu_beam%stage(i)%ix_ele_lr_wake
    call multipass_chain (lat%ele(j), lat, ix_pass, n_links = n)
    if (ix_pass /= 1 .and. n /= 0) cycle
    n_ele = n_ele + 1
  enddo

  print *, 'Number of physical lr wake elements:', n_ele
  print *, 'Number of elements in lattice:      ', lat%n_ele_track

  ! Loop over calculation repetitions
  do irep = 1,bbu_param%nrep

     beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

    if (bbu_param%stable_orbit_anal) then

      print *,' Analyzing stable orbit for repetition ', irep

      lat = lat0 ! Restore lr wakes
      call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_power_gain, growth_rate, lost)

    ! Print output for stable orbit analysis
      do i = 1, size(bbu_beam%stage)
        write(56,'(i10,2(1x,e15.6),i10/,4(6(1x,e15.6)/))')i,bbu_beam%stage(i)%time_at_wake_ele, &
                                               bbu_beam%stage(i)%hom_power_max,bbu_beam%stage(i)%n_orb, & 
                      bbu_beam%stage(i)%ave_orb, bbu_beam%stage(i)%rms_orb, &
                      bbu_beam%stage(i)%min_orb, bbu_beam%stage(i)%max_orb

      enddo

    ! Re-randomize HOM frequencies

      do i = 1, lat0%n_ele_max
        call randomize_lr_wake_frequencies (lat0%ele(i))
      enddo
! Skip threshold calculation
      cycle

    endif

    charge0 = 0
    charge1 = -1      ! Mark as not set yet 
    charge_old = -1   ! Mark as not set yet 
    charge_try = -1   ! Mark as not set yet 

    Print *, 'Searching for a current where the tracking is unstable...'

    do
      lat = lat0 ! Restore lr wakes
      call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_power_gain, growth_rate, lost)
      call calc_next_charge_try
      if (hom_power_gain > 1) exit
      if (lost) then
        print *, 'Particle(s) lost. Assuming unstable...'
        exit
      endif
    enddo

     ! Track to bracket threshold

      print *, 'Now converging on the threshold...'

      do
        lat = lat0 ! Restore lr wakes
        call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_power_gain, growth_rate, lost)
        if (lost) print *, 'Particle(s) lost. Assuming unstable...'
        call calc_next_charge_try
        if (charge1 - charge0 < charge1 * bbu_param%rel_tol) exit
      enddo

      beam_init%bunch_charge = (charge0 + charge1) / 2
      print *, 'Threshold Current (A):', beam_init%bunch_charge / beam_init%dt_bunch 
      i = bbu_beam%ix_stage_power_max
      j = bbu_beam%stage(i)%ix_ele_lr_wake
      ele => lat%ele(j)
      ele2 => ele
      if (ele2%slave_status == multipass_slave$) then
        ele2 => pointer_to_multipass_lord (ele, lat)
      endif
      i_lr = bbu_beam%stage(i)%ix_hom_max
      print *, 'Element with critical HOM:', ele2%ix_ele, ':   ', ele2%name
      print *, 'Critical HOM: Input Frequency: ', ele%wake%lr(i_lr)%freq_in 
      print *, 'Critical HOM: Actual Frequency:', ele%wake%lr(i_lr)%freq
      print *, 'Critical HOM: R_overQ:         ', ele%wake%lr(i_lr)%r_over_q
      print *, 'Critical HOM: Q:               ', ele%wake%lr(i_lr)%q
      print *, 'Critical HOM: Angle:           ', ele%wake%lr(i_lr)%angle

      if (bbu_param%nrep.gt.1)write(55,'(i6,e14.6,2i7,6(e14.6,1x))') &
          irep, beam_init%bunch_charge / beam_init%dt_bunch, ele%ix_ele, ele2%ix_ele, ele%s, &
          ele%wake%lr(i_lr)%freq_in, ele%wake%lr(i_lr)%freq, ele%wake%lr(i_lr)%r_over_q, &
          ele%wake%lr(i_lr)%q, ele%wake%lr(i_lr)%angle

    ! Re-randomize HOM frequencies

    do i = 1, lat0%n_ele_max
      call randomize_lr_wake_frequencies (lat0%ele(i))
    enddo

  enddo  ! End of repetition loop

  if (bbu_param%drscan) write(50,*) trtb, currth, beam_init%bunch_charge / beam_init%dt_bunch 

enddo  ! End of DRSCAN loop

if (bbu_param%drscan) close(50)
if (bbu_param%nrep.gt.1)close(55)
if (bbu_param%stable_orbit_anal.eq..true.)close(56)

call run_timer ('STOP', time)
print *
print *, 'Time for calculation (min): ', time/60

!-------------------------------------------------------------------------------
contains

subroutine calc_next_charge_try()

real(rp) c, min_delta, dc0, dc1, c0, c1, g0, g1

! Print info

if (growth_rate > 0 .or. lost) then
  print *, '  Unstable at (mA):', 1e3 * beam_init%bunch_charge / beam_init%dt_bunch 
else
  print *, '  Stable at (mA):', 1e3 * beam_init%bunch_charge / beam_init%dt_bunch 
endif

print *, '         Head bunch index: ', bbu_beam%bunch(bbu_beam%ix_bunch_head)%ix_bunch
print *, '         Growth rate: ', growth_rate

!

if (growth_rate > 0 .or. lost) then  ! unstable
  charge1 = beam_init%bunch_charge
  growth_rate1 = growth_rate

else
  charge0 = beam_init%bunch_charge
  growth_rate0 = growth_rate
endif


! Cannot print threshold prediction until there are two trackings.

if (charge_old > 0) then  ! If we have two trackings.

  !---------------------------------------------------
  ! If Have bracked the threshold

  if (charge0 > 0 .and. charge1 > 0) then  

    c0 = charge0;  g0 = growth_rate0
    c1 = charge1;  g1 = growth_rate1

    ! the tracking may not have calculated a valid growth rate. 
    ! If so, use the average current as the threshold

    if (g0 == real_garbage$ .or. g1 == real_garbage$ .or. g0 == g1) then
      charge_threshold = (c0 + c1) / 2
    else
      charge_threshold = (c0 * g1 - c1 * g0) / (g1 - g0)
    endif

    ! current to use in the next tracking must be significantly different from c0 and c1.

    charge_try = charge_threshold

    min_delta = max(c0, c1) * bbu_param%rel_tol
    dc0 = charge_threshold - c0
    dc1 = charge_threshold - c1

    ! If closer to c0...
    if (abs(dc0) < abs(dc1) .and. abs(dc0) < min_delta) then
      c = charge_threshold + sign(min_delta, dc0)
      if (abs(c-c1) < abs(c-c0)) then ! If moved closer to c1 then just use the average
        charge_try = (c1 + c0) / 2
      else
        charge_try = c
      endif

    ! If closer to c1...
    else if (dc1 < dc0 .and. dc1 < min_delta) then
      c = charge_threshold + sign(min_delta, dc1)
      if (abs(c-c0) < abs(c-c1)) then  ! If moved closer to c0 then just use the average
        charge_try = (c1 + c0) / 2
      else
        charge_try = c
      endif
    endif

  !---------------------------------------------------
  ! If Have stable points but no unstable point

  elseif (charge1 < 0) then  
    c0 = charge0;  g0 = growth_rate0
    c1 = charge_old;  g1 = growth_rate_old

    ! the tracking may not have calculated a valid growth rate. 
    ! If so, use twice the last current as the threshold

    if (g0 == real_garbage$ .or. g1 == real_garbage$ .or. g0 == g1) then
      charge_threshold = 2 * c0   
    else
      charge_threshold = (c0 * g1 - c1 * g0) / (g1 - g0)
    endif

    ! If the threshold is less than charge0 then there must be a lot of noise so
    ! assume we are far from the threshold and increase the current by a factor of 10.
    ! In any case, Demand at least a 10% change.

    if (charge_threshold < c0) then
      charge_try = 10 * c0
    else
      charge_try = max(charge_threshold, 1.1 * c0)
    endif

  !---------------------------------------------------
  ! If Have unstable points but no stable point

  else                       
    c0 = charge_old;  g0 = growth_rate_old
    c1 = charge1;  g1 = growth_rate1

    ! the tracking may not have calculated a valid growth rate. 
    ! If so, use half the last current as the threshold

    if (g0 == real_garbage$ .or. g1 == real_garbage$ .or. g0 == g1) then
      charge_threshold = c0 / 2   
    else
      charge_threshold = (c0 * g1 - c1 * g0) / (g1 - g0)
    endif

    ! Demand at least a 10% change but if negative just set to 0.

    charge_try = min(charge_threshold, 0.9 * c1)
    if (charge_try < 0) charge_try = 0
  endif


  print *, '         Predicted threshold:', 1d3 * charge_threshold / beam_init%dt_bunch 

endif

charge_old = beam_init%bunch_charge
growth_rate_old = growth_rate

! Set new try value

if (bbu_param%use_interpolated_threshold .and. charge_try >= 0) then
  beam_init%bunch_charge = charge_try

else
  if (charge1 > 0) then  ! Has been set so have stable and unstable points
    beam_init%bunch_charge = (charge0 + charge1) / 2
  else                   ! Still searching for an unstable point.
    beam_init%bunch_charge = beam_init%bunch_charge * 2
  endif
endif

print *, '         Current to try next:', 1d3 * beam_init%bunch_charge / beam_init%dt_bunch 

end subroutine

end program
