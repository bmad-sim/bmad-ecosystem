program bbu_program

use bbu_track_mod_old

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
integer n_loc
logical err

type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in, lat0
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele, ele2
type (wake_lr_struct), pointer :: lr(:)

integer i, ix, j, n_hom, n, nn, n_ele, ix_pass, i_lr, ie
integer irep
integer istep
integer :: nstep = 100

real(rp) deldr,dr, charge_try
real(rp) hom_voltage_gain, charge0, charge1, charge_old, charge_threshold
real(rp) trtb, currth, growth_rate, growth_rate0, growth_rate1, growth_rate_old
real(rp) time

logical lost

!New variables needed for regression test output
character(15) :: proto_th_str
character(15) :: proto_num_str
character(1) ::  nstep_str
character(13) ::  tail_str
character(29) :: final_th_str
character(29) :: final_num_str

namelist / bbu_params / bbu_param, beam_init, bmad_com

! Read in parameters

beam_init%n_particle = 1

print *,'================================================================'
print *,'================         BBU_PROGRAM        ===================='
print *,'================================================================'
print *,'Reading input file BBU.INIT'
print *

open (1, file = 'bbu.init', status = 'old')
read (1, nml = bbu_params)
close (1)

if (bbu_param%stable_orbit_anal) bbu_param%nstep = 1
if (bbu_param%regression) bbu_param%drscan = .true.
if (bbu_param%regression) bbu_param%nstep = 2

if (bbu_param%verbose) write(*,'(a,f6.1/,a,l/,a,i5/,a,l/)') &
        ' Maximum number of turns: ',bbu_param%simulation_turns_max, &
        ' DRSCAN analysis: ', bbu_param%drscan,&
        ' Number of repetitions: ',bbu_param%nrep,&
        ' Stable orbit analysis: ',bbu_param%stable_orbit_anal

! Define distance between bunches

beam_init%dt_bunch = 1 / bbu_param%bunch_freq
call ran_seed_put (bbu_param%ran_seed)
if (bbu_param%verbose) print *, 'Random number seed:', bbu_param%ran_seed

if (bbu_param%ran_gauss_sigma_cut > 0) then
  call ran_gauss_converter (set_sigma_cut = bbu_param%ran_gauss_sigma_cut)
  if (bbu_param%verbose) print *, 'ran_gauss sigma cut: ', bbu_param%ran_gauss_sigma_cut 
endif

! Init

if (bbu_param%verbose) print *, 'Lattice file: ', trim(bbu_param%lat_file_name)
call bmad_parser (bbu_param%lat_file_name, lat_in)
call twiss_propagate_all (lat_in)
call lat_make_mat6(lat_in)  ! Necessary if a match lattice element is present.
call run_timer ('START')

! Remove HOM's of higher order

if (bbu_param%hom_order_cutoff > 0) then
  do i = 1, lat_in%n_ele_max
    ele => lat_in%ele(i)
    if (.not. associated(ele%wake)) cycle
    if (.not. allocated(ele%wake%lr)) cycle
    n = count(ele%wake%lr(:)%m > bbu_param%hom_order_cutoff)
    if (n == 0) cycle  ! Nothing to remove
    if (n == size(ele%wake%lr)) then
      deallocate (ele%wake%lr)
      cycle
    endif
    lr => ele%wake%lr
    nn = size(ele%wake%lr) - n
    allocate(ele%wake%lr(nn))
    n = 0
    do j = 1, size(lr)
      if (lr(j)%m > bbu_param%hom_order_cutoff) cycle
      n = n + 1; ele%wake%lr(n) = lr(j)
    enddo
    deallocate(lr)
  enddo
endif

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
      if (bbu_param%verbose) write(6,'(a,i6,a,a40)')&
            ' Element index ',i,' found for DR scan element ',bbu_param%elname
! Open DRSCAN output file for threshold current calculation comparison
      if(.not. bbu_param%regression) open (50, file = 'drscan.out', status = 'unknown')
      if(bbu_param%regression) open (40, file = 'output.now', status = 'unknown') 
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
  if (bbu_param%verbose) write(6,'(a,i10,a)')&
        ' Opening output file REP.OUT for',bbu_param%nrep,' repetitions'
  open (55, file = 'rep.out', status = 'unknown')
endif

! Open file for stable orbit analysis data
if (bbu_param%stable_orbit_anal)then
  if (bbu_param%verbose) write(6,'(a,i10,a)')&
        ' Opening output file STABLE_ORBIT.OUT for stable orbit analysis'
  open (56, file = 'stable_orbit.out', status = 'unknown')

! Open file for HOM voltage stability data
  if (bbu_param%verbose) write(6,'(a,i10,a)')&
        ' Opening output file HOM_VOLTAGE.OUT for stable orbit analysis'
  open (57, file = 'hom_voltage.out', status = 'unknown')

endif

!-------------------------------------------------------
! DRSCAN loop. NSTEP=1 if no DRSCAN.

do istep = 1, nstep

  if (bbu_param%drscan) then
    ! Change length of taylor element
    deldr = 0.0
    if(nstep > 1) deldr = (bbu_param%enddr - bbu_param%begdr)/(nstep-1)
    dr = bbu_param%begdr + (istep-1) * deldr
    ie = bbu_param%elindex
    lat_in%ele(ie)%value(l$) = dr * c_light / bbu_param%bunch_freq
    if (bbu_param%verbose) write(6,'(a,2f8.3)')' DRSCAN analysis step: dr, scan element length = ', &
                 dr, lat_in%ele(bbu_param%elindex)%value(l$)
    call set_flags_for_changed_attribute (lat_in%ele(ie), lat_in%ele(ie)%value(l$))
    call lattice_bookkeeper(lat_in)
  endif

  if (bbu_param%hybridize) then
    if (bbu_param%verbose) print *, 'Note: Hybridizing lattice...'
    do i = 1, lat_in%n_ele_max
      ! Keep element if defined as end of tracking
      ele => lat_in%ele(i)
      ele%select = .false.
      if (ele%name == bbu_param%ele_track_end) then
        ele%select = .true.
        cycle
      endif
      if (ele%key /= lcavity$) cycle
      if (.not. bbu_param%keep_all_lcavities) then
        if (.not. associated (ele%wake)) cycle
        if (size(ele%wake%lr) == 0) cycle
      endif
      ele%select = .true.
    enddo
    call make_hybrid_lat (lat_in, lat, bbu_param%use_taylor_for_hybrids)
  else
    lat = lat_in
  endif

  lat0 = lat

  ! Define element at which tracking ends

  if (bbu_param%ele_track_end.ne.' ')then
      call lat_ele_locator(bbu_param%ele_track_end,lat, eles, n_loc, err)
      if(err)call err_exit
      if(n_loc.eq.0)then
       if (bbu_param%verbose) print '(2a)', 'No matching element found for ',bbu_param%ele_track_end  
       call err_exit
      elseif(n_loc.gt.1) then
       if (bbu_param%verbose) print '(2a)', 'Multiple matching elements found for ',bbu_param%ele_track_end  
       if (bbu_param%verbose) print '(a)', 'Will use the first instance as the end of the tracking'
      endif
      ele => eles(1)%ele
      if (eles(1)%ele%lord_status == super_lord$) ele => pointer_to_slave(ele, ele%n_slave)
      ix = ele%ix_ele
      if (ix > lat%n_ele_track) then
         if (bbu_param%verbose) print *, 'STOPPING ELEMENT IS A LORD! ', bbu_param%ele_track_end
         call err_exit
      endif
      bbu_param%ix_ele_track_end = ix
      if (bbu_param%verbose) print *,' Tracking will be halted after element ',bbu_param%ele_track_end
  endif

  !

  if (bbu_param%write_hom_info) then
   call write_homs(lat, bbu_param%bunch_freq, trtb, currth,bbu_param%verbose)
  ! Update starting current according to analytic approximation
   if (currth.gt.0.)bbu_param%current = currth
  else
   trtb=dr
  endif

  call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)

  ! Print some information and get the analytic approximation result for the threshold current

  if (bbu_param%verbose) print '(2a)', ' Lattice File: ', trim(lat%input_file_name)
  if (bbu_param%verbose) print '(2a)', ' Lattice Name: ', trim(lat%lattice)
  if (bbu_param%verbose) print '(a, i7)', ' Num Elements to Track: ', bbu_param%ix_ele_track_end
  if (bbu_param%verbose) print '(a, i7)', ' Num Elements Elements: ', lat%n_ele_track
  if (bbu_param%verbose) print '(a, es12.2)', ' Beam Energy: ', lat%ele(0)%value(e_tot$)

  if (bbu_param%verbose) print *, 'Number of lr wake elements in tracking lattice:', size(bbu_beam%stage)

  n_ele = 0
  do i = 1, size(bbu_beam%stage)
    j = bbu_beam%stage(i)%ix_ele_lr_wake
    call multipass_chain (lat%ele(j), ix_pass, n_links = n)
    if (ix_pass /= 1 .and. n /= 0) cycle
    n_ele = n_ele + 1
  enddo

  if (bbu_param%verbose) print *, 'Number of physical lr wake elements:', n_ele
  if (bbu_param%verbose) print *, 'Number of elements in lattice:      ', lat%n_ele_track

  ! Loop over calculation repetitions
  do irep = 1,bbu_param%nrep

    beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

    if (bbu_param%verbose) print *, '  Bunch charge (pC): ', 1e12_rp*beam_init%bunch_charge

    if (bbu_param%stable_orbit_anal) then

      if (bbu_param%verbose) print *,' Analyzing stable orbit for repetition ', irep

      lat = lat0 ! Restore lr wakes
      call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)

      if (lost) then
         if (bbu_param%verbose) print *, 'PARTICLE(S) LOST. ASSUMING UNSTABLE...'
      else
    ! Print output for stable orbit analysis
                  do i = 1, size(bbu_beam%stage)
                     j = bbu_beam%stage(i)%ix_stage_pass1
          write(56,'(i10,5(1x,e15.8),i10/,4(6(1x,e15.6)/))')i,bbu_beam%stage(i)%time_at_wake_ele, &
                      real(bbu_beam%stage(j)%hom_voltage_max), &
                      real(lat%ele(bbu_beam%stage(i)%ix_ele_lr_wake)%wake%lr(bbu_beam%stage(j)%ix_hom_max)%freq), &
                      real(lat%ele(bbu_beam%stage(i)%ix_ele_lr_wake)%wake%lr(bbu_beam%stage(j)%ix_hom_max)%R_over_Q), &
                      real(lat%ele(bbu_beam%stage(i)%ix_ele_lr_wake)%wake%lr(bbu_beam%stage(j)%ix_hom_max)%Q), &
                      bbu_beam%stage(i)%n_orb, & 
                      bbu_beam%stage(i)%ave_orb, bbu_beam%stage(i)%rms_orb, &
                      bbu_beam%stage(i)%min_orb, bbu_beam%stage(i)%max_orb
        enddo

      endif

    ! Re-randomize HOM frequencies

      do i = 1, lat0%n_ele_max
        call randomize_lr_wake_frequencies (lat0%ele(i))
      enddo
      ! Skip threshold calculation
      cycle

    endif

    ! If not bbu_param%stable_orbit_anal...

    charge0 = 0
    charge1 = -1      ! Mark as not set yet 
    charge_old = -1   ! Mark as not set yet 
    charge_try = -1   ! Mark as not set yet 

    if (bbu_param%verbose) Print *, 'Searching for a current where the tracking is unstable...'

    do
      lat = lat0 ! Restore lr wakes
      call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
      call calc_next_charge_try
      if (hom_voltage_gain > 1) exit
      if (lost) then
        if (bbu_param%verbose) print *, 'PARTICLE(S) LOST. ASSUMING UNSTABLE...'
        exit
      endif
    enddo

     ! Track to bracket threshold

      if (bbu_param%verbose) print *, 'Now converging on the threshold...'

      do
        lat = lat0 ! Restore lr wakes
        call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
        if (lost .and. bbu_param%verbose) print *, 'Particle(s) lost. Assuming unstable...'
        call calc_next_charge_try
        if (charge1 - charge0 < charge1 * bbu_param%rel_tol) exit
      enddo

      beam_init%bunch_charge = (charge0 + charge1) / 2
      if (bbu_param%verbose) print *, 'Threshold Current (A):', beam_init%bunch_charge / beam_init%dt_bunch 
      i = bbu_beam%ix_stage_voltage_max
      j = bbu_beam%stage(i)%ix_ele_lr_wake
      ele => lat%ele(j)
      ele2 => ele
      if (ele2%slave_status == multipass_slave$) then
        ele2 => pointer_to_multipass_lord (ele)
      endif
      i_lr = bbu_beam%stage(i)%ix_hom_max
      if (bbu_param%verbose) print *, 'Element with critical HOM:', ele2%ix_ele, ':   ', ele2%name
      if (bbu_param%verbose) print *, 'Critical HOM: Input Frequency: ', ele%wake%lr(i_lr)%freq_in 
      if (bbu_param%verbose) print *, 'Critical HOM: Actual Frequency:', ele%wake%lr(i_lr)%freq
      if (bbu_param%verbose) print *, 'Critical HOM: R_overQ:         ', ele%wake%lr(i_lr)%r_over_q
      if (bbu_param%verbose) print *, 'Critical HOM: Q:               ', ele%wake%lr(i_lr)%q
      if (bbu_param%verbose) print *, 'Critical HOM: Angle:           ', ele%wake%lr(i_lr)%angle

      if (bbu_param%nrep.gt.1)write(55,'(i6,e14.6,2i7,6(e14.6,1x))') &
          irep, beam_init%bunch_charge / beam_init%dt_bunch, ele%ix_ele, ele2%ix_ele, ele%s, &
          ele%wake%lr(i_lr)%freq_in, ele%wake%lr(i_lr)%freq, ele%wake%lr(i_lr)%r_over_q, &
          ele%wake%lr(i_lr)%q, ele%wake%lr(i_lr)%angle

    ! Re-randomize HOM frequencies

    do i = 1, lat0%n_ele_max
      call randomize_lr_wake_frequencies (lat0%ele(i))
    enddo

  enddo  ! End of repetition loop

  if (bbu_param%drscan .and. .not. bbu_param%regression) write(50,*) trtb, currth, beam_init%bunch_charge / beam_init%dt_bunch
  if (bbu_param%regression) then
     proto_th_str = '"Theor_current('
     proto_num_str = '"Numer_current('
     if (istep == 1) nstep_str = '1'
     if (istep == 2) nstep_str = '2'
     tail_str = ')" REL  1E-10'
     final_th_str = proto_th_str//nstep_str//tail_str
     final_num_str = proto_num_str//nstep_str//tail_str
     
     write(40, '(a, es22.15)') final_th_str, currth
     write(40, '(a, es22.15)') final_num_str, beam_init%bunch_charge / beam_init%dt_bunch
  end if
  
enddo  ! End of DRSCAN loop

if (bbu_param%drscan .and. .not. bbu_param%regression) close(50)
if (bbu_param%regression) close(40)
if (bbu_param%nrep.gt.1)close(55)
if (bbu_param%stable_orbit_anal) then
 close(56)
 close(57)
endif

call run_timer ('STOP', time)
print *
print *, 'Time for calculation (min): ', time/60
print *
print *,'================================================================'
print *,'================      Exit BBU_PROGRAM       ==================='
print *,'================================================================'

!-------------------------------------------------------------------------------
contains

subroutine calc_next_charge_try()

real(rp) c, min_delta, dc0, dc1, c0, c1, g0, g1

! Print info

if (growth_rate > 0 .or. lost) then
  if (bbu_param%verbose) print *, '  Unstable at (mA):', 1e3 * beam_init%bunch_charge / beam_init%dt_bunch 
else
  if (bbu_param%verbose) print *, '  Stable at (mA):', 1e3 * beam_init%bunch_charge / beam_init%dt_bunch 
endif

if (bbu_param%verbose) print *, '         Head bunch index: ', bbu_beam%bunch(bbu_beam%ix_bunch_head)%ix_bunch
if (bbu_param%verbose) print *, '         Growth rate: ', growth_rate

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


  if (bbu_param%verbose) print *, '         Predicted threshold:', 1d3 * charge_threshold / beam_init%dt_bunch 

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

if (bbu_param%verbose) print *, '         Current to try next:', 1d3 * beam_init%bunch_charge / beam_init%dt_bunch 

end subroutine

end program
