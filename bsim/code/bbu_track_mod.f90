module bbu_track_mod

use bmad
use beam_mod

type bbu_stage_struct
  integer :: ix_ele_lr_wake   ! Element index of element with the wake 
  integer :: ix_ele_stage_end ! Element at end of stage.
  integer :: ix_pass          ! Pass index when in multipass section
  integer :: ix_stage_pass1   ! Index of corresponding stage on first pass
  integer :: ix_head_bunch
  integer :: ix_hom_max
  real(rp) :: hom_voltage_max
  real(rp) time_at_wake_ele 
  real(rp) :: ave_orb(6) = 0, rms_orb(6) = 0, min_orb(6) = 0, max_orb(6) = 0
  integer n_orb
end type

type bbu_beam_struct
  type (bunch_struct), allocatable :: bunch(:)  ! Bunches in the lattice
  type (bbu_stage_struct), allocatable :: stage(:)
  integer, allocatable :: ix_ele_bunch(:)       ! element where bunch is 
  integer ix_bunch_head       ! Index to head bunch(:)
  integer ix_bunch_end        ! Index of the end bunch(:). -1 -> no bunches.
  integer n_bunch_in_lat      ! Number of bunches transversing the lattice.
  integer ix_stage_voltage_max
  real(rp) hom_voltage_max
  real(rp) time_now
  real(rp) one_turn_time
  real(rp) rf_wavelength_max
end type

type bbu_param_struct
  character(500) :: lat_filename = 'erl.lat'     ! Bmad lattice file name
  character(500) :: lat2_filename = ''     ! Bmad lattice2 file name for secondary parser
  character(100) :: bunch_by_bunch_info_file = '' ! For outputting bunch-by-bunch info.
  logical :: hybridize = .true.                  ! Combine non-hom elements to speed up simulation?
  logical :: write_digested_hybrid_lat = .false. ! For debugging purposes.
  logical :: write_voltage_vs_time_dat = .false. ! For debugging purposes.
  logical :: keep_overlays_and_groups = .false.  ! Keep when hybridizing?
  logical :: keep_all_lcavities  = .false.       ! Keep when hybridizing?
  logical :: use_taylor_for_hybrids = .false.    ! Use taylor map for hybrids when true. Otherwise tracking method is linear.
  logical :: stable_orbit_anal = .false.         ! Write stable_orbit.out and hom_voltage.out?
  real(rp) :: limit_factor = 2                   ! Init_hom_amp * limit_factor = simulation unstable limit
  real(rp) :: simulation_turns_max = 20          ! Sets the duration of the simulation.
  real(rp) :: bunch_freq = 1.3e9                 ! Freq in Hz.
  real(rp) :: init_particle_offset = 1e-8        ! Initial particle offset for particles born in the first turn period.
  real(rp) :: current = 20e-3                    ! Starting current (amps)
  real(rp) :: rel_tol = 1e-2                     ! Final threshold current accuracy. 
  logical :: drscan = .true.                     ! If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3.
  logical :: use_interpolated_threshold = .true.
  logical :: write_hom_info = .true.             ! Write HOM parameters to main output file?
  integer :: elindex
  character*40 :: elname = 'T1'                  ! Element to step length for DRSCAN
  integer :: nstep = 100                         ! Number of steps for DRSCAN.      
  real(rp) :: begdr = 5.234                      ! Beginning DR value for DRSCAN.          
  real(rp) :: enddr = 6.135                      ! End DR value for DRSCAN.
  integer :: nrep  = 1                           ! Number of times to repeat threshold calculation
  integer :: ran_seed = 0                        ! If set to 0, the output results will vary from run to run.
  integer :: hom_order_cutoff = 0                ! If positive -> ignore HOM's with order greater than this.
  real(rp) :: ran_gauss_sigma_cut = -1
  character(40) :: ele_track_end = ' '
  integer :: ix_ele_track_end = -1               ! Default: set to last element with a wake
  logical :: regression = .false.                ! Do regression test?
  logical :: normalize_z_to_rf = .false.         ! make starting z = mod(z, rf_wavelength)?

  ! Ramp parameters
  logical :: ramp_on = .false.
  real(rp) :: ramp_pattern(1000) = real_garbage$
  integer :: ramp_n_start = 0                    ! Index of start of ramp

  ! Internal parameters
  integer :: n_ramp_pattern = -1                       ! Number of valid ramp_pattern
end type

contains

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_setup (lat, dt_bunch, bbu_param, bbu_beam)

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct) bbu_beam
type (beam_init_struct) beam_init
type (bbu_param_struct) bbu_param
type (ele_struct), pointer :: wake_ele
type (ele_pointer_struct), allocatable, save :: chain_ele(:)

integer i, j, k, ih, ixs, ix_pass, n_links, n_stage, ix_track_end

real(rp) dt_bunch, rr(4)

character(16) :: r_name = 'bbu_setup'

! Size bbu_beam%bunch
if (dt_bunch == 0) then
  call out_io (s_fatal$, r_name, 'DT_BUNCH IS ZERO!')
  call err_exit
endif

! Find all elements that have a lr wake.
! This determines the number of stages.

n_stage = 0
do i = 1, lat%n_ele_track
  wake_ele => lat%ele(i)
  if (wake_ele%slave_status == super_slave$) then
    wake_ele => pointer_to_lord(wake_ele, 1, ix_slave_back = ixs)
    if (ixs /= wake_ele%n_slave) cycle  ! Only use the last slave of a super_lord
  endif
  if (.not. associated(wake_ele%wake)) cycle
  if (size(wake_ele%wake%lr%mode) == 0) cycle
  n_stage = n_stage + 1
enddo

if (n_stage == 0) then
  call out_io (s_fatal$, r_name, 'NO LR WAKES FOUND IN LATTICE!')
  call err_exit
endif

if (allocated(bbu_beam%stage)) deallocate (bbu_beam%stage)
allocate (bbu_beam%stage(n_stage))

! Bunches have to go through a given physical lcavity in the correct order. 
! To do this correctly, the lattice is divided up into stages.
! A given stage has one and only one lcavity.
! The first stage starts at the beginning of the lattice.
! The last stage ends at the end of the lattice.
! All stages, except the last one, end at an lcavity.

! bbu_beam%stage(i)%ix_ele_lr_wake holds the index in lat%ele(:) of the lcavity of the i^th stage.
! bbu_beam%stage(i)%ix_ele_stage_end holds the index in lat%ele(:) of the last element of the i^th stage.
! Normally these two are the same except when the lcavity is a super_lord. In this case %ix_ele_end_stage
! will be the last super_slave.

! bbu_beam%stage(i)%ix_head_bunch holds the index in bbu_beam%bunch(:) of the head 
! bunch for the i^th stage.
! The head bunch is the next bunch that will be propagated through the stage.

bbu_beam%stage%ix_head_bunch = -1    ! Indicates there are no bunches ready for a stage.
bbu_beam%stage(1)%ix_head_bunch = 1

j = 0   !The "stage" index
do i = 1, lat%n_ele_track
  wake_ele => lat%ele(i)
  if (wake_ele%slave_status == super_slave$) then
    wake_ele => pointer_to_lord(wake_ele, 1, ix_slave_back = ixs)
    if (ixs /= wake_ele%n_slave) cycle  ! Only consider the last slave of a super_lord
  endif
  if (.not. associated(wake_ele%wake)) cycle
  if (size(wake_ele%wake%lr%mode) == 0) cycle
  j = j + 1
  
  !stage j holds the index of the "wake element", i
  bbu_beam%stage(j)%ix_ele_stage_end = i
  bbu_beam%stage(j)%ix_ele_lr_wake = wake_ele%ix_ele

  ! Returns the chain of multipass elements that represent the SAME physical element 
  ! ix_pass: multipass number of the element, -1 if not a multipass element
  ! n_links: number of pass through the element
  call multipass_chain (wake_ele, ix_pass, n_links, chain_ele)
  bbu_beam%stage(j)%ix_pass = ix_pass

  ! This fills out bbu_beam%stage(j)%ix_stage_pass1, the indes of the
  ! corresponding STAGE on it's first pass 
  ! That is, a multipass stage with pass number >1 has ix_stage_pass1 set to
  ! the corresponding stage with pass number = 1
  if (ix_pass > 0) then
    do k = 1, j
      if (bbu_beam%stage(k)%ix_ele_lr_wake == chain_ele(1)%ele%ix_ele) then
        bbu_beam%stage(j)%ix_stage_pass1 = k
        exit
      endif
    enddo
    if (k > j) then
      print *, 'BOOKKEEPING ERROR'
      call out_io (s_fatal$, r_name, 'BOOKKEEPING ERROR.')
      call err_exit
    endif
  else
    bbu_beam%stage(j)%ix_stage_pass1 = j
  endif
enddo

! Bunch init
if (bbu_param%ix_ele_track_end > 0) bbu_beam%stage(n_stage)%ix_ele_stage_end = bbu_param%ix_ele_track_end 
ix_track_end = bbu_beam%stage(n_stage)%ix_ele_stage_end

bbu_beam%n_bunch_in_lat = (lat%ele(ix_track_end)%ref_time / dt_bunch) + 1
bbu_beam%one_turn_time = lat%ele(ix_track_end)%ref_time

if (allocated(bbu_beam%bunch)) deallocate (bbu_beam%bunch)
allocate(bbu_beam%bunch(bbu_beam%n_bunch_in_lat+10))
bbu_beam%ix_bunch_head = 1
bbu_beam%ix_bunch_end = -1  ! Indicates: No bunches

call re_allocate (bbu_beam%ix_ele_bunch, bbu_beam%n_bunch_in_lat+10)

! RF_wavelength_max calc



end subroutine bbu_setup

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_normalized, growth_rate, lost, irep) 

implicit none

type (lat_struct) lat
type (ele_struct), pointer :: end_ele
type (ele_struct), pointer :: ele
type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (beam_init_struct) beam_init
type (coord_struct), allocatable :: orbit(:)

real(rp) hom_voltage0, hom_voltage_sum, r_period, r_period0, voltage
real(rp) hom_voltage_normalized, growth_rate, orb(6)
real(rp) max_x,rms_x,max_y,rms_y

integer i, n, n_period, n_count, n_period_old, irep
integer, save :: m = 0
integer n_settling_period
integer final_count, ix_stage_tracked

real(rp) hom_voltage_normalized_final_sum, hom_voltage_normalized_final_avg


logical lost

character(16) :: r_name = 'bbu_track_all'

! Setup.
call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)
print *, 'bbu_setup complete.'

!! Populate bunches into the lattice
!print *,'Initial number of bunches:', size(bbu_beam%bunch)
do i = 1, size(bbu_beam%bunch)
  call bbu_add_a_bunch (lat, bbu_beam, bbu_param, beam_init)
enddo
print *, 'initial bunch population complete.'

call reallocate_coord (orbit, lat%n_ele_track)

! init time at hom and hom_voltage
bbu_beam%stage%hom_voltage_max = 0
bbu_beam%hom_voltage_max = 0
bbu_beam%ix_stage_voltage_max = 1
growth_rate = real_garbage$
r_period = 0

if (bbu_param%write_voltage_vs_time_dat) then
  !open(28, file = 'volt_v_turn2.txt', status = 'unknown', access = 'append')
  open(29, file = 'volt_v_turn.txt', status = 'unknown', access = 'append')
  write(29,'(a15, es12.4)') 'current:', bbu_param%current
  write(29,'(a15, es12.4)') 'bunch_freq:', bbu_param%bunch_freq
  write(29,'(a15, es12.4)') 'simulation_turn_max:', bbu_param%simulation_turns_max
  write(29,'(a15, i8)') 'ran_seed:', bbu_param%ran_seed
  write(29,'(a15, a15, a15, a15)') 'n_period_old', 'time_now', 'hom_vol_max'
endif

bbu_beam%stage%time_at_wake_ele = 1e30  ! something large
end_ele => lat%ele(bbu_beam%stage(1)%ix_ele_stage_end)
bbu_beam%stage(1)%time_at_wake_ele = bbu_beam%bunch(1)%t_center + end_ele%ref_time

! Track...
! To decide whether things are stable or unstable we look at the HOM voltage integrated
! over one turn period. The voltage is integrated over one period since there can 
! be a large variation in HOM voltage within a period.
! And since it is in the first period that all the stages are getting populated with bunches,
! We don't start integrating until the second period.

n_period_old = 1
hom_voltage_sum = 0
n_count = 0
final_count = 0
hom_voltage_normalized_final_sum = 0 
hom_voltage_normalized_final_avg = 0

print *, 'bbu_track_all loop begins, TEST CURRENT: ', bbu_param%current

do
  ! Track the beam until beam loss or gain too small or large 
  call bbu_track_a_stage (lat, bbu_beam, bbu_param, lost, ix_stage_tracked)
  if (lost) then
    print *, 'Beam is lost after tracking a stage. (HOM voltage computed might be garbage)'
    exit
  endif

  ! If the head bunch is finished then remove it and seed a new one.
  if (ix_stage_tracked == size(bbu_beam%stage)) then
    call bbu_remove_head_bunch (bbu_beam)
    call bbu_add_a_bunch (lat, bbu_beam, bbu_param, beam_init)
  endif

  !print *, bbu_beam%time_now,bbu_beam%one_turn_time
  ! Test for the end of the period
  r_period = bbu_beam%time_now / bbu_beam%one_turn_time  !! Ticking up "continuously" (dt based on bbu_track_a_stage) 
  n_period = int(r_period) + 1                           !! Ticking up at integer steps

  n_settling_period = 10

  !!! This huge IF statement is activated only at the beginning of a new period, when n_period ticks up
  ! n_period_old is updated to n_period at the end of the statement, so the 
  ! IF statement is not entered until n_period ticks up by 1 
  if (n_period /= n_period_old) then
    !print *, n_period
    
    ! computing orb for this period ?
    !! Commeted out Mar 31 2017
    !do i = 1, size(bbu_beam%stage)
    !  if (bbu_beam%stage(i)%n_orb == 0) cycle
    !  bbu_beam%stage(i)%ave_orb  = bbu_beam%stage(i)%ave_orb / bbu_beam%stage(i)%n_orb 
    !  orb  = bbu_beam%stage(i)%rms_orb/bbu_beam%stage(i)%n_orb - bbu_beam%stage(i)%ave_orb**2
    !  where (orb < 0) orb = 0
    !  bbu_beam%stage(i)%rms_orb = sqrt(orb)
    !enddo
    !print *, "sqrt:",size(orb),sqrt(orb)

    if (n_period_old == n_settling_period) then
      ! The baseline/reference hom voltage is computed over the 2nd period after the offset bunches 
      ! have passed through the lattice.
      hom_voltage0 = hom_voltage_sum / n_count
      r_period0 = r_period
    endif

    if (n_period_old > n_settling_period) then
      hom_voltage_normalized = (hom_voltage_sum / n_count) / hom_voltage0
      !!! Showing detailed HOM voltage for all periods after settling
      !! For debug
      !print *, '-------------------------------------------'
      !print *, 'period #: ', n_period_old
      !print *, 'HOM voltage0: ', hom_voltage0
      !print *, 'HOM voltage sum: ', hom_voltage_sum
      !print *, 'Normalized HOM voltage', hom_voltage_normalized
      
      if (bbu_param%simulation_turns_max - n_period_old < 5) then
        final_count = final_count + 1
        hom_voltage_normalized_final_sum = hom_voltage_normalized_final_sum + hom_voltage_normalized  
        hom_voltage_normalized_final_avg = hom_voltage_normalized_final_sum / final_count 
        print *, ' Avg over the final ', final_count,' periods: ', hom_voltage_normalized_final_avg 
      endif

      !!growth_rate = log(hom_voltage_normalized) / (r_period - r_period0)
      growth_rate = log(hom_voltage_normalized) 
    else
      hom_voltage_normalized = -1   ! Dummy value for output
    endif

!    write(*,'(a,i8,a,es15.6,5x,a,es15.6/,a,es13.6,a,i9,a/,a,es13.6)')&
!         ' Period ', n_period_old,'   Time: ',bbu_beam%time_now, &
!         ' Beam current(A): ',beam_init%bunch_charge / beam_init%dt_bunch, &
!         ' Sum of max HOM wake amplitudes (V): ', hom_voltage_sum, &
!         '   over ', n_count, ' bunch passages', &
!         ' HOM normalized voltage (gain factor): ', hom_voltage_normalized

    if (n_period_old > n_settling_period) then
      if (.not. bbu_param%stable_orbit_anal) then  ! default stable_orbit_anal is false
        if (hom_voltage_normalized < 1/bbu_param%limit_factor) then
          print *, 'Tracking stops b/c hom_voltage_normalized BELOW 1/limit_factor'
          exit  ! default limit_factor is 2
        endif
        if (hom_voltage_normalized > bbu_param%limit_factor) then
          print *, 'Tracking stops b/c hom_voltage_normalized ABOVE limit_factor'
          exit
        endif  
      else
!  write(57,'(2i10,e13.6,x,i8,x,e15.6)')irep,n_period_old,hom_voltage_sum,n_count,hom_voltage_normalized
      endif
    endif

    !if (r_period > bbu_param%simulation_turns_max - n_settling_period) then  
      !! loop ends when r_period > simulation_turns_max
      !print *, 'r_period: ', r_period
      !print *, 'Stage #:', bbu_beam%ix_stage_voltage_max, 'HOM #: ', bbu_beam%stage(bbu_beam%ix_stage_voltage_max)%ix_hom_max  
    !endif
    
    ! Exit loop over periods of max period reached
    if (r_period > bbu_param%simulation_turns_max) then
      !! These 2 lines are optional. It allows the final hom_voltage_normalized to be
      !! taken as the average over the last few periods, instead of just the final period
      !! for which the value can be slightly >1.0 (numerical noise), 
      !! but the test current is actually physically stable. 
      !! In short, this prevents the program to be fooled by numerical noise. 
      hom_voltage_normalized = hom_voltage_normalized_final_avg 
      growth_rate = log(hom_voltage_normalized) 
      !! This line exits the loop, and the hom_voltage_normalized (or its log) will be used to
      !! determine the stability of the test current in an outside (python) program.
      call out_io (s_warn$, r_name, 'Simulation turns max exceeded -- Ending Tracking  \f10.6\ ', &
                                                   r_array = (/ hom_voltage_normalized /) )
      exit
    endif

    ! Reset parameters for the next n_period 
    hom_voltage_sum = 0
    n_count = 0
    n_period_old = n_period
    do i = 1, size(bbu_beam%stage)
      bbu_beam%stage(i)%ave_orb = 0
      bbu_beam%stage(i)%rms_orb = 0
      bbu_beam%stage(i)%min_orb = 1     ! Something large
      bbu_beam%stage(i)%max_orb = -1    ! Something small
      bbu_beam%stage(i)%n_orb = 0  
    enddo
  endif

  ! Compute integrated voltage.
  ! Specifically, update bbu_beam%voltage_max after every track_a_stage
  ! and sum them over for each n_period
  ! The sum is reset when n_period ticks up by 1
  call bbu_hom_voltage_calc (lat, bbu_beam, n_period, ix_stage_tracked)
  hom_voltage_sum = hom_voltage_sum + bbu_beam%hom_voltage_max
  n_count = n_count + 1 ! Number of track_a_stage in one period 
  !print *, 'period:', n_period   !For debug
  !print *, 'count:', n_count   !For debug
  !print *, 'hom_volt_max:',  bbu_beam%hom_voltage_max   !For debug

  if (bbu_param%write_voltage_vs_time_dat) then
    write(29,'(i4, es12.4, es12.4)') n_period_old, bbu_beam%time_now, bbu_beam%hom_voltage_max
  endif

!  if (bbu_param%write_voltage_vs_time_dat) then
!    do i=1, lat%n_ele_track
!      ele => lat%ele(i)
!      if (ele%key == lcavity$) then
!       !write(28,'(a12, i10, es15.6)') ele%name, n_period_old, hom_voltage_sum
!       !write(29,'(i10, es15.6)') n_period_old, hom_voltage_sum
!       !write(29,'(i4, i4, es15.6)') n_count, n_period_old, bbu_beam%hom_voltage_max
!        write(29,'(i4, es12.4, es12.4)') n_period_old, bbu_beam%time_now, bbu_beam%hom_voltage_max
!        exit
!      endif
!    enddo
!  endif

enddo

end subroutine bbu_track_all

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! This routine tracks one bunch through one stage.

subroutine bbu_track_a_stage (lat, bbu_beam, bbu_param, lost, ix_stage_tracked)

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct), target :: bbu_beam
type (bbu_param_struct) bbu_param
type (wake_lr_mode_struct), pointer :: lr
type (bbu_stage_struct), pointer :: this_stage
type (ele_struct), pointer :: end_ele

real(rp) min_time_at_wake_ele, time_at_wake_ele

integer i, ix_stage_tracked, ix_bunch, ix_ele, ie, ip
integer ix_ele_start, j, ib, ib2

character(*), parameter :: r_name = 'bbu_track_a_stage'

logical err, lost

! Look at each stage, track the bunch with the earliest time to finish, and track this stage.

ix_stage_tracked = minloc(bbu_beam%stage%time_at_wake_ele, 1)
this_stage => bbu_beam%stage(ix_stage_tracked)
bbu_beam%time_now = this_stage%time_at_wake_ele
end_ele => lat%ele(this_stage%ix_ele_stage_end)

ib = this_stage%ix_head_bunch
ix_ele_start = bbu_beam%bunch(ib)%ix_ele

!! Print the bunch at the first stage to observe bunch_pattern ( for debug )
!if (ix_stage_tracked == 1) then
!  print *, '------------------------------------------------- \n'
!  print *,'time_now(*1.3GHz): ',bbu_beam%time_now*(1.3*10**9)
!endif

! Track the bunch
open(999, file = 'bunch_vec.txt', status = 'unknown')
do j = ix_ele_start+1, end_ele%ix_ele

  call track1_bunch(bbu_beam%bunch(ib), lat%ele(j), err)

  ! For debug, to extract the bunch orbit.
  !print *, 'Orbit information saved to bunch_vec.txt '
  !write(999,'(a8, i8, a10, i8, a11, 6es11.3)') 'ele_ix:',j, 'bunch_ix:', ib, ' orbit:', bbu_beam%bunch(ib)%particle(1)%vec


  if (.not. all(bbu_beam%bunch(ib)%particle%state == alive$)) then
    print *, 'PARTICLE(S) LOST WHILE TRACKING ELEMENT: ', trim(lat%ele(j)%name), '  (', j, ')'
    lost = .true.
    return
  endif
enddo

! Write info to file if needed
if (bbu_param%bunch_by_bunch_info_file /= '') then
  call write_bunch_by_bunch_info (lat, bbu_beam, bbu_param, this_stage)
endif

! If the next stage does not have any bunches waiting to go through then the
! tracked bunch becomes the head bunch for that stage.

if (ix_stage_tracked /= size(bbu_beam%stage)) then  ! If not last stage
  if (bbu_beam%stage(ix_stage_tracked+1)%ix_head_bunch == -1) &
                      bbu_beam%stage(ix_stage_tracked+1)%ix_head_bunch = ib
endif

! If the bunch upstream from the tracked bunch is at the same stage as was tracked through,
! then this bunch becomes the new head bunch for the stage. Otherwise there is no head bunch
! for the stage.

if (ib == bbu_beam%ix_bunch_end) then
  call out_io (s_fatal$, r_name, 'NO BUNCHES FOR THE FIRST STAGE. GET HELP!')
  call err_exit
else
  ib2 = modulo (ib, size(bbu_beam%bunch)) + 1 ! Next bunch upstream
  if (bbu_beam%bunch(ib2)%ix_ele == ix_ele_start) then
    this_stage%ix_head_bunch = ib2
  else
    this_stage%ix_head_bunch = -1  ! No one waiting to go through this stage
  endif
endif

! Now correct min_time array
ix_bunch = this_stage%ix_head_bunch
if (ix_bunch < 0) then
  this_stage%time_at_wake_ele = 1e30  ! something large
else
  this_stage%time_at_wake_ele = bbu_beam%bunch(ix_bunch)%t_center + end_ele%ref_time
endif

if (ix_stage_tracked /= size(bbu_beam%stage)) then  ! If not last stage
  i = ix_stage_tracked + 1
  ix_bunch = bbu_beam%stage(i)%ix_head_bunch
  end_ele => lat%ele(bbu_beam%stage(i)%ix_ele_stage_end)
  bbu_beam%stage(i)%time_at_wake_ele = bbu_beam%bunch(ix_bunch)%t_center + end_ele%ref_time
endif

lost = .false.

end subroutine bbu_track_a_stage

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine bbu_add_a_bunch (lat, bbu_beam, bbu_param, beam_init)

implicit none

type (lat_struct) lat
type (bbu_beam_struct), target :: bbu_beam
type (beam_init_struct) beam_init
type (bunch_struct), pointer :: bunch
type (bbu_param_struct) bbu_param

integer i, ixb, ix0, ix_bunch
real(rp) r(2), t_rel, d_charge
real(rp) charge_diff, slope

character(20) :: r_name = 'bbu_add_a_bunch'

! Init bunch
if (bbu_beam%ix_bunch_end == -1) then ! if no bunches
  ixb = bbu_beam%ix_bunch_head
else
  ixb = modulo (bbu_beam%ix_bunch_end, size(bbu_beam%bunch)) + 1 ! Next empty slot
  if (ixb == bbu_beam%ix_bunch_head) then
    call out_io (s_fatal$, r_name, 'BBU_BEAM%BUNCH ARRAY OVERFLOW')
    call err_exit
  endif
endif

bunch => bbu_beam%bunch(ixb)

! If this is not the first bunch (most of time),
! we need to correct some of the bunch information

ix_bunch = 0
if (ixb /= bbu_beam%ix_bunch_head) then
  ix0 = bbu_beam%ix_bunch_end
  ix_bunch = bbu_beam%bunch(ix0)%ix_bunch + 1
endif

call init_bunch_distribution (lat%ele(0), lat%param, beam_init, ix_bunch, bunch)

!! Bunch pattern is achieved by multiplying bunch%charge_tot with d_charge(time) 
!! The "correct" current is not computed here 

if (bbu_param%ramp_on) then
  ixb = nint(bunch%t_center/beam_init%dt_bunch)
  if (ixb >= bbu_param%ramp_n_start) then
    i = mod(ixb - bbu_param%ramp_n_start, max(1, bbu_param%n_ramp_pattern - 1)) + 1
    d_charge = bbu_param%ramp_pattern(i)
    bunch%charge_tot = bunch%charge_tot * d_charge
    bunch%particle%charge = bunch%particle%charge * (d_charge / size(bunch%particle))
  endif
endif

bbu_beam%ix_bunch_end = ixb

! To excite HOMs initially, offset particles if the particle is born within 
!"the first turn period" (initial population of bunches in the lattice).

if (bunch%t_center < bbu_beam%one_turn_time) then
  do i = 1, size(bunch%particle)
    call ran_gauss (r)
    bunch%particle%vec(1) = bbu_param%init_particle_offset * r(1)
    bunch%particle%vec(3) = bbu_param%init_particle_offset * r(2)
  enddo
endif

end subroutine bbu_add_a_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine bbu_remove_head_bunch (bbu_beam)

implicit none

type (bbu_beam_struct) bbu_beam
character(20) :: r_name = 'bbu_remove_head_bunch'

! If there are no bunches then there is an error
if (bbu_beam%ix_bunch_end == -1) then
  call out_io (s_fatal$, r_name, 'TRYING TO REMOVE NON-EXISTANT BUNCH!')
  call err_exit
endif

! Mark ix_bunch_end if we are poping the last bunch.
if (bbu_beam%ix_bunch_end == bbu_beam%ix_bunch_head) bbu_beam%ix_bunch_end = -1

! Update ix_bunch_head 
bbu_beam%ix_bunch_head = modulo(bbu_beam%ix_bunch_head, size(bbu_beam%bunch)) + 1

end subroutine bbu_remove_head_bunch

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Calculates voltage in mode with maximal amplitude
! Specifically, update bbu_beam%ix_stage_voltage_max and bbu_beam%hom_voltage_max 

subroutine bbu_hom_voltage_calc (lat, bbu_beam, n_period, ix_stage_last_tracked)

implicit none

type (lat_struct), target :: lat
type (bbu_beam_struct) bbu_beam
type (wake_lr_mode_struct), pointer :: lr
type (ele_struct), pointer :: wake_ele

real(rp) hom_voltage_max, hom_voltage2
real(rp) hom_voltage_test

type (wake_lr_mode_struct), pointer :: lr_test

integer i, j, i1, ixm, ix_stage_last_tracked
integer o, ix_test, n_period
character(len=20) :: filename

! Only need to update the last stage tracked
hom_voltage_max = -1
i = ix_stage_last_tracked
i1 = bbu_beam%stage(i)%ix_stage_pass1    ! Find the corresponding stage of 1st pass 
wake_ele => lat%ele(bbu_beam%stage(i1)%ix_ele_lr_wake)   ! Find the wake element (cavity) of that corresponding stage
!! Find which wake of the stage has the max voltage
do j = 1, size(wake_ele%wake%lr%mode)  ! Number of lr wakes assigned to the cavity
  lr => wake_ele%wake%lr%mode(j)
  hom_voltage2 = max(lr%b_sin**2 + lr%b_cos**2, lr%a_sin**2 + lr%a_cos**2) 
  if (hom_voltage_max < hom_voltage2) then
    hom_voltage_max = hom_voltage2                 ! store the max voltage
    bbu_beam%stage(i1)%ix_hom_max = j              ! store the target wake index 
  endif
enddo

! Update the new max hom voltage in that cavity only
hom_voltage_max = sqrt(hom_voltage_max)
bbu_beam%stage(i1)%hom_voltage_max = hom_voltage_max        

! Compare the voltage with the stored universal max-voltage:  bbu_beam%hom_voltage_max
! Find the maximum hom voltage in any element.
!if (hom_voltage_max > bbu_beam%hom_voltage_max) then
!  bbu_beam%ix_stage_voltage_max = i1
!  bbu_beam%hom_voltage_max = hom_voltage_max

!! If not, get the maximum 
!! This only matters when the stage tracked IS the stage holding the old maximum
!elseif (i1 == bbu_beam%ix_stage_voltage_max) then
!  ixm = maxloc(bbu_beam%stage%hom_voltage_max, 1)
!  bbu_beam%ix_stage_voltage_max = ixm
!  bbu_beam%hom_voltage_max = bbu_beam%stage(ixm)%hom_voltage_max
!endif

ixm = maxloc(bbu_beam%stage%hom_voltage_max, 1)
bbu_beam%ix_stage_voltage_max = ixm
bbu_beam%hom_voltage_max = bbu_beam%stage(ixm)%hom_voltage_max
  

!! This section is for Voltage_evolution analysis
!! This will slow down the program a lot.
!! Make sure the only the 1st cavity is assigned with HOMs
ix_test = bbu_beam%stage(1)%ix_ele_lr_wake   ! Find the wake element (cavity) of that corresponding stage
lr_test => lat%ele(ix_test)%wake%lr%mode(1)
hom_voltage_test = sqrt( max(lr_test%b_sin**2 + lr_test%b_cos**2, lr_test%a_sin**2 + lr_test%a_cos**2)) 
!o = lunget()
!write (filename,"('Volt',I3.3,'.dat')") n_period  
!! Save voltage data every period
!open(o, file = filename, Access='append', status = 'unknown')
!write(o,'(es14.4E3,  es18.8E3)')  bbu_beam%time_now, hom_voltage_test
!close(o)

!print *,  bbu_beam%ix_stage_voltage_max, bbu_beam%stage(bbu_beam%ix_stage_voltage_max)%ix_hom_max  
end subroutine bbu_hom_voltage_calc


!!! If user didn't pre-assign the HOMs --------
!!! output 'hom_info.txt' so that HOMs can be assigned starting next bbu run --------
subroutine rf_cav_names (lat)
implicit none
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
integer i

!

open(20, file = 'hom_info.txt', status = 'unknown')
write(20,'(a15)') 'cavity_name'

do i = 1, lat%n_ele_max
  ele => lat%ele(i)
  if (ele%key /= lcavity$) cycle
  if (ele%slave_status == super_slave$) cycle
  if (ele%slave_status == multipass_slave$) cycle
  write(20,'(a)') trim(ele%name)
enddo

end subroutine rf_cav_names

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

!!!! output 'hom_info.txt' so that HOMs can be assigned starting next bbu run --------
subroutine check_rf_freq (lat, fb)
implicit none
type (bbu_param_struct) bbu_param
type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
integer i
character(25) cav_name
real(rp) fRF, fb

do i = 1, lat%n_ele_track
  ele => lat%ele(i)
  if (ele%key /= lcavity$) cycle
  if (ele%slave_status == multipass_slave$) ele => pointer_to_lord(ele, 1)
  cav_name = ele%name 
  fRF = ele%value(rf_frequency$) 
  
  if( abs(floor(fRF/fb) - fRF/fb) > 10E-6 ) then
    !print *, abs(floor(fRf/fb) - fb/fRF)
    print *,  'WARNING: For cavity ', cav_name, ' the RF frequency is not an integer &
    multiple of the input bunch frequency!! Bunches might be lost due to this!'  
  endif
enddo

end subroutine check_rf_freq

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

subroutine write_bunch_by_bunch_info (lat, bbu_beam, bbu_param, this_stage)

implicit none
type (lat_struct), target :: lat  !(input)
type (bbu_beam_struct), target :: bbu_beam !(input)
type (bunch_struct), pointer :: bunch
type (bbu_param_struct) bbu_param          !(input)
type (bbu_stage_struct), pointer :: this_stage    !(input)
type (wake_lr_struct), pointer :: lr

integer :: i, ios
integer, save :: iu = 0

!

if (iu == 0) then
  iu = lunget()
  open(iu, file = bbu_param%bunch_by_bunch_info_file, iostat = ios)
endif

bunch => bbu_beam%bunch(this_stage%ix_head_bunch)

! All HOMs:
do i=1, size(lat%ele(this_stage%ix_ele_lr_wake)%wake%lr%mode(:) )
  lr => lat%ele(this_stage%ix_ele_lr_wake)%wake%lr
  write(iu, '(4es15.7)') lr%t_ref, hom_voltage(lr%mode(i)), bunch%charge_live, bunch%particle(1)%vec(1)
enddo
!

end subroutine write_bunch_by_bunch_info

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

function hom_voltage(lr_wake) result(voltage)
implicit none
type(wake_lr_mode_struct) lr_wake
real(rp) :: voltage
!
voltage = max(hypot(lr_wake%a_sin,lr_wake%a_cos), hypot(lr_wake%b_sin,lr_wake%b_cos) )

end function

!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

function logical_to_python (logic) result (string)
implicit none

logical logic  
character(5) string

if (logic) then
  string = 'True'
else
  string = 'False'
endif

end function

end module
