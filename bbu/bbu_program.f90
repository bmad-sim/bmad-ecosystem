program bbu_program

use bbu_track_mod

implicit none

type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in, lat0
type (beam_init_struct) beam_init

integer i, ix, j, n_hom, n, n_ele, ix_pass

integer istep
integer nstep /100/
real(rp) deldr,dr

real(rp) hom_power0, hom_power1, charge0, charge1

real(rp) currth

logical lost
logical, allocatable :: keep_ele(:)

namelist / bbu_params / bbu_param, beam_init 

! Read in parameters

bbu_param%lat_file_name = 'erl.lat'  ! Bmad lattice file name
bbu_param%simulation_time = 1e-4     ! Time in sec.
bbu_param%bunch_freq = 1.3e9         ! Freq in Hz.
bbu_param%init_hom_amp = 1e-6        ! Initial wake amplitude        
bbu_param%limit_factor = 1e1         ! Init_hom_amp * limit_factor = simulation unstable limit
bbu_param%hybridize = .true.         ! Combine non-hom elements to speed up simulation?
bbu_param%keep_overlays_and_groups = .false. ! keep when hybridizing?
bbu_param%current = 20e-3            ! Starting current (amps)
bbu_param%rel_tol = 1e-3             ! Final threshold current accuracy.
bbu_param%write_hom_info = .true.  
bbu_param%num_stages_tracked_per_power_calc = 100 
bbu_param%drscan = .true.        ! If true, scan DR variable as in PRSTAB 7 (2004) Fig. 3.
bbu_param%nstep = 100
bbu_param%begdr = 5.234
bbu_param%enddr = 6.135

  beam_init%n_particle = 1

open (1, file = 'bbu.init', status = 'old')
read (1, nml = bbu_params)
close (1)

! Define distance between bunches
beam_init%ds_bunch = c_light / bbu_param%bunch_freq

! Init

print *, 'Lattice file: ', trim(bbu_param%lat_file_name)
call bmad_parser (bbu_param%lat_file_name, lat_in)
call twiss_propagate_all (lat_in)

! Initialize DR scan

nstep = 1
if (bbu_param%drscan) then
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

do istep = 1, nstep

  if (bbu_param%drscan) then
    ! Change length of taylor element
    deldr = 0.0
    if(nstep > 1) deldr = (bbu_param%enddr - bbu_param%begdr)/(nstep-1)
    dr = bbu_param%begdr + (istep-1) * deldr
    lat_in%ele(bbu_param%elindex)%value(l$) = dr * c_light / bbu_param%bunch_freq
    write(6,'(a,2f8.3)')' DRSCAN analysis step: tr/tb, scan element length = ', &
                 dr, lat_in%ele(bbu_param%elindex)%value(l$)
    call lattice_bookkeeper(lat_in)
  endif

  if (bbu_param%hybridize) then
    print *, 'Note: Hybridizing lattice...'
    allocate (keep_ele(lat_in%n_ele_max))
    keep_ele = .false.
    do i = 1, lat_in%n_ele_max
      if (lat_in%ele(i)%key /= lcavity$) cycle
      if (.not. associated (lat_in%ele(i)%wake)) cycle
      if (size(lat_in%ele(i)%wake%lr) == 0) cycle
      call update_hybrid_list (lat_in, i, keep_ele, bbu_param%keep_overlays_and_groups)
    enddo
    call make_hybrid_lat (lat_in, keep_ele, .true., lat)
    deallocate (keep_ele)
  else
    lat = lat_in
  endif

  lat0 = lat

  call bbu_setup (lat, beam_init%ds_bunch, bbu_param%init_hom_amp, bbu_beam)
  call bbu_hom_power_calc (lat, bbu_beam, hom_power0)

  bbu_param%high_power_lim = hom_power0 * bbu_param%limit_factor
  bbu_param%low_power_lim  = hom_power0 / bbu_param%limit_factor

  ! Print some information

  if (bbu_param%write_hom_info) call write_homs(lat,bbu_param%bunch_freq,currth)

  ! Update starting current according to analytic approximation
  if (currth.gt.0.)bbu_param%current = currth

  print *, 'Number of lr wake elements in tracking lattice:', size(bbu_beam%stage)

  n_ele = 0
  do i = 1, size(bbu_beam%stage)
    call multipass_chain (bbu_beam%stage(i)%ix_ele_lr_wake, lat, ix_pass, n_links = n)
    if (ix_pass /= 1 .and. n /= 0) cycle
    n_ele = n_ele + 1
  enddo

  print *, 'Number of physical lr wake elements:', n_ele

  ! Track to find upper limit

  beam_init%bunch_charge = bbu_param%current * beam_init%ds_bunch / c_light
  charge0 = 0

  Print *, 'Searching for a current where the tracking is unstable...'

  do
    lat = lat0 ! Restore lr wakes
    call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_power1, lost)
    if (hom_power1 > hom_power0) exit
    if (lost) then
      print *, 'Particle(s) lost stopping here.'
      stop
    endif
    charge0 = beam_init%bunch_charge
    print *, '  Stable at (mA):', 1e3 * charge0 * c_light / beam_init%ds_bunch 
    print *, '         Head bunch index: ', bbu_beam%bunch(bbu_beam%ix_bunch_head)%ix_bunch
    beam_init%bunch_charge = beam_init%bunch_charge * 2
  enddo

  charge1 = beam_init%bunch_charge
  print *, '  Unstable at (mA):', 1e3 * charge1 * c_light / beam_init%ds_bunch 
  print *, '         Head bunch index: ', bbu_beam%bunch(bbu_beam%ix_bunch_head)%ix_bunch

  ! Track to bracket threshold

  print *, 'Now converging on the threshold...'

  do
    beam_init%bunch_charge = (charge0 + charge1) / 2
    lat = lat0 ! Restore lr wakes
    call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_power1, lost)
    if (lost) then
      print *, 'Particle(s) lost stopping here.'
      stop
    endif
    if (hom_power1 > hom_power0) then
      charge1 = beam_init%bunch_charge
      print *, '  Unstable at (mA):', 1e3 * charge1 * c_light / beam_init%ds_bunch 
      print *, '         Head bunch index: ', bbu_beam%bunch(bbu_beam%ix_bunch_head)%ix_bunch
    else
      charge0 = beam_init%bunch_charge
      print *, '  Stable at (mA):', 1e3 * charge0 * c_light / beam_init%ds_bunch 
      print *, '         Head bunch index: ', bbu_beam%bunch(bbu_beam%ix_bunch_head)%ix_bunch
    endif

    if (charge1 - charge0 < charge1 * bbu_param%rel_tol) exit
  enddo

  beam_init%bunch_charge = (charge0 + charge1) / 2
  print *, 'Threshold Current (A):', beam_init%bunch_charge * c_light / beam_init%ds_bunch 

  if (bbu_param%drscan) write(50,*) dr, currth,beam_init%bunch_charge * c_light / beam_init%ds_bunch 

enddo

if(bbu_param%drscan)close(50)

end program
