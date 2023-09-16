program bbu_program

use bbu_track_mod

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in, lat0
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele
type (wake_lr_mode_struct), pointer :: lr(:)
type (coord_struct), allocatable :: orb(:) 

integer i, ix, j, n, nn, n_ele, ix_pass, o
integer irep
integer n_loc
real(rp) dr
real(rp) time, trtb, currth
real(rp) hom_voltage_gain
real(rp) growth_rate
logical err
logical lost

integer :: file_unit, file_unit2
integer k, status

character(200) init_file

namelist / bbu_params / bbu_param, beam_init, bmad_com


! Defaults for namelist
beam_init%n_particle = 1

select case (command_argument_count())
case (0)
  init_file = 'bbu.init'
case (1)
  call get_command_argument(1, init_file)
case default
  print *, 'CONFUSED: MULTIPLE COMMAND LINE ARGUMENTS!'
  stop
end select

print *,'Reading input file: ' // trim(init_file)
print *
open (1, file = init_file, status = 'old')
read (1, nml = bbu_params)
close (1)

! Define distance between bunches
beam_init%dt_bunch = 1 / bbu_param%bunch_freq

! Seed random number generator
call ran_seed_put (bbu_param%ran_seed)
if (bbu_param%ran_gauss_sigma_cut > 0) then
  call ran_gauss_converter (set_sigma_cut = bbu_param%ran_gauss_sigma_cut)
endif

! Init and parse
print *, 'Lattice file: ', trim(bbu_param%lat_filename)
call bmad_parser (bbu_param%lat_filename, lat_in) 

bmad_com%auto_bookkeeper = .false. ! To speed things up.

!For DR-scan, parse additional lattice (lat2) 
if (bbu_param%lat2_filename /= '') then
  print *, 'DR-scan or Phase-scan, parsing: ',bbu_param%lat2_filename
  call bmad_parser2 (bbu_param%lat2_filename, lat_in)
endif

!! Closed orbit computed
call twiss_and_track(lat_in,orb,status,0,.true.)

! Remove HOMs of higher order
if (bbu_param%hom_order_cutoff > 0) then
  do i = 1, lat_in%n_ele_max
    ele => lat_in%ele(i)
    ! Find cavity element with lr_wake
    if (.not. associated(ele%wake)) cycle
    if (.not. allocated(ele%wake%lr%mode)) cycle
    n = count(ele%wake%lr%mode(:)%m > bbu_param%hom_order_cutoff)
    if (n == 0) cycle  !All HOMs order <= cutoff m, nothing to remove
    if (n == size(ele%wake%lr%mode)) then  ! All HOMs order > m, remove this lcavity  
      deallocate (ele%wake%lr%mode)
      cycle
    endif
    !! If some (not all) HOMs order > m, extract the HOMs with order <= m 
    lr => ele%wake%lr%mode
    nn = size(ele%wake%lr%mode) - n    ! nn = number of HOMs to be kept
    allocate(ele%wake%lr%mode(nn))
    n = 0
    do j = 1, size(lr)
      if (lr(j)%m > bbu_param%hom_order_cutoff) cycle
      n = n + 1; ele%wake%lr%mode(n) = lr(j)
    enddo
    deallocate(lr)
  enddo
endif

!!!!!!!!!!!!!!!!!! HYBRIDIZATION !!!!!!!!!!!!!!!!!!!!!!!!
!! ele%select == true means the element will be kept, NOT hybridized 

if (bbu_param%hybridize) then
  do i = 1, lat_in%n_ele_max
    ele => lat_in%ele(i)
    ele%select = .false.     ! F => Hybridize
    ! Keep the element at the end of tracking, if specified by user
    if (ele%name == bbu_param%ele_track_end) then
      ele%select = .true.
      cycle
    endif

    if(ele%key == taylor$) then
      ele%select = .true.
      cycle
    endif    
    
    if (ele%key /= lcavity$) cycle
    if (.not. bbu_param%keep_all_lcavities) then
      if (.not. associated (ele%wake)) cycle
      if (size(ele%wake%lr%mode) == 0) cycle
    endif
    ele%select = .true.
  enddo

  call make_hybrid_lat (lat_in, lat, bbu_param%use_taylor_for_hybrids)
  print *, 'Hybridization complete !!!'

  if (bbu_param%write_digested_hybrid_lat) then
    call write_digested_bmad_file('hybrid.digested', lat)
    print *, 'Wrote hybrid lattice: hybrid.digested'
    call write_bmad_lattice_file('hybrid.lat', lat)
  endif

else ! If not hybridizing, keep the original lattice
  lat = lat_in
endif

! Keep the lattice ready to use?
lat0 = lat 

! Define element at which tracking ends, if user didn't specify one
if (bbu_param%ele_track_end.ne.' ') then
  call lat_ele_locator(bbu_param%ele_track_end,lat, eles, n_loc, err)
  if(err) call err_exit
  if(n_loc.eq.0)then
    print '(2a)', 'No matching element found for ',bbu_param%ele_track_end  
    call err_exit
  elseif(n_loc.gt.1) then
    print '(2a)', 'Multiple matching elements found for ',bbu_param%ele_track_end  
    print '(a)', 'Will use the first instance as the end of the tracking'
  endif
  ele => eles(1)%ele
  if (eles(1)%ele%lord_status == super_lord$) ele => pointer_to_slave(ele, ele%n_slave)
  ix = ele%ix_ele
  if (ix > lat%n_ele_track) then
    print *, 'STOPPING ELEMENT IS A LORD! ', bbu_param%ele_track_end
    call err_exit
  endif
  bbu_param%ix_ele_track_end = ix
endif

!! Record cavity info in hom_info.txt
if (bbu_param%write_hom_info) then
  call rf_cav_names(lat)
endif

call check_rf_freq(lat, bbu_param%bunch_freq)

call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)
print *, 'bbu_setup complete !!!'

beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

print '(a, 2i10)', 'Number of stages and elements in the tracking lattice: ' , size(bbu_beam%stage),  lat%n_ele_track

lat = lat0 ! Restore lr wakes

call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
print *, 'bbu_track_all complete !!!'
 
print *, 'HOM VOLT GAIN: ', hom_voltage_gain
print *, 'growth_rate: ', growth_rate

! Output the BBU results, store them in "for_py.txt" to be analyzed by Python
o = lunget() 
open(o, file = 'for_py.txt', status = 'unknown')
write(o,'(2a)') 'lostbool = ', logical_to_python(lost)  
write(o,'(a, es18.8E3)') 'v_gain = ', hom_voltage_gain
write(o,'(a,es14.6)') 'bunch_dt = ', beam_init%dt_bunch
write(o,'(2a)') 'growth_rate_set = ', logical_to_python( .NOT.(growth_rate == real_garbage$))
write(o,'(a, es14.6)') 'growth_rate = ', growth_rate
close(o)
 
end program
