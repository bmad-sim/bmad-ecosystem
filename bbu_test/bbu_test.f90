program bbu_test

use bbu_track_mod
!use bmad

!use beam_mod

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in, lat0
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele
type (wake_lr_struct), pointer :: lr(:)

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
integer k
logical :: ok

type (coord_struct), allocatable :: orb(:) 
namelist / bbu_params / bbu_param, beam_init, bmad_com


! Defaults for namelist
beam_init%n_particle = 1

print *,'BBU regression test begins. Reading input file BBU.INIT...'
print *
open (1, file = 'bbu.init', status = 'old')
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
call bmad_parser (bbu_param%lat_filename, lat_in) !! lat_in is the parsed lattice



!For DR-scan, parse additional lattice (lat2) 
if (bbu_param%lat2_filename /= '') then
  print *, 'DR-scan or Phase-scan, parsing: ',bbu_param%lat2_filename
  call bmad_parser2 (bbu_param%lat2_filename, lat_in)
endif

!! Closed orbit computed
call twiss_and_track(lat_in,orb,ok,0,.true.)

!!!!!!!!!!!!!!!!!! HYBRIDIZATION !!!!!!!!!!!!!!!!!!!!!!!!
!! ele%select == true means the element will be kept, NOT hybridized 
!if (bbu_param%hybridize) then
!  do i = 1, lat_in%n_ele_max
!    ele => lat_in%ele(i)
!    ele%select = .false.
!    ! Keep the element at the end of tracking, if specified by user
!    if (ele%name == bbu_param%ele_track_end) then
!      ele%select = .true.
!      cycle
!    endif
!    
!    if(ele%key==8) then !! Taylor
!      ele%select = .true.
!      cycle
!    endif    
!    
!    ! Any non-cavity elements are hybridized (unless user gives exemption above) 
!    ! If user chooses not to keep all cavities,
!    ! check and keep cavities with lr_wake
!    
!    if (ele%key /= lcavity$) cycle
!    if (.not. bbu_param%keep_all_lcavities) then
!      if (.not. associated (ele%wake)) cycle
!      if (size(ele%wake%lr) == 0) cycle
!    endif
!    ele%select = .true.
!  enddo
!  call make_hybrid_lat (lat_in, lat, bbu_param%use_taylor_for_hybrids)
!  print *, 'Hybridization complete !!!'
!else
  ! If not hybridizing, keep the original lattice
!  lat = lat_in
!endif

lat = lat_in  !! Comment this out if hybridization is applied.

!!!!!!!!!!!!!!!!!!!!!  END OF HYBRIDIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Keep the lattice ready to use?
lat0 = lat 

call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)
print *, 'bbu_setup complete !!!'

lat = lat0 ! Restore lr wakes

!print *, 'Number of stages and elements in the (hybridized) lattice: ', size(bbu_beam%stage),  lat%n_ele_track

open (2, file = 'output.now', recl = 200)


beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

!print *, 'bbu_track_all running...'
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
!print *, 'bbu_track_all complete !!!'

write (2, '(a, i3)') '"lost_boolean_1A"      ABS 0', lost  
write (2, '(a, es20.12E3)') '"hom_voltage_gain_1A"  ABS 1E-9', hom_voltage_gain
write (2, '(a, es20.12E3)') '"growth_rate_1A"       ABS 1E-9', growth_rate

!print *, 'LOST', logical_to_python(lost)
!print *, 'HOM VOLT GAIN: ', hom_voltage_gain
!print *, 'growth_rate: ', growth_rate
 
!bbu_track_mod_mp_bbu_track_all_param%current = 0.001
bbu_param%current = 0.001
beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

!print *, 'bbu_track_all running...'
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)

write (2, '(a, i3)') '"lost_boolean_1mA"      ABS 0', lost  
write (2, '(a, es20.12E3)') '"hom_voltage_gain_1mA"  ABS 1E-9', hom_voltage_gain
write (2, '(a, es20.12E3)') '"growth_rate_1mA"       ABS 1E-9', growth_rate
!print *, 'LOST', logical_to_python(lost)
!print *, 'HOM VOLT GAIN: ', hom_voltage_gain
!print *, 'growth_rate: ', growth_rate


bbu_param%current = 100
beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

!print *, 'bbu_track_all running...'
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)

write (2, '(a, i3)') '"lost_boolean_100A"      ABS 0', lost  
write (2, '(a, es20.12E3)') '"hom_voltage_gain_100A"  ABS 1E-9', hom_voltage_gain
write (2, '(a, es20.12E3)') '"growth_rate_100A"       ABS 1E-9', growth_rate
!print *, 'LOST', logical_to_python(lost)
!print *, 'HOM VOLT GAIN: ', hom_voltage_gain
!print *, 'growth_rate: ', growth_rate


close(2)

end program
