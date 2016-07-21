program bbu_program

use bbu_track_mod
use write_lat_file_mod

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in, lat0
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele
type (wake_lr_mode_struct), pointer :: lr(:)

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

print *,'Reading input file BBU.INIT'
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


call run_timer ('START')

!print *, 'lat2 file name is:', bbu_param%lat2_filename
!For DR-scan, parse additional lattice (lat2) 
if (bbu_param%lat2_filename /= '') then
  print *, 'DR-scan or Phase-scan, parsing: ',bbu_param%lat2_filename
  call bmad_parser2 (bbu_param%lat2_filename, lat_in)
endif

!! Closed orbit computed
call twiss_and_track(lat_in,orb,ok,0,.true.)

! Remove HOMs of higher order
if (bbu_param%hom_order_cutoff > 0) then
  do i = 1, lat_in%n_ele_max
    ele => lat_in%ele(i)
    ! Find cavity element with lr_wake
    if (.not. associated(ele%wake)) cycle
    if (.not. allocated(ele%wake%lr_mode)) cycle
    n = count(ele%wake%lr_mode(:)%m > bbu_param%hom_order_cutoff)
    if (n == 0) cycle  !All HOMs order <= cutoff m, nothing to remove
    if (n == size(ele%wake%lr_mode)) then  ! All HOMs order > m, remove this lcavity  
      deallocate (ele%wake%lr_mode)
      cycle
    endif
    !! If some (not all) HOMs order > m, extract the HOMs with order <= m 
    lr => ele%wake%lr_mode
    nn = size(ele%wake%lr_mode) - n    ! nn = number of HOMs to be kept
    allocate(ele%wake%lr_mode(nn))
    n = 0
    do j = 1, size(lr)
      if (lr(j)%m > bbu_param%hom_order_cutoff) cycle
      n = n + 1; ele%wake%lr_mode(n) = lr(j)
    enddo
    deallocate(lr)
  enddo
endif

!!!!!!!!!!!!!!!!!! HYBRIDIZATION !!!!!!!!!!!!!!!!!!!!!!!!
!! ele%select == true means the element will be kept, NOT hybridized 
if (bbu_param%hybridize) then
!  print *, 'CANT HYBRIDIZE -- SET THIS TO FALSE'
!  print *, 'WILL TRY ANYWAY'
!  call err_exit
  do i = 1, lat_in%n_ele_max
    ele => lat_in%ele(i)
    ele%select = .false.
    ! Keep the element at the end of tracking, if specified by user
    if (ele%name == bbu_param%ele_track_end) then
      ele%select = .true.
      cycle
    endif

    !! Specify the element (names) to be kept from hybridization
    !! Avoid choosing the patch elements 

    !if(ele%name=='taylorW') then
    !  ele%select = .true.
    !  cycle
    !endif    
    
    if(ele%key==8) then !! Taylor
      ele%select = .true.
      cycle
    endif    
    
    ! Any non-cavity elements are hybridized (unless user gives exemption above) 
    ! If user chooses not to keep all cavities,
    ! check and keep cavities with lr_wake
    
    if (ele%key /= lcavity$) cycle
    if (.not. bbu_param%keep_all_lcavities) then
      if (.not. associated (ele%wake)) cycle
      if (size(ele%wake%lr_mode) == 0) cycle
    endif
    ele%select = .true.
  enddo
  call make_hybrid_lat (lat_in, lat, bbu_param%use_taylor_for_hybrids)
  print *, 'Hybridization complete !!!'
else
  ! If not hybridizing, keep the original lattice
  lat = lat_in
endif

!!!!!!!!!!!!!!!!!!!!!  END OF HYBRIDIZATION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!   Helpful functions to investigate the (hybridized) lattice !!!!!!!!!!!!!!!
!!!!!!!!   Only use these for testing purpose, since they can slow down or even stop the program !!!!!!!!!! 

!!! Print the element names of the hybridized lattice 
!if (bbu_param%hybridize) then
!  do i = 1, lat%n_ele_track
!     print *, lat%ele(i)%name
!  enddo
!!!endif

!!! Output the lattice file for the  hybridize lattice
!call write_bmad_lattice_file ('/home/wl528/nfs/linux_lib/bsim/bbu/h0.dat', lat)
!call write_bmad_lattice_file ('/home/wl528/nfs/linux_lib/bsim/bbu/pscan.dat', lat)

!!! Output properties of  hybridized elements, or the mat6s between all elements
!!! of the hybridized lattice
!! This line will stop the program during its 2nd run
!open(newunit = file_unit, file = '/home/wl528/nfs/linux_lib/bsim/bbu/mat6.dat', status = "new", action= "write")
!! This line will NOT stop the program. the file will be overwritten every BBU run
!open(newunit = file_unit, file = '/home/wl528/nfs/linux_lib/bsim/bbu/mat6.dat')
! do k = 1, lat%n_ele_track
!   if (lat%ele(k)%key == 16) then     !!key=16 means the element is a hybrid 
!     write(file_unit,"(A20)") "haha"
!     write(file_unit, '(es18.8E2)') lat%ele(k)%value(L$)
!     write(file_unit, '(es18.8E2)') lat%ele(k)%value(E_TOT_START$)
!     write(file_unit, '(es18.8E2)') lat%ele(k)%value(DELTA_E$)
!     write(file_unit, '(es18.8E2)') lat%ele(k)%value(delta_ref_time$)
!
!   endif
!   !do i =1,6
!   !  write(file_unit,'(6F14.7,1es18.8)')(lat%ele(k)%mat6(i,j),j=1,6), lat%ele(k)%vec0(i)
   !enddo
 !enddo  
!close (file_unit)



! Keep the lattice ready to use?
lat0 = lat 

! Define element at which tracking ends, if user didn't specify one
if (bbu_param%ele_track_end.ne.' ') then
  call lat_ele_locator(bbu_param%ele_track_end,lat, eles, n_loc, err)
  ! eles: (output) Ele_pointer_struct, allocatable: Array of matching elements
  ! n_loc: (output) Number of locations found
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

!! Record cavity names in hom_info.txt
!! hom_info.txt can be useful if the user intends to assign HOM files randomly
!! to the cavities
!if (bbu_param%write_hom_info) then
!  call rf_cav_names (lat)
!endif

!print *, 'bbu_setup running...'
call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)
print *, 'bbu_setup complete !!!'

!! This computes n_ele, which is not used at all?
!n_ele = 0
!do i = 1, size(bbu_beam%stage)
!  j = bbu_beam%stage(i)%ix_ele_lr_wake ! j = index of THE lr_wake in stage i
!  call multipass_chain (lat%ele(j), ix_pass, n_links = n)
!  ! ix_pass  -- Integer: Multipass pass number of the input element ( can be an lr wake ). 
!  !                        -1 if input element is not in a multipass section.
!  ! n_links  -- Integer: Number of times the physical element is  passed through.
!  !print *, ix_pass, n
!  
!  if (ix_pass /= 1 .and. n /= 0) cycle
!  n_ele = n_ele + 1
!enddo
!  print *, n_ele

beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

print *, 'Number of stages and elements in the (hybridized) lattice: ' &
, size(bbu_beam%stage),  lat%n_ele_track

lat = lat0 ! Restore lr wakes

!print *, 'bbu_track_all running...'
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
print *, 'bbu_track_all complete !!!'
 
!! these are not accurate
!! These only tell the greatest hom voltage at the end period of tracking
!print *, 'Target cavity number with max volt:', bbu_beam%ix_stage_voltage_max                             ! Target cavity (stage)
!print *, 'Target hom wake number with max volt:', bbu_beam%stage(bbu_beam%ix_stage_voltage_max)%ix_hom_max  ! Target hom wake
!print *,  bbu_beam%hom_voltage_max  

!print *, 'LostBool:', lost                                            
print *, 'HOM VOLT GAIN: ', hom_voltage_gain
print *, 'growth_rate: ', growth_rate
! Output the BBU results, store them in "for_py.txt" to be analyzed by Python
o = lunget() 
open(o, file = 'for_py.txt', status = 'unknown')
write(o,'(2a)') 'lostbool = ', logical_to_python(lost)  
write(o,'(a, es18.8E3)') 'v_gain = ', hom_voltage_gain
write(o,'(a,es14.6)') 'rel_tol = ', bbu_param%rel_tol 
write(o,'(a,es14.6)') 'bunch_dt = ', beam_init%dt_bunch
write(o,'(2a)') 'growth_rate_set = ', logical_to_python( .NOT.(growth_rate == real_garbage$))
write(o,'(a, es14.6)') 'growth_rate = ', growth_rate
close(o)
 
call run_timer ('STOP', time)

end program
