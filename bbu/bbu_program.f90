program bbu_program

use bbu_track_mod

implicit none

type (ele_pointer_struct), allocatable :: eles(:)
type (bbu_beam_struct) bbu_beam
type (bbu_param_struct) bbu_param
type (lat_struct) lat, lat_in, lat0
type (beam_init_struct) beam_init
type (ele_struct), pointer :: ele
type (wake_lr_struct), pointer :: lr(:)

integer, allocatable :: ix_out(:)
integer i, ix, j, n, nn, n_ele, ix_pass, o
integer irep
integer n_loc
real(rp) dr
real(rp) time
real(rp) hom_voltage_gain
real(rp) growth_rate
logical, allocatable :: keep_ele(:)
logical err
logical lost

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
call ran_seed_put (bbu_param%ran_seed)

if (bbu_param%ran_gauss_sigma_cut > 0) then
  call ran_gauss_converter (set_sigma_cut = bbu_param%ran_gauss_sigma_cut)
endif

! Init and parse
print *, 'Lattice file: ', trim(bbu_param%lat_filename)
call bmad_parser (bbu_param%lat_filename, lat_in)
call twiss_propagate_all (lat_in)
call run_timer ('START')

!Parse additional settings for drscan
if (bbu_param%lat2_filename /= '') then
  print *, 'Parsing: ',bbu_param%lat2_filename
  call bmad_parser2 (bbu_param%lat2_filename, lat_in)
endif

! Remove HOMs of higher order
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

if (bbu_param%hybridize) then
  print *, 'Note: Hybridizing lattice...'
  allocate (keep_ele(lat_in%n_ele_max))
  allocate (ix_out(lat_in%n_ele_max))
  keep_ele = .false.
  do i = 1, lat_in%n_ele_max
    ! Keep element if defined as end of tracking
    if(lat_in%ele(i)%name == bbu_param%ele_track_end) then
       call update_hybrid_list (lat_in, i, keep_ele, bbu_param%keep_overlays_and_groups)
       cycle
    endif
    if (lat_in%ele(i)%key /= lcavity$) cycle
    if (.not. bbu_param%keep_all_lcavities) then
      if (.not. associated (lat_in%ele(i)%wake)) cycle
      if (size(lat_in%ele(i)%wake%lr) == 0) cycle
    endif
    call update_hybrid_list (lat_in, i, keep_ele, bbu_param%keep_overlays_and_groups)
  enddo
  call make_hybrid_lat (lat_in, keep_ele, .false., lat, ix_out, bbu_param%use_taylor_for_hybrids)
  deallocate (keep_ele)
  deallocate (ix_out)
else
  lat = lat_in
endif

lat0 = lat

! Define element at which tracking ends
if (bbu_param%ele_track_end.ne.' ') then
  call lat_ele_locator(bbu_param%ele_track_end,lat, eles, n_loc, err)
  if(err)call err_exit
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

if (bbu_param%write_hom_info) then
   call write_homs(lat, bbu_param%bunch_freq)
endif

call bbu_setup (lat, beam_init%dt_bunch, bbu_param, bbu_beam)


n_ele = 0
do i = 1, size(bbu_beam%stage)
  j = bbu_beam%stage(i)%ix_ele_lr_wake
  call multipass_chain (lat%ele(j), ix_pass, n_links = n)
  if (ix_pass /= 1 .and. n /= 0) cycle
  n_ele = n_ele + 1
enddo

beam_init%bunch_charge = bbu_param%current * beam_init%dt_bunch

print *, 'Number of lr wake elements in tracking lattice:', size(bbu_beam%stage)

print *, 'Number of elements in lattice:      ', lat%n_ele_track
lat = lat0 ! Restore lr wakes
call bbu_track_all (lat, bbu_beam, bbu_param, beam_init, hom_voltage_gain, growth_rate, lost, irep)
o = lunget() 
open(o, file = 'for_py.txt', status = 'unknown')
write(o,'(2a)') 'lostbool = ', logical_to_python(lost)  
write(o,'(a, es14.6)') 'v_gain = ', hom_voltage_gain
write(o,'(a,es14.6)') 'rel_tol = ', bbu_param%rel_tol 
write(o,'(a,es14.6)') 'bunch_dt = ', beam_init%dt_bunch
write(o,'(2a)') 'growth_rate_set = ', logical_to_python( .NOT.(growth_rate == real_garbage$))
write(o,'(a, es14.6)') 'growth_rate = ', growth_rate
close(o)
 
call run_timer ('STOP', time)

end program
