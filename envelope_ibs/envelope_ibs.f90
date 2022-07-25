!+
! This is a front end for bmad/multiparticle/envelope_mod.f90
!
! The module and this front end currently implement the beam envelope matrix
! based synchrotron radiation and damping calculation contained in
! "From the beam-envelope matrix to synchrotron-radiation integrals" by
! K. Ohmi, K. Hirata, and K. Oide.
!-

program envelope_ibs

use bmad
use transfer_map_mod
use mode3_mod
use envelope_mod

implicit none

character(200) lat_file, in_file
character(5) ix_str

type(lat_struct) lat
type(ele_struct) temp_ele
type(ele_struct), allocatable :: eles(:)
type(coord_struct), allocatable :: co(:)
type(coord_struct), allocatable :: coos(:)
type(normal_modes_struct) mode

integer i,j
integer nturns, nslices, ns, six
integer status

logical err_flag, include_ibs, tail_cut, user_supplied_tunes, regression_test

real(rp) tune_x, tune_y, tune_z
real(rp) one_turn_mat(6,6), one_turn_vec(6)
real(rp) one_turn_check(6,6), vec_check(6)
real(rp) alpha(3), emit(3)
real(rp) mat6(6,6), vec0(6)
real(rp) Sigma_ent(6,6), Sigma_exit(6,6)
real(rp) normal(3)
real(rp) current, npart, now_time, last_time
real(rp) bend_slice_length, slice_length
real(rp) tau(3), tau_max, ey0, Ykick_strength
real(rp) starting_a_emit, starting_b_emit, starting_c_emit
real(rp), allocatable :: M(:,:,:), Bbar(:,:,:), Ybar(:,:,:)

complex(rp) Lambda(6,6), Theta(6,6), Iota_base(6,6), Iota(6,6)
complex(rp) eval(6), evec(6,6)

namelist /envelope_tracker/ lat_file, starting_a_emit, starting_b_emit, starting_c_emit, ey0, nturns, include_ibs, &
                            current, tail_cut, bend_slice_length, slice_length, user_supplied_tunes, tune_x, &
                            tune_y, tune_z, regression_test

!set defaults
user_supplied_tunes = .false.
ey0 = -1.0 !negative number to disable
bend_slice_length = 100.0d0
slice_length = 100.0d0  !large value results in element length being used
starting_a_emit = 1.0d-9
starting_b_emit = 10.0d-12
starting_c_emit = 0.01 * 1.0d-4
nturns = 40000
include_ibs = .false.
tail_cut = .true.
current = 0.001
regression_test = .false.

call getarg(1, in_file)
if (in_file == '') in_file = 'envelope_ibs.init'
open (unit = 20, file = in_file, action='read')
read (20, nml = envelope_tracker)
close (20)

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call bmad_parser(lat_file, lat)
call twiss_and_track(lat,co,status)
call calc_z_tune(lat%branch(0))

npart = current / e_charge * lat%param%total_length / c_light

! Make_smat_from_abc uses mode tunes to label modes

if (user_supplied_tunes) then
  mode%a%tune = twopi * tune_x
  mode%b%tune = twopi * tune_y
  mode%z%tune = twopi * tune_z

  if (tune_x == 0 .or. tune_y == 0 .or. tune_z == 0) then
    print *, 'USER_SUPPLIED_TUNES = TRUE BUT NOT ALL TUNE_X, TUNE_Y, OR TUNE_Z SET!'
    stop
  endif

  if (abs(tune_x) > 1 .or. abs(tune_y) > 1 .or. abs(tune_z) > 1) then
    print *, 'UNITS OF FRACTIONAL TUNES TUNE_X, TUNE_Y, AND TUNE_Z ARE NOW RAD/2PI AND MUST BE IN THE INTERVAL (0:1).'
    stop
  endif

else
  call transfer_matrix_calc(lat, mat6, vec0)
  call mat_eigen (mat6, eval, evec, err_flag)

  mode%a%tune = modulo(atan2(aimag(eval(1)), real(eval(1))), twopi)
  mode%b%tune = modulo(atan2(aimag(eval(3)), real(eval(3))), twopi)
  mode%z%tune = modulo(atan2(aimag(eval(5)), real(eval(5))), twopi)
  if (mode%z%tune > pi) mode%z%tune = twopi - mode%z%tune
endif

mode%a%emittance = starting_a_emit
mode%b%emittance = starting_b_emit
mode%z%emittance = starting_c_emit

print '(a, 3f10.6)', 'Tunes: ', mode%a%tune/twopi, mode%b%tune/twopi, mode%z%tune/twopi
print '(a, 3es14.5)', "Initial Emittances: ", mode%a%emittance, mode%b%emittance, mode%z%emittance

call transfer_matrix_calc(lat, mat6, vec0, ix1=0, one_turn=.true.)
call make_smat_from_abc(mat6, mode, Sigma_ent, err_flag)

open(10,file='sigma_start.out')
do i=1,6
  write(10,'(6es14.4)') Sigma_ent(i,:)
enddo
close(10)

! Build element slices and damping matrix D and diffusion matrix B for each slice.

nslices = 0
do i=1,lat%n_ele_track
  if(any(lat%ele(i)%key==(/marker$,monitor$/))) then
    cycle
  else
    if(any(lat%ele(i)%key==(/sbend$,rbend$,wiggler$/))) then
      ns = max(1,ceiling(lat%ele(i)%value(l$) / bend_slice_length))
    else
      ns = max(1,ceiling(lat%ele(i)%value(l$) / slice_length))
    endif
    nslices = nslices + ns
  endif
enddo
write(*,'(a,i8)') "Number of slices: ", nslices
allocate(eles(nslices))
allocate(coos(0:nslices))
allocate(M(6,6,nslices))
allocate(Bbar(6,6,nslices))
allocate(Ybar(6,6,nslices))
six = 0
coos(0)=co(0)
do i=1,lat%n_ele_track
  if(any(lat%ele(i)%key== [marker$,monitor$])) then
    cycle
  else
    if(any(lat%ele(i)%key==(/sbend$,rbend$,wiggler$/))) then
      ns = max(1,ceiling(lat%ele(i)%value(l$) / bend_slice_length))
    else
      ns = max(1,ceiling(lat%ele(i)%value(l$) / slice_length))
    endif
    if(ns .gt. 1) then
      do j=1,ns
        six = six + 1
        call element_slice_iterator(lat%ele(i),lat%param,j,ns,temp_ele)
        call make_mat6(temp_ele, lat%param, coos(six-1), coos(six))
        eles(six) = temp_ele
        call make_SR_mats(eles(six),coos(six),M(:,:,six),Bbar(:,:,six))
        call make_Ykick_mat(eles(six),Ybar(:,:,six))
      enddo
    else
      six = six + 1
      coos(six) = co(i) !co at element exit
      eles(six) = lat%ele(i)
      call make_SR_mats(eles(six),coos(six),M(:,:,six),Bbar(:,:,six))
      call make_Ykick_mat(eles(six),Ybar(:,:,six))
    endif
  endif
enddo

one_turn_mat = I6
one_turn_vec = 0.0d0
do i=1,size(eles)
  call concat_transfer_mat(eles(i)%mat6,eles(i)%vec0,one_turn_mat,one_turn_vec,one_turn_mat,one_turn_vec)
enddo
open(50,file='one_turn_mat_from_slices.diag')
write(50,*) "mat6 from slices:"
do i=1,6
  write(50,'(6es17.8)') one_turn_mat(i,:)
enddo 
write(50,*) "vec0 from slices:"
write(50,'(6es17.8)') one_turn_vec(:)
call transfer_matrix_calc(lat, one_turn_check, vec_check)
write(50,*)
write(50,*) "These matrices should be identical."
write(50,*)
write(50,*) "mat6 from canned routine:"
do i=1,6
  write(50,'(6es17.8)') one_turn_check(i,:)
enddo 
write(50,*) "vec0 from canned routine:"
write(50,'(6es17.8)') vec_check(:)
close(50)

write(*,*) "Calculating integrated matrices ..."
call integrated_mats(eles,coos,Lambda,Theta,Iota_base,mode)
write(*,*) "... done."
Ykick_strength = ey0 * 2.0d0 * real(Lambda(3,3)) / Iota_base(3,3)
Iota = Iota_base * Ykick_strength !for integral calculations
Ybar = Ybar * Ykick_strength !for tracking

write(*,*) "Calculating envelope radiation integrals without IBS ..."
call envelope_radints(Lambda,Theta,Iota,alpha,emit)
tau = lat%param%total_length / c_light / alpha
write(*,*) "... done."
write(*,*) "From envelope-based radiation integrals without IBS: "
write(*,'(a,3es14.5)') "  Emittance (SR only):        ", emit
write(*,'(a,3es14.5)') "  Damping decrements (alpha): ", alpha
write(*,'(a,3es14.5)') "  Damping times:              ", tau

if(include_ibs) then
  write(*,*) "Calculating envelope radiation integrals with IBS ..."
  call envelope_radints_ibs(Lambda,Theta,Iota,eles,alpha,emit,mode,tail_cut,npart, lat%param%particle)
  tau = lat%param%total_length / c_light / alpha
  write(*,*) "... done."
  write(*,*) "From envelope-based radiation integrals with IBS: "
  write(*,'(a,3es14.5)') "  Emittance (SR only):        ", emit
  write(*,'(a,3es14.5)') "  Emittance (with IBS):       ", mode%a%emittance, mode%b%emittance, mode%z%emittance
  write(*,'(a,3es14.5)') "  Damping decrements (alpha): ", alpha
  write(*,'(a,3es14.5)') "  Damping times:              ", tau
endif

open(20,file='emit_from_integrals.out')
write(20,*) "From envelope-based radiation integrals and IBS calculations: "
write(20,'(a,3es14.5)') "  Emittance (SR only):        ", emit
if(include_ibs) write(20,'(a,3es14.5)') "  Emittance (with IBS):       ", mode%a%emittance, mode%b%emittance, mode%z%emittance
write(20,'(a,3es14.5)') "  Damping decrements (alpha): ", alpha
write(20,'(a,3es14.5)') "  Damping times:              ", tau
close(20)

if(nturns .gt. 0) then
  write(*,*) "Starting tracking..."
  tau_max = maxval(tau)
  ! Track element-slice by element-slice for number of turns.
  open(11,file='emit_vs_turn.out')
  write(11,'(a8,3a14)') "# turn", "emit_a", "emit_b", "emit_c"

  last_time = 0
  call run_timer ('START')

  do i=1,nturns
    do j=1,nslices
      if(include_ibs) then
        call transport_with_sr_and_ibs(eles(j),M(:,:,j),Bbar(:,:,j),Ybar(:,:,j),Sigma_ent,Sigma_exit,&
                                                                              tail_cut,tau_max,npart, lat%param%particle)
      else
        call transport_with_sr(eles(j),M(:,:,j),Bbar(:,:,j),Ybar(:,:,j),Sigma_ent,Sigma_exit) 
      endif
      Sigma_ent = Sigma_exit
    enddo
    if( mod(i,10) == 0 ) then
      call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)
      write(11,'(i8,3es14.5)') i, normal(1:3)
    endif

    call run_timer ('READ', now_time)
    if (now_time - last_time > 100) then
      write(*,'(a, f10.2, i8, a, i8, a, i8)') 'Time (min): ', now_time/60.0, i, " turns ", i, ' of ', nturns
      last_time = now_time
    endif
  enddo
  close(11)

  ! Extract normal mode emittances from final sigma matrix.
  call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)

  write(*,*)
  write(*,'(a)') "From tracking beam enevelope: "
  write(*,'(a,3es14.5)') "  Emittance:        ", normal
  write(*,*) "Beam sigma matrix saved as sigma_finish.out"

  open(10,file='sigma_finish.out')
  do i=1,6
    write(10,'(6es14.4)') Sigma_exit(i,:)
  enddo
  close(10)
else
  write(*,*) "Tracking disabled in .in file.  nturns < 0"
endif

deallocate(eles)
deallocate(coos)

! Create output appropriate for Bmad regression test suite.

if (regression_test) then
  open (1, file = 'output.now')

  write (1, '(a, 3es16.8)') '"init-emit" REL 1E-8', mode%a%emittance, mode%b%emittance, mode%z%emittance

  do i = 1, 6
    write (1, '(a, 6es16.8)') '"sigma-ent-' // int_str(i) // '" ABS 1E-12', sigma_ent(i,:)
  enddo

  do i = 1, 6
    write(1, '(a, 6es16.8)') '"one-turn-mat-' // int_str(i) // '" ABS 1e-12', one_turn_mat(i,:)
  enddo 

  write(1, '(a, 3es16.8)') '"emit" REL 1e-8', emit
  write(1, '(a, 3es16.8)') '"alpha" REL 1e-8', alpha
  write(1, '(a, 3es16.8)') '"tao" REL 1e-8', tau
  write(1, '(a, 3es16.8)') '"normal" REL 1e-5', normal

  do i = 1, 6
    write(1, '(a, 6es16.8)') '"sigma-exi-' // int_str(i) // '" ABS 1e-12', sigma_exit(i,:)
  enddo 

  close (1)
endif

end program
