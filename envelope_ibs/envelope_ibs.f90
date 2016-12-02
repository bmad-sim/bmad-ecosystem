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

character lat_file*200
character in_file*200
character ix_str*5

type(lat_struct) lat
type(ele_struct) temp_ele
type(ele_struct), allocatable :: eles(:)
type(coord_struct), allocatable :: co(:)
type(coord_struct) co_at_s
type(coord_struct), allocatable :: coos(:)
type(normal_modes_struct) mode

integer i,j
integer nturns, nslices, ns, six
integer status
integer pct_complete, last_display

logical err_flag, include_ibs, tail_cut

real(rp) alpha(3), emit(3)
real(rp) mat6(6,6), vec0(6)
real(rp) Sigma_ent(6,6), Sigma_exit(6,6)
real(rp) normal(3)
real(rp) current, npart, lasts
real(rp) bend_slice_length, slice_length
real(rp) tau(3), tau_max, etay
real(rp) starting_a_emit, starting_b_emit, starting_c_emit
real(rp), allocatable :: M(:,:,:), Bbar(:,:,:)

namelist /envelope_tracker/ lat_file, starting_a_emit, starting_b_emit, starting_c_emit, etay, nturns, include_ibs, current, tail_cut, bend_slice_length, slice_length

!set defaults
etay = -1.0 ! negative number to disable
bend_slice_length = 0.01d0
slice_length = 100.0d0  !large value results in element length being used
starting_a_emit = 1.0d-9
starting_b_emit = 10.0d-12
starting_c_emit = 0.01 * 1.0d-4
nturns = 40000
include_ibs = .false.
tail_cut = .true.
current = 0.001

call getarg(1, in_file)
open (unit = 20, file = in_file, action='read')
read (20, nml = envelope_tracker)
close (20)

bmad_com%radiation_damping_on = .true.
bmad_com%radiation_fluctuations_on = .false.

call bmad_parser(lat_file, lat)
call twiss_and_track(lat,co,status)
call calc_z_tune(lat)
call transfer_matrix_calc(lat, mat6, vec0, ix1=0, one_turn=.true.)

npart = current / e_charge * lat%param%total_length / c_light

!make_smat_from_abc uses mode tunes to label modes
mode%a%tune = lat%a%tune
mode%b%tune = lat%b%tune
mode%z%tune = lat%z%tune
mode%a%emittance = starting_a_emit
mode%b%emittance = starting_b_emit
mode%z%emittance = starting_c_emit
call make_smat_from_abc(mat6, mode, Sigma_ent, err_flag)

write(*,'(a,3es14.5)') "Initial Emittances: ", mode%a%emittance, mode%b%emittance, mode%z%emittance
open(10,file='sigma_start.out')
do i=1,6
  write(10,'(6es14.4)') Sigma_ent(i,:)
enddo
close(10)

! Build element slices and damping matrix D and diffusion matrix B for each slice.
co_at_s = co(0)
nslices = 0
do i=1,lat%n_ele_track
  if(lat%ele(i)%value(l$) .gt. 1.0d-6) then
    if(any(lat%ele(i)%key==(/sbend$,rbend$/))) then
      ns = ceiling(lat%ele(i)%value(l$) / bend_slice_length)
      nslices = nslices+ns
    else
      ns = ceiling(lat%ele(i)%value(l$) / slice_length)
      nslices = nslices+ns
    endif
  endif
enddo
write(*,'(a,i8)') "Number of slices: ", nslices
allocate(eles(nslices))
allocate(coos(nslices))
allocate(M(6,6,nslices))
allocate(Bbar(6,6,nslices))
lasts = 0.0d0
six = 0
do i=1,lat%n_ele_track
  if(lat%ele(i)%value(l$) .gt. 1.0d-6) then
    if(any(lat%ele(i)%key==(/sbend$,rbend$/))) then
      ns = ceiling(lat%ele(i)%value(l$) / bend_slice_length)
    else
      ns = ceiling(lat%ele(i)%value(l$) / slice_length)
    endif
    if(ns .gt. 1) then
      co_at_s = co(i)
      do j=1,ns
        six = six + 1
        call create_uniform_element_slice(lat%ele(i),lat%param,j,ns,temp_ele)
        eles(six) = temp_ele
        call mat6_from_s_to_s(lat, mat6, vec0, lasts, eles(six)%s, co_at_s)
        coos(six) = co_at_s
        eles(six)%mat6 = mat6
        eles(six)%vec0 = vec0
        call make_SR_mats(eles(six),coos(six),etay,M(:,:,six),Bbar(:,:,six))
        lasts = eles(six)%s
      enddo
    else
      six = six + 1
      coos(six) = co(i)
      eles(six)%value(E_TOT$) = lat%ele(i)%value(E_TOT$)
      eles(six)%value(l$) = lat%ele(i)%value(l$)
      eles(six)%s = lasts + eles(six)%value(l$)
      eles(six)%mat6 = lat%ele(i)%mat6
      eles(six)%vec0 = lat%ele(i)%vec0
      call make_SR_mats(eles(six),coos(six),etay,M(:,:,six),Bbar(:,:,six))
      lasts = lat%ele(i)%s
    endif
  endif
enddo

call envelope_radints(eles,coos,alpha,emit)
tau = lat%param%total_length / c_light / alpha
write(*,*) "From envelope-based radiation integrals calculations: "
write(*,'(a,3es14.5)') "  Emittance:                  ", emit
write(*,'(a,3es14.5)') "  Damping decrements (alpha): ", alpha
write(*,'(a,3es14.5)') "  Damping times:              ", tau
write(*,*) "Starting tracking..."

tau_max = maxval(tau)

! Track element-slice by element-slice for number of turns.
open(11,file='emit_vs_turn.out')
write(11,'(a8,3a14)') "# turn", "emit_a", "emit_b", "emit_c"
last_display = -1.0
do i=1,nturns
  do j=1,nslices
    if(include_ibs) then
      call transport_with_sr_and_ibs(eles(j),M(:,:,j),Bbar(:,:,j),Sigma_ent,Sigma_exit,tail_cut,tau_max,npart)
    else
      call transport_with_sr(eles(j),M(:,:,j),Bbar(:,:,j),Sigma_ent,Sigma_exit) 
    endif
    Sigma_ent = Sigma_exit
  enddo
  if( mod(i,10) == 0 ) then
    call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)
    write(11,'(i8,3es14.5)') i, normal(1:3)
  endif
  pct_complete = floor((i*100.0d0)/nturns)
  if( pct_complete .gt. last_display ) then
    write(*,'(a,i3,a,$)') " ...",pct_complete,"%"
    last_display = pct_complete
  endif
enddo
close(11)

! Extract normal mode emittances from final sigma matrix.
call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)

write(*,*)
write(*,'(a,3es14.5)') "  Final Emittances: ", normal

open(10,file='sigma_finish.out')
do i=1,6
  write(10,'(6es14.4)') Sigma_exit(i,:)
enddo
close(10)

deallocate(eles)
deallocate(coos)

end program












