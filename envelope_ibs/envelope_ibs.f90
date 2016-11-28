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
type(ele_struct) ele_save, ele_at_s
type(ele_struct), allocatable :: eles(:)
type(coord_struct), allocatable :: co(:)
type(coord_struct) co_at_s
type(coord_struct), allocatable :: coos(:)
type(normal_modes_struct) mode

integer i,j,ix,tix
integer nturns, nslices
integer status

logical err_flag

real(rp) mat6(6,6), vec0(6)
real(rp) Sigma_ent(6,6), Sigma_exit(6,6)
real(rp) normal(3)
real(rp) s, delta_s
real(rp) starting_a_emit, starting_b_emit, starting_c_emit
real(rp), allocatable :: M(:,:,:), Bbar(:,:,:)

namelist /envelope_tracker/ lat_file, starting_a_emit, starting_b_emit, starting_c_emit, nturns, nslices

!set defaults
starting_a_emit = 1.0e-9
starting_b_emit = 10.0e-12
starting_c_emit = 0.01 * 1.0e-4
nturns = 40000
nslices = 30000

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

delta_s = lat%param%total_length/(nslices*1.0d0)
write(*,'(a,f11.4,a)') "Slice length is ", delta_s, " m"

! Build element slices and damping matrix D and diffusion matrix B for each slice.
co_at_s = co(0)
s=0.0d0
allocate(eles(nslices))
allocate(coos(nslices))
allocate(M(6,6,nslices))
allocate(Bbar(6,6,nslices))
do i=1,nslices
  call twiss_and_track_at_s(lat, s+delta_s, eles(i), co)
  call mat6_from_s_to_s(lat, mat6, vec0, s, s+delta_s, co_at_s, ele_save=ele_save)
  eles(i)%value(l$) = delta_s
  coos(i) = co_at_s
  eles(i)%mat6 = mat6
  eles(i)%vec0 = vec0
  call make_SR_mats(eles(i),coos(i),M(:,:,i),Bbar(:,:,i))
  s = s + delta_s
enddo

! Track element-slice by element-slice for number of turns.
open(11,file='emit_vs_turn.out')
write(11,'(a8,3a14)') "# turn", "emit_a", "emit_b", "emit_c"
do i=1,nturns
  do j=1,nslices
    call transport_with_sr_ele(eles(j),M(:,:,j),Bbar(:,:,j),Sigma_ent,Sigma_exit) 
    Sigma_ent = Sigma_exit
  enddo
  if( mod(i,10) == 0 ) then
    call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)
    write(11,'(i8,3es14.5)') i, normal(1:3)
  endif
enddo
close(11)

! Extract normal mode emittances from final sigma matrix.
call get_emit_from_sigma_mat(Sigma_exit, normal, err_flag = err_flag)

write(*,'(a,3es14.5)') "Final Emittances: ", normal

open(10,file='sigma_finish.out')
do i=1,6
  write(10,'(6es14.4)') Sigma_exit(i,:)
enddo
close(10)

deallocate(eles)
deallocate(coos)

end program












