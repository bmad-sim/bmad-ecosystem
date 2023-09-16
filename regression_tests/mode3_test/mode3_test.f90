program one_turn_mat

use bmad
use mode3_mod

implicit none

type (lat_struct) lat
type (ele_pointer_struct), allocatable :: eles(:)
type (coord_struct), allocatable :: orbit(:)

real(rp) tmat(6,6)
real(rp) N_mat(6,6)
real(rp) tunes_out(3)
real(rp) abz_tunes(3)
real(rp), allocatable :: dk1(:)

integer ix
integer i, k, n, status

logical err_flag, ok

character(6) :: tmat_file = 't6.dat'

open(10,file=tmat_file)
do i=1,6
  read(10,*) tmat(i,:)
enddo
close(10)

call make_N(tmat,N_mat,err_flag,tunes_out=tunes_out)

open(20,file='output.now')
write(20, '(a,6es17.8)') '"make_N:ord_by_dom_column_sums"   REL   1.0E-6 ', (sum(N_mat(:,k)),k=1,6) 
write(20, '(a,3es17.8)') '"make_N:ord_by_dom_tunes"         REL   1.0E-6 ', tunes_out

abz_tunes = tunes_out((/3,1,2/))

call make_N(tmat,N_mat,err_flag,abz_tunes=abz_tunes,tunes_out=tunes_out)
write(20, '(a,6es17.8)') '"make_N:ord_by_tune_column_sums"  REL   1.0E-6 ', (sum(N_mat(:,k)),k=1,6) 
write(20,' (a,3es17.8)') '"make_N:ord_by_tune_tunes"        REL   1.0E-6 ', tunes_out

!

call bmad_parser ('small_ring.bmad', lat)
call twiss_and_track(lat, orbit)
call choose_quads_for_set_tune(lat%branch(0), dk1, eles)
ok = set_tune(2.4_rp*twopi, 1.7_rp*twopi, dk1, eles, lat%branch(0), orbit)
call set_z_tune(lat%branch(0), -0.2_rp)

n = lat%n_ele_track
write (20, '(a, 2es17.8)') '"tune"     ABS  1.0E-6', lat%ele(n)%a%phi/twopi, lat%ele(n)%b%phi/twopi
write (20, '(a, 2es17.8)') '"z_tune"   ABS  1.0E-6', lat%z%tune

!

close(20)

end program



