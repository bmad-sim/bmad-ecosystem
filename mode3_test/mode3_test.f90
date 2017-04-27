program one_turn_mat

use bmad
use mode3_mod

implicit none

real(rp) tmat(6,6)
real(rp) N(6,6)
real(rp) tunes_out(3)
real(rp) abz_tunes(3)

integer ix
integer i, k, status

logical err_flag

character(6) :: tmat_file = 't6.dat'

open(10,file=tmat_file)
do i=1,6
  read(10,*) tmat(i,:)
enddo
close(10)

call make_N(tmat,N,err_flag,tunes_out=tunes_out)

open(20,file='output.now')
write(20,'(a,6es17.8)') '"make_N:ord_by_dom_column_sums"   REL   1.0E-6 ', (sum(N(:,k)),k=1,6) 
write(20,'(a,3es17.8)') '"make_N:ord_by_dom_tunes"         REL   1.0E-6 ', tunes_out

abz_tunes = tunes_out((/3,1,2/))

call make_N(tmat,N,err_flag,abz_tunes=abz_tunes,tunes_out=tunes_out)
write(20,'(a,6es17.8)') '"make_N:ord_by_tune_column_sums"  REL   1.0E-6 ', (sum(N(:,k)),k=1,6) 
write(20,'(a,3es17.8)') '"make_N:ord_by_tune_tunes"        REL   1.0E-6 ', tunes_out

close(20)

end program



