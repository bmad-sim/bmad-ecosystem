!+
! Regression tests for simulation utility routines
!-

program math_test

use bmad
use str_find_first_substring_module
use random_mod
use nr

implicit none

type (coord_struct) orbit
real(rp) array(4), dE
integer i, which, where
logical match
character(40) str, sub1, sub2, sub3

!

open (1, file = 'output.now')

! Random test

call ran_engine ('quasi')

do i = 1, 10
  call ran_uniform_vector (array)
enddo

write (1, '(a, 4es20.10)') '"QuasiRan" ABS  0', array

! str matching test

str = ' Hello world testing s3 '
sub1 = 's1' ; sub2 = 's2' ; sub3 = 's3'
match = str_find_first_substring(str, where, which, sub1, sub2, sub3)

write (1, *)
write (1, '(a, i0)') '"Which-Find" ABS 0   ', which
write (1, '(a, i0)') '"Where-Find" ABS 0   ', where

! 

write (1, *)

orbit%p0c = 1e6
orbit%vec = 0
orbit%vec(6) = 0.1
orbit%species = positron$

call convert_pc_to (orbit%p0c * (1 + orbit%vec(6)), positron$, beta = orbit%beta)
call apply_energy_kick (1d2, orbit)
write (1, '(a, 2es20.12)') '"apply_energy_kick:0" REL 1E-12  ', orbit%beta, orbit%vec(6)

call convert_pc_to (orbit%p0c * (1 + orbit%vec(6)), positron$, beta = orbit%beta)
call apply_energy_kick (1d6, orbit)
write (1, '(a, 2es20.12)') '"apply_energy_kick:1" REL 1E-12  ', orbit%beta, orbit%vec(6)

!

close(1)

end program
