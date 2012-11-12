!+
! Regression tests for math routines
!-

program math_test

use sim_utils
use random_mod
use nr

implicit none

real(rp) array(4)
integer i

!

call ran_engine ('quasi')

do i = 1, 10
  call ran_uniform_vector (array)
enddo

open (1, file = 'output.now')
write (1, '(a, 4es20.10)') '"Quasi" ABS  0', array
close(1)

end program
