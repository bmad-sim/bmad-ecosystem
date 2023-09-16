program beam_test

use sim_utils

implicit none

real(rp) sig(6,6), sig_s(6,6)
complex(rp) :: eval(6) = 0.0, evec(6,6)

integer i

logical err

!

open (1, file = 'output.now', recl = 200)

sig(1,:) = [ 1.19931384E-05_rp, -3.14623816E-06_rp,  3.50412581E-11_rp, -1.41948045E-11_rp, -1.42521693E-05_rp,  1.37765015E-14_rp]
sig(2,:) = [-3.14623816E-06_rp,  1.34651100E-06_rp,  8.13328978E-12_rp,  3.81984869E-12_rp,  4.08589489E-06_rp,  6.40870364E-15_rp]
sig(3,:) = [ 3.50412581E-11_rp,  8.13328978E-12_rp,  8.84304935E-07_rp, -2.94145634E-07_rp, -3.18540607E-11_rp,  1.44886279E-15_rp]
sig(4,:) = [-1.41948045E-11_rp,  3.81984869E-12_rp, -2.94145634E-07_rp,  1.68526113E-07_rp, -5.80361508E-13_rp,  8.10684695E-16_rp]
sig(5,:) = [-1.42521693E-05_rp,  4.08589489E-06_rp, -3.18540607E-11_rp, -5.80361508E-13_rp,  1.71678462E-05_rp,  5.10790004E-14_rp]
sig(6,:) = [ 1.37765015E-14_rp,  6.40870364E-15_rp,  1.44886279E-15_rp,  8.10684695E-16_rp,  5.10790004E-14_rp,  1.00000000E-14_rp]

sig_s(:,1) = -sig(:,2)
sig_s(:,2) =  sig(:,1)
sig_s(:,3) = -sig(:,4)
sig_s(:,4) =  sig(:,3)
sig_s(:,5) = -sig(:,6)
sig_s(:,6) =  sig(:,5)


call mat_eigen (sig_s, eval, evec, err)

do i = 1, 6
  write (1, '(a, i0, a, 2es16.8)') '"EigenVal', i, '" ABS 1E-20', eval(i)
  write (1, '(a, i0, a, 6f14.10)') '"EigenVec-re', i, '" ABS 1E-6', real(evec(i,:), rp)
  write (1, '(a, i0, a, 6f14.10)') '"EigenVec-im', i, '" ABS 1E-6', aimag(evec(i,:))
  write (1, *)
enddo

close (1)

end program
