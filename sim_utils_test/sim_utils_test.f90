!+
! Regression tests for simulation utility routines
!-

program math_test

use bmad
use random_mod
use nr
use naff_mod

implicit none

type (coord_struct) orbit

real(rp) array(4), dE, freq(3)
complex(rp) cdata(32)
complex(rp) amp(3)
integer i, which, where, n_freq
logical match
character(40) str, sub1, sub2, sub3

real(rp) sig1, sig2, sig3
real(rp) phi1, phi2, phi3
real(rp) vec3(3), vec3a(3), vec3b(3), axis(3), angle, w_mat(3,3)
complex(rp) amp1, amp2, amp3

!

open (1, file = 'output.now')

! rotation test

axis = [3, 4, 5] / sqrt(50.0_rp)
vec3 = [-2, 3, -4]
angle = 0.67

!axis = [1, 0, 0]
!vec3 = [-2, 3, -4]
!angle = pi/2

vec3a = vec3
call rotate_vec_given_axis_angle (vec3a, axis, angle)

call axis_angle_to_w_mat (axis, angle, w_mat)
vec3b = matmul(w_mat, vec3)

write (1, '(a, 3f11.6)') '"rot vecA" ABS 1E-10  ', vec3a
write (1, '(a, 3f11.6)') '"rot vecB" ABS 1E-10  ', vec3b
write (1, '(a, 3es10.2)') '"drot vec" ABS 1E-10  ', vec3a - vec3b

! naff test

amp1 = cmplx(1.8000,0.0000)
sig1 = 0.753262
amp2 = cmplx(0.3000,0.3000)
sig2 = 0.423594
amp3 = cmplx(0.01230,0.1545)
sig3 = 0.173

do i = 1, size(cdata)
  phi1 = twopi*(sig1*(i-1))
  phi2 = twopi*(sig2*(i-1))
  phi3 = twopi*(sig3*(i-1))
  cdata(i) = amp1*exp(cmplx(0.0d0,-phi1)) + amp2*exp(cmplx(0.0d0,-phi2)) + amp3*exp(cmplx(0.0d0,-phi3))
enddo

call naff (cdata, freq, amp)
write (1, '(a, 3es16.8)') '"naff-freq1" REL 2E-6   ', freq(1), real(amp(1)), aimag(amp(1))
write (1, '(a, 3es16.8)') '"naff-freq2" REL 2E-6   ', freq(2), real(amp(2)), aimag(amp(2))
write (1, '(a, 3es16.8)') '"naff-freq3" REL 3E-6   ', freq(3), real(amp(3)), aimag(amp(3))

! Random test

call ran_engine ('quasi')

do i = 1, 10
  call ran_uniform_vector (array)
enddo

write (1, '(a, 4es20.10)') '"QuasiRan" ABS  0', array

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
