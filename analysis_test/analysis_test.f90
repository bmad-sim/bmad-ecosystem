program analysis_test

use bmad
implicit none

type (coord_struct), allocatable :: track(:)
type (coord_struct) track0

real(rp) map1(6,6), vec0(6), chi, m_in(6,6), v_const(6), v_osc(6), v0(6)
integer i

!

open (1, file = 'output.now')

!

m_in(1,:) = [-1.82011930,  10.42829553,   0.26136852,   1.81672992,   0.11951783,   3.31755608]
m_in(2,:) = [-0.29169861,   1.17887090,  -0.00175668,   0.29248149,   0.01872841,   0.72175051]
m_in(3,:) = [-0.13133050,   1.65002161,  -1.12832722,  -1.36793276,  -0.01614562,  -0.93887564]
m_in(4,:) = [-0.06126760,   0.17321056,   0.09119966,  -0.70837788,   0.00262819,   0.30036736]
m_in(5,:) = [ 0.30041831,  -3.77149768,  -0.01739667,   0.73226876,   0.64433077, -15.20374345]
m_in(6,:) = [ 0.01636467,  -0.12082477,   0.00460790,   0.01722688,   0.03695975,   0.64613703]

allocate(track(100))
v_osc  = [0.001, 0.002, 0.003, 0.004, 0.005, 0.006]
v_const = [0.006, 0.005, 0.004, 0.003, 0.002, 0.001]
v0 = v_osc

do i = 1, 100
  track(i)%vec = v_osc + v_const
  v_osc = matmul(m_in, v_osc)
enddo


call multi_turn_tracking_to_mat (track, 6, map1, vec0, track0, chi)

do i = 1, 6
  write (1, '(a, i0, a, 6es12.4)') '"mt_mat', i, '" ABS 2E-13', map1(i,:) - m_in(i,:)
enddo

write (1, '(a, 6es12.4)') '"mt_vec" ABS 1E-13', track0%vec - v_const

end program
