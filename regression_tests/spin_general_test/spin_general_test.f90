program spin_general_test

use bmad
use pointer_lattice, only: c_linear_map, operator(*), assignment(=)
use tao_interface

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (ele_struct) t_ele
type (coord_struct) orb0, orb_start, orb_end, orb1, orb2
type (spin_orbit_map1_struct) map1, inv_map1
type (tao_lattice_branch_struct), target :: tao_branch
type (tao_spin_map_struct) :: sm

real(rp) spin_a(3), spin_b(3), spin0(3), dr(6), a_quat(0:3), n_vec(3)
real(rp) mat6(6,6), smap(0:3,0:6), q0(0:3), q1(0:3), t, q, vec(6)
real(rp) xi_sum, xi_diff, n0(3)

complex(rp) orb_eval(6), orb_evec(6,6), spin_evec(6,3)
integer i, j, nargs, n_eigen
logical print_extra, err, err_flag

character(40) :: lat_file = 'spin_general_test.bmad'
character(100) excite_zero(3)

namelist / param / dr

!                  

global_com%exit_on_error = .false.

print_extra = .false.
nargs = command_argument_count() 
if (nargs > 1) then  
  print *, 'Only one command line arg permitted.'
  call err_exit                                  

elseif (nargs > 0)then
  call get_command_argument(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.             
endif                              

!

open (1, file = 'output.now')

!---------------------------------
t = 0.1
n0 = [0.2, 0.3, 0.5]
n0 = n0 / norm2(n0)
map1%spin_q(0:3,0) = [cos(t), sin(t)*n0(1), sin(t)*n0(2), sin(t)*n0(3)]
map1%spin_q(0,1:) = [1, -2, 3, -4, 5, 6]
map1%spin_q(0,1:) = [2, 4, -6, 8, 0, -6]
map1%spin_q(0,1:) = [8, 7, -6, 5, 4, 3]
map1%spin_q(0,1:) = [5, 7, 9, 2, 3, 5]
map1%spin_q(0,1:) = [6, -4, 2, -4, 8, 6]
map1%spin_q(0,1:) = [-7, 1, 9, 2, -4, -3]

call mat_make_unit(map1%orb_mat)
map1%orb_mat(1,2) = 2
map1%orb_mat(4,3) = 3

inv_map1 = map1_inverse(map1)
map1 = map1 * inv_map1

do i = 0, 3
  write (1, '(a, 7f14.8)') '"sq*sqinv:' // int_str(i) // '" ABS 1e-6', map1%spin_q(i,:)
enddo

do i = 1, 6
  write (1, '(a, 6f14.8)') '"m*minv:' // int_str(i) // '" ABS 1e-6', map1%orb_mat(i,:)
enddo

!---------------------------------
! Eigen anal with and without RF.

mat6(1,:) = [-1.45026813625181_rp, 9.65474735831485_rp, -0.31309608791633_rp, -0.171199120220268_rp, 0.0_rp, 3.37037875356424_rp]
mat6(2,:) = [-0.328934364801585_rp, 1.500272137146_rp, -0.0503772356507922_rp, -0.0276056627287384_rp, 0.0_rp, 0.912803122319566_rp]
mat6(3,:) = [-0.114957344655893_rp, -0.0631990019710067_rp, 1.21287214837796_rp, 4.50079805390412_rp, 0.0_rp, 0.000321081798152189_rp]
mat6(4,:) = [0.00316047124266308_rp, 0.00157494979268561_rp, -0.293681197933547_rp, -0.265335840335757_rp, 0.0_rp, -8.76011184832822e-06_rp]
mat6(5,:) = [0.215175896257439_rp, -3.75639824621782_rp, 0.115921051143903_rp, 0.063183785322055_rp, 1.0_rp, -0.929978144660634_rp]
mat6(6,:) = [0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 0.0_rp, 1.0_rp]

smap(0,:) = [0.629042006322742_rp, -0.418211907218498_rp, 0.635632897145065_rp, -0.0567487273391509_rp, -0.0114013118658281_rp, 0.0_rp, 0.466973604342706_rp]
smap(1,:) = [-0.00685250940958695_rp, -0.0220406297379937_rp, -0.0269930570480571_rp, 0.667012158147191_rp, 2.17195796294257_rp, 0.0_rp, -0.00293561427302652_rp]
smap(2,:) = [-0.777294845209733_rp, -0.338517412818692_rp, 0.5143587138526_rp, -0.0476276423794947_rp, -0.0307172897350696_rp, 0.0_rp, 0.377861514185347_rp]
smap(3,:) = [-0.00848062536777657_rp, 0.0242925995028888_rp, 0.0255153607860429_rp, -0.382910379812211_rp, 0.214739411520287_rp, 0.0_rp, 0.00664126972595546_rp]

call spin_mat_to_eigen(mat6, smap, orb_eval, orb_evec, n0, spin_evec, err)
write (1, '(a, 3f14.10)') '"n0 noRF" ABS 1E-9', n0
do i = 1, 3
  write (1, '(a, 6f14.10)') '"spin_evec' // int_str(i) // ' RE noRF" ABS 1E-9', real(spin_evec(:,i), rp)
  write (1, '(a, 6f14.10)') '"spin_evec' // int_str(i) // ' IM noRF" ABS 1E-9', aimag(spin_evec(:,i))
enddo

mat6(1,:) = [-1.86188861871179_rp, 11.6415650949441_rp, -0.379398600528061_rp, -0.207432122016991_rp, 0.00304692758210734_rp, 4.04419952113488_rp]
mat6(2,:) = [-0.315028290447258_rp, 1.42969286047843_rp, -0.0480313327017454_rp, -0.0263103137420259_rp, 0.00203307819643807_rp, 0.889365901699835_rp]
mat6(3,:) = [-0.113800985432133_rp, -0.0627208591970316_rp, 0.90520668432316_rp, 4.29461375890543_rp, 0.0_rp, 0.000220836368260771_rp]
mat6(4,:) = [0.0031371439833448_rp, 0.00156512517121112_rp, -0.287935721862071_rp, -0.261367231459817_rp, 0.0_rp, 3.62949397953724e-05_rp]
mat6(5,:) = [0.383008304690509_rp, -4.56790453302225_rp, 0.142987043987865_rp, 0.0778167332865311_rp, 0.999599101473338_rp, -1.20504653880018_rp]
mat6(6,:) = [-0.000819009327739449_rp, -0.00462386992979728_rp, 0.000124247737021488_rp, 6.61961003550175e-05_rp, 0.00523961427182708_rp, 0.999599099135173_rp]

smap(0,:) = [0.627585943823847_rp, -0.400652158344974_rp, 0.548132807366073_rp, -0.0552283272173976_rp, -0.0136685745371593_rp, 0.000769662570518589_rp, 0.43832877203216_rp]
smap(1,:) = [-0.00577912096536385_rp, -0.022433622147659_rp, -0.023934389475195_rp, 0.643992369281101_rp, 2.15609605612843_rp, -7.07063211731558e-06_rp, -0.00228775603796164_rp]
smap(2,:) = [-0.778492991161124_rp, -0.323048089979131_rp, 0.441845517473186_rp, -0.0458892036514809_rp, -0.0290461360214925_rp, 0.000620467256904083_rp, 0.353328140106095_rp]
smap(3,:) = [-0.00715175423263742_rp, 0.0247012868571641_rp, 0.0230889442582705_rp, -0.371630255499699_rp, 0.220039460098232_rp, 5.71356858448481e-06_rp, 0.00541359028621752_rp]

call spin_mat_to_eigen(mat6, smap, orb_eval, orb_evec, n0, spin_evec, err)
write (1, '(a, 3f14.10)') '"n0 wRF" ABS 1E-9', n0
do i = 1, 3
  write (1, '(a, 6f14.10)') '"spin_evec' // int_str(i) // ' RE wRF" ABS 1E-9', real(spin_evec(:,i), rp)
  write (1, '(a, 6f14.10)') '"spin_evec' // int_str(i) // ' IM wRF" ABS 1E-9', aimag(spin_evec(:,i))
enddo

do i = 1, 3
  j = 2 * i - 1
  q = atan2(aimag(orb_eval(j)), real(orb_eval(j),rp)) / twopi
  call spin_quat_resonance_strengths(orb_evec(j,:), smap, xi_sum, xi_diff)
  write (1, '(a, i0, a, 3es14.6)') '"Spin-Res-', i, '" ABS 1E-9', q, xi_sum, xi_diff
enddo

!---------------------------------

call bmad_parser ('small_line.bmad', lat)

call spin_concat_linear_maps (err, map1, lat%branch(0), 0, 0)

!---------------------------------

call bmad_parser (lat_file, lat)

open (2, file = lat_file)
read (2, nml = param)
close (2)

call init_coord (orb0, lat%particle_start, lat%ele(0), downstream_end$)
call ptc_transfer_map_with_spin (lat%branch(0), t_ele%taylor, t_ele%spin_taylor, orb0, err_flag, 0, 1)

orb_start = orb0
orb_start%vec = orb_start%vec + dr

a_quat = track_taylor(orb_start%vec, t_ele%spin_taylor, t_ele%taylor%ref)
spin_a = quat_rotate (a_quat, orb0%spin)

bmad_com%spin_tracking_on = .true.
call track1 (orb_start, lat%ele(1), lat%param, orb_end)
spin_b = orb_end%spin

!

write (1, '(a, 3f14.9)') '"dPTC-Quad"   ABS 0   ', spin_a - orb0%spin
write (1, '(a, 3f14.9)') '"dBmad-Quad"  ABS 0   ', spin_b - orb0%spin


if (print_extra) then
  call type_taylors (t_ele%taylor)
  print *, '--------------------------------'
  call type_taylors (t_ele%spin_taylor)

  print '(a, 3f12.6)', 'Init:      ', orb0%spin
  print '(a, 3f12.6)', 'dPTC_Quad: ', spin_a - orb0%spin
  print '(a, 3f12.6)', 'dBmad-Quad:', spin_b - orb0%spin
endif

!

n_vec = [1.0_rp, 2.0_rp, 3.0_rp] / sqrt(14.0_rp)
orb2 = orb0
call rotate_spin (n_vec * 0.12_rp, orb2%spin)

write (1, '(a, 4es10.2)') '"dRot"   ABS 1e-10   ', orb2%spin
if (print_extra) then
  write (*, '(a, 4es10.2)') '"dRot" ABS 1e-10   ', orb2%spin
endif

!

ele => lat%ele(2)
bmad_com%spin_tracking_on = .true.

orb_start = orb0
ele%spin_tracking_method = taylor$
call track1 (orb_start, lat%ele(2), lat%param, orb_end)
write (1, '(a, 3f12.8)') '"Taylor-Taylor" ABS 1e-10  ', orb_end%spin

orb_start = orb0
ele%spin_tracking_method = symp_lie_ptc$
call track1 (orb_start, lat%ele(2), lat%param, orb_end)
write (1, '(a, 3f12.8)') '"PTC-Taylor" ABS 1e-10  ', orb_end%spin

!---------------------------------

call bmad_parser('small_ring.bmad', lat)
branch => lat%branch(0)
call closed_orbit_calc(lat, tao_branch%orbit, 6)

excite_zero = ''
call tao_spin_polarization_calc (branch, tao_branch, excite_zero, '')
call spin_concat_linear_maps(err, sm%map1, branch, 0, 0, orbit = tao_branch%orbit, excite_zero = excite_zero)

call spin_mat_to_eigen (sm%map1%orb_mat, sm%map1%spin_q, orb_eval, orb_evec, n0, spin_evec, err)
call spin_quat_resonance_strengths(orb_evec(j,:), sm%map1%spin_q, xi_sum, xi_diff)

write(1, '(a, f12.8)')   '"Polarization Limit ST" ABS 1e-8                   ', tao_branch%spin%pol_limit_st
write(1, '(a, f12.8)')   '"Polarization Limit DK" ABS 1e-8                   ', tao_branch%spin%pol_limit_dk
write(1, '(a, 3f12.8)')  '"Polarization Limits DK (a,b,c-modes)" ABS 1e-8    ', tao_branch%spin%pol_limit_dk_partial
write(1, '(a, 3f12.8)')  '"Polarization Limits DK (bc,ac,ab-modes)" ABS 1e-8 ', tao_branch%spin%pol_limit_dk_partial2

write(1, '(a, es16.8)')  '"Polarization Rate BKS" REL 1e-8       ', tao_branch%spin%pol_rate_bks
write(1, '(a, es16.8)')  '"Depolarization Rate" REL 1e-8         ', tao_branch%spin%depol_rate
write(1, '(a, 3es16.8)') '"Depolarization Rate Partial" REL 1e-8 ', tao_branch%spin%depol_rate_partial
write(1, '(a, 3es16.8)') '"Depolarization Rate Partial2" REL 1e-8', tao_branch%spin%depol_rate_partial2

write(1, '(a, es16.8)')  '"Integral g^3 * b_hat * n_0" REL 1e-8         ', tao_branch%spin%integral_bn
write(1, '(a, es16.8)')  '"Integral g^3 * b_hat * dn/ddelta" ABS 5E-17  ', tao_branch%spin%integral_bdn
write(1, '(a, es16.8)')  '"Integral g^3 (1 - 2(n * s_hat)/9)" REL 1e-8  ', tao_branch%spin%integral_1ns
write(1, '(a, es16.8)')  '"Integral g^3 * 11 (dn/ddelta)^2 / 9" REL 1E-8', tao_branch%spin%integral_dn2

do i = 1, 6
  vec = abs(orb_evec(i, :))
  select case (i)
  case (1, 2); vec(1:4) = vec([1,2,5,6])
  case (3, 4); vec(1:4) = vec([1,2,5,6])
  case (5, 6); vec(1:4) = vec([1,2,5,6])
  end select
  write (1, '(a, 6f12.6)')  '"orb-evec-' // int_str(i) // '" ABS 1e-8', vec(1:4)
  write (1, '(a, 3es16.8)') '"spin-re-evec-' // int_str(i) // '" REL 1e-8',  real(spin_evec(i, 1:3:2), 8)
  write (1, '(a, 3es16.8)') '"spin-im-evec-' // int_str(i) // '" REL 1e-8',  aimag(spin_evec(i, 1:3:2))
enddo

do i = 1, 3
  j = 2 * i - 1
  call spin_quat_resonance_strengths(orb_evec(j,:), sm%map1%spin_q, xi_sum, xi_diff)
  write (1, '(a, f13.7, 2(f17.7, es13.5))') '"Res-Strength-' // int_str(i) // '" REL 1e-8', xi_sum, xi_diff
enddo

!----------------------------

close (1)

end program
