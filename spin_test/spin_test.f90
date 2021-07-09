program spin_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct) t_ele
type (coord_struct) orb0, orb_start, orb_end, orb1, orb2

real(rp) spin_a(3), spin_b(3), spin0(3), dr(6), a_quat(0:3), n_vec(3)
real(rp) mat6(6,6)

integer nargs

character(40) :: lat_file = 'spin_test.bmad'
logical print_extra, err_flag

namelist / param / dr

!                  

global_com%exit_on_error = .false.

print_extra = .false.
nargs = cesr_iargc() 
if (nargs > 1) then  
  print *, 'Only one command line arg permitted.'
  call err_exit                                  

elseif (nargs > 0)then
  call cesr_getarg(1, lat_file)
  print *, 'Using ', trim(lat_file)
  print_extra = .true.             
endif                              

call bmad_parser (lat_file, lat)

open (1, file = lat_file)
read (1, nml = param)
close (1)

!

open (1, file = 'output.now')

!

!mat6(1,:) = [-1.795314806127,  11.169503607591,  -0.072040449858,  -0.309254514765,   0.139536527474,   2.500915635320]
!mat6(2,:) = [-0.307551607326,   1.376714853191,  -0.004285833502,  -0.049235846905,   0.022015858160,   0.639803917413]
!mat6(3,:) = [ 0.025950938219,  -0.270713920515,  -1.173200609674,  -2.134837975277,   0.003863406582,   0.161444030610]
!mat6(4,:) = [ 0.009574724184,  -0.014681768689,   0.142691836094,  -0.590953501728,  -0.000296313025,  -0.049128246447]
!mat6(5,:) = [ 0.176693701798,  -2.992991984786,   0.018015388299,  -0.118086642355,   0.570940266616, -14.232737739618]
!mat6(6,:) = [ 0.019747369685,  -0.140791393136,  -0.000671030912,  -0.004067544028,   0.044879709376,   0.573027752937]

! 

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

close (1)

end program
