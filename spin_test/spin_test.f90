program spin_test

use bmad
use ptc_spin_mod
use spin_mod
use taylor_mod

implicit none

type (lat_struct) lat
type (ele_struct) ele
type (coord_struct) orb0, orb_start, orb_end, orb1, orb2

real(rp) spin_a(3), spin_b(3), spin0(3), dr(6), a_quat(0:3), n_vec(3)

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

call init_coord (orb0, lat%beam_start, lat%ele(0), downstream_end$)
call transfer_map_calc_with_spin (lat, ele%taylor, ele%spin_taylor, orb0, err_flag, 0, 1)

orb_start = orb0
orb_start%vec = orb_start%vec + dr

spin0 = spinor_to_vec(orb0%spin)
spin_a = matmul (spin_taylor_to_mat(orb_start%vec, ele%taylor%ref, ele%spin_taylor), spin0)

bmad_com%spin_tracking_on = .true.
call track1 (orb_start, lat%ele(1), lat%param, orb_end)
spin_b = spinor_to_vec(orb_end%spin)

!

write (1, '(a, 3f14.9)') '"dPTC"   ABS 0   ', spin_a - spin0
write (1, '(a, 3f14.9)') '"dBmad"  ABS 0   ', spin_b - spin0


if (print_extra) then
  call type_taylors (ele%taylor)
  print *, '--------------------------------'
  call type_spin_taylors (ele%spin_taylor)

  print '(a, 3f12.6)', 'Init:   ', spin0
  print '(a, 3f12.6)', 'dPTC:   ', spin_a - spin0
  print '(a, 3f12.6)', 'dBmad:  ', spin_b - spin0
endif

!

orb1 = orb0
n_vec = [1.0_rp, 2.0_rp, 3.0_rp] / sqrt(14.0_rp)
a_quat = calc_rotation_quaternion (n_vec, 0.12_rp)
call quaternion_track (a_quat, orb1%spin)

orb2 = orb0
call rotate_spinor (n_vec * 0.12_rp, orb2%spin)

write (1, '(a, 4es10.2)') '"dRot"   ABS 1e-10   ', orb2%spin - orb1%spin
if (print_extra) then
  write (*, '(a, 4es10.2)') '"dRot" ABS 1e-10   ', orb2%spin - orb1%spin
endif

close (1)

end program
