!+
! Program reverse_test
!
! The basic idea is to: 
!   0) Start with some given particle coordinates.
!   1) Track the particle through an element, 
!   2) Reverse the particle's momentum.
!   3) Track the particle in reverse through the element.
!   4) Reverse the particle's momentum.
! Given the right conditions, the particle's final coordinates should equal the particle's 
! initial coordinates. This gives a check as to whether the reverse tracking is correct.
!
! The needed conditions for guaranteeing final = initial are:
!   A) For static magnectic fields: The particle charge is also reversed before reverse tracking.
!   B) For static magnectic fields: The particle charge is the same in reverse tracking.
!   C) For rfcavities with  longitudinal mirror symmetry: The particle charge is the same.
!      Mirror symmetry is valid for standaing wave cavities.
!-

program reverse_test

use bmad
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
character(40) :: lat_file  = 'reverse.bmad'
type (ele_struct), pointer :: ele_f, ele_r
type (coord_struct) orb_0f, orb_1f, orb_0r, orb_1r

real(rp) mat_f(6,6), m(6,6), vec1(6)
real(rp) dz, dpc, dct
logical :: err_flag
logical :: verbosity = .false.
integer nargs

nargs = cesr_iargc()
if (nargs == 1)then
   call cesr_getarg(1, lat_file)
   print *, 'Using ', trim(lat_file)
   verbosity = .true.
elseif (nargs > 1) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

! Init

open (1, file = 'output.now')

call bmad_parser (lat_file, lat)
!call write_bmad_lattice_file ('lat.bmad', lat)
!call bmad_parser ('lat.bmad', lat)

lat%branch(1)%ele(0)%value(e_tot$) = lat%branch(0)%ele(1)%value(e_tot$)
lat%branch(1)%ele(0)%value(p0c$) = lat%branch(0)%ele(1)%value(p0c$)

call lat_compute_ref_energy_and_time (lat, err_flag)

orb_0f = lat%beam_start
ele_f => lat%branch(0)%ele(1)
call init_coord (orb_0f, orb_0f%vec, ele_f, .false.)

call track1 (orb_0f, ele_f, ele_f%branch%param, orb_1f)
call make_mat6 (ele_f, ele_f%branch%param, orb_0f)

if (verbosity == .true.) then
   print '(a, 6es10.2, 5x, es10.2)', '0: ', orb_0f%vec, orb_0f%t
   print '(a, 6es10.2, 5x, es10.2)', '1: ', orb_1f%vec, orb_1f%t
end if

orb_0r        = orb_1f
orb_0r%vec(2) = -orb_1f%vec(2)
orb_0r%vec(4) = -orb_1f%vec(4)

orb_0r%vec(5) = orb_0f%vec(5)
orb_0r%vec(6) = orb_0f%vec(6)
orb_0r%t      = orb_0f%t
orb_0r%p0c    = orb_0f%p0c
orb_0r%beta   = orb_0f%beta

ele_r => lat%branch(1)%ele(1)
if (ele_r%key == rfcavity$ .or. ele_r%key == elseparator$) then
  ele_r%branch%param%rel_tracking_charge = ele_f%branch%param%rel_tracking_charge
endif

call track1 (orb_0r, ele_r, ele_r%branch%param, orb_1r)
call make_mat6 (ele_r, ele_r%branch%param, orb_0r)

if (verbosity == .true.) then
   print '(a, 6es10.2, 5x, es10.2)', '1: ', orb_0r%vec, orb_0r%t
   print '(a, 6es10.2, 5x, es10.2)', '2: ', orb_1r%vec, orb_1r%t
end if

!

mat_f(1, 1:6) =  ele_f%mat6(1, 1:6)
mat_f(2, 1:6) = -ele_f%mat6(2, 1:6)
mat_f(3, 1:6) =  ele_f%mat6(3, 1:6)
mat_f(4, 1:6) = -ele_f%mat6(4, 1:6)
mat_f(5, 1:6) = [0, 0, 0, 0, 1, 0]
mat_f(6, 1:6) = [0, 0, 0, 0, 0, 1]

call mat_inverse (mat_f, mat_f)

m(1, 1:6) = [1,  0, 0,  0, 0, 0]
m(2, 1:6) = [0, -1, 0,  0, 0, 0]
m(3, 1:6) = [0,  0, 1,  0, 0, 0]
m(4, 1:6) = [0,  0, 0, -1, 0, 0]
m(5, 1:6) = ele_f%mat6(5, 1:6)
m(6, 1:6) = ele_f%mat6(6, 1:6)
mat_f = matmul(m, mat_f)

!

dz  = (orb_1r%vec(5) - orb_0r%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5)) 
dpc = (orb_1r%vec(6) - orb_0r%vec(6)) - (orb_1f%vec(6) - orb_0f%vec(6)) 
dct = ((orb_1r%t - orb_0r%t) - (orb_1f%t - orb_0f%t)) * c_light

orb_1r%vec(2) = -orb_1r%vec(2)
orb_1r%vec(4) = -orb_1r%vec(4)

write (1, '(a, es11.3)') '"dorb(1)" ABS 1d-14 ', orb_1r%vec(1) - orb_0f%vec(1)
write (1, '(a, es11.3)') '"dorb(2)" ABS 1d-14 ', orb_1r%vec(2) - orb_0f%vec(2)
write (1, '(a, es11.3)') '"dorb(3)" ABS 1d-14 ', orb_1r%vec(3) - orb_0f%vec(3)
write (1, '(a, es11.3)') '"dorb(4)" ABS 1d-14 ', orb_1r%vec(4) - orb_0f%vec(4)
write (1, '(a, es11.3)') '"dorb(5)" ABS 1d-14 ', dz
write (1, '(a, es11.3)') '"dorb(6)" ABS 1d-14 ', dpc
write (1, '(a, es11.3)') '"c*dt"    ABS 1d-14 ', dct 

if (verbosity == .true.) then
   write (1, *)
   write (1, '(a, 6es11.3)') '"mat_f(1,:)" ABS 1d-14 ', ele_f%mat6(1,:)
   write (1, '(a, 6es11.3)') '"mat_f(2,:)" ABS 1d-14 ', ele_f%mat6(2,:)
   write (1, '(a, 6es11.3)') '"mat_f(3,:)" ABS 1d-14 ', ele_f%mat6(3,:)
   write (1, '(a, 6es11.3)') '"mat_f(4,:)" ABS 1d-14 ', ele_f%mat6(4,:)
   write (1, '(a, 6es11.3)') '"mat_f(5,:)" ABS 1d-14 ', ele_f%mat6(5,:)
   write (1, '(a, 6es11.3)') '"mat_f(6,:)" ABS 1d-14 ', ele_f%mat6(6,:)

   write (1, *)
   write (1, '(a, 6es11.3)') '"mat_r(1,:)" ABS 1d-14 ', ele_r%mat6(1,:)
   write (1, '(a, 6es11.3)') '"mat_r(2,:)" ABS 1d-14 ', ele_r%mat6(2,:)
   write (1, '(a, 6es11.3)') '"mat_r(3,:)" ABS 1d-14 ', ele_r%mat6(3,:)
   write (1, '(a, 6es11.3)') '"mat_r(4,:)" ABS 1d-14 ', ele_r%mat6(4,:)
   write (1, '(a, 6es11.3)') '"mat_r(5,:)" ABS 1d-14 ', ele_r%mat6(5,:)
   write (1, '(a, 6es11.3)') '"mat_r(6,:)" ABS 1d-14 ', ele_r%mat6(6,:)
end if

write (1, *)
write (1, '(a, 6es11.3)') '"dmat(1,:)" ABS 1d-14 ', ele_r%mat6(1,:) - mat_f(1,:)
write (1, '(a, 6es11.3)') '"dmat(2,:)" ABS 1d-14 ', ele_r%mat6(2,:) - mat_f(2,:)
write (1, '(a, 6es11.3)') '"dmat(3,:)" ABS 1d-14 ', ele_r%mat6(3,:) - mat_f(3,:)
write (1, '(a, 6es11.3)') '"dmat(4,:)" ABS 1d-14 ', ele_r%mat6(4,:) - mat_f(4,:)
write (1, '(a, 6es11.3)') '"dmat(5,:)" ABS 1d-14 ', ele_r%mat6(5,:) - mat_f(5,:)
write (1, '(a, 6es11.3)') '"dmat(6,:)" ABS 1d-14 ', ele_r%mat6(6,:) - mat_f(6,:)

if (verbosity == .true.) then
   write (1, *)
   write (1, '(a,  3es11.3)') '"max(dvec, dmat)" ', maxval(abs([orb_1r%vec(1:4)-orb_0f%vec(1:4), dz, dpc, dct])), &
                                     maxval(abs(ele_r%mat6 - mat_f))
end if

! And close

close (1)

end program
