

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
!   B) For static electric fields: The particle charge is the same in reverse tracking.
!
! For rfcaities life is more complicated since the fields are time dependent and the arrow 
! of time cannot be reversed. In this case, longitudinal mirror symmetry must be used.
! The procedure for rfcavities is then:
!   0) Start with some given particle coordinates.
!   1) Track the particle through an element, 
!   2) Start with the same initial coords and same charge and track in reverse.
!      Note: Must shift the time and z to be at the correct phase of the accelerating backward wave.
!-

program reverse_test

use bmad
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
real(rp) max_diff_vec_r, max_diff_vec_b, max_diff_mat
integer nargs, i
logical :: verbosity = .false.

character(40) :: lat_file  = 'reverse.bmad'
character(100) :: str

!

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

max_diff_vec_r = 0
max_diff_vec_b = 0
max_diff_mat = 0

!

do i = 1, lat%n_ele_max - 1

  ele => lat%ele(i)
  call test_this (ele)

  if (valid_tracking_method (ele, lat%param%particle, runge_kutta$)) then
    ele%tracking_method = runge_kutta$
    call test_this (ele)
  endif

!  if (valid_tracking_method (ele, lat%param%particle, symp_lie_ptc$)) then
!    ele%tracking_method = symp_lie_ptc$
!    call test_this (ele)
!  endif

  if (valid_tracking_method (ele, lat%param%particle, symp_lie_bmad$)) then
    ele%tracking_method = symp_lie_bmad$
    call test_this (ele)
  endif

enddo

! And close

close (1)
print '(a, 3es12.3)', 'Largest Max Diff: ', max_diff_vec_r, max_diff_vec_b, max_diff_mat

!--------------------------------------------------------------------------
contains

subroutine test_this (ele)

type (ele_struct) ele, ele2
type (coord_struct) orb_0f, orb_1f, orb_0r, orb_1r, dorb_r, orb_0b, orb_1b, dorb_b
type (coord_struct) orb_1r_sav, orb_1b_sav

real(rp) dmat(6,6), m(6,6), vec1(6)
real(rp) diff_vec_r, diff_vec_b, diff_mat

integer n
logical :: err_flag

!
orb_0f = lat%beam_start
call init_coord (orb_0f, orb_0f%vec, ele, upstream_end$)

select case (ele%tracking_method)
case (runge_kutta$)
  ele%mat6_calc_method = tracking$
case (symp_lie_ptc$, symp_lie_bmad$)
  ele%mat6_calc_method = ele%tracking_method
end select

call track1 (orb_0f, ele, ele%branch%param, orb_1f)
call make_mat6(ele, ele%branch%param, orb_0f, orb_1f, .true.)

str = trim(ele%name) // '@' // tracking_method_name(ele%tracking_method)

if (verbosity) then
  print *, '!--------------------------------'
  print *, str
  print '(a, 6es12.4, 5x, es12.4)', '0: ', orb_0f%vec, orb_0f%t
  print '(a, 6es12.4, 5x, es12.4)', '1: ', orb_1f%vec, orb_1f%t
end if

! Reverse tracking

orb_0r           = orb_1f
orb_0r%vec(2)    = -orb_1f%vec(2)
orb_0r%vec(4)    = -orb_1f%vec(4)  

ele2 = ele
if (ele2%key == elseparator$) then
elseif (ele2%key == rfcavity$) then
  orb_0r%species   = flip_species_charge(orb_0r%species)
  orb_0r%vec(5) = orb_0f%vec(5)
  orb_0r%t      = orb_0f%t

elseif (ele2%key == patch$) then
   ele2%value(upstream_ele_dir$) = -1
   ele2%value(downstream_ele_dir$) = -1
else
  orb_0r%species   = flip_species_charge(orb_0r%species)
endif

ele2%orientation = -1

call track1(orb_0r, ele2, ele2%branch%param, orb_1r)
orb_1r_sav = orb_1r

orb_1r%vec(2)    = -orb_1r%vec(2)
orb_1r%vec(4)    = -orb_1r%vec(4)  

dorb_r%vec    = orb_1r%vec - orb_0f%vec
dorb_r%vec(5) = (orb_1r%vec(5) - orb_0r%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5))
dorb_r%t      = (orb_1r%t - orb_0r%t) - (orb_1f%t - orb_0f%t)

! Backwards tracking

orb_0b = orb_0r
orb_0b%direction = -1

call track1(orb_0b, ele, ele%branch%param, orb_1b)
orb_1b_sav = orb_1b

orb_1b%vec(2)    = -orb_1b%vec(2)
orb_1b%vec(4)    = -orb_1b%vec(4)  

dorb_b%vec    = orb_1b%vec - orb_0f%vec
dorb_b%vec(5) = (orb_1b%vec(5) - orb_0b%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5))
dorb_b%t      = (orb_1b%t - orb_0b%t) - (orb_1f%t - orb_0f%t)

! Matrix

call make_mat6(ele2, ele%branch%param, orb_0r, orb_1r, .true.)

call mat_inverse(ele2%mat6, dmat)
dmat(2,:) = -dmat(2,:)
dmat(4,:) = -dmat(4,:)
dmat(5,:) = -dmat(5,:)
dmat(:,2) = -dmat(:,2)
dmat(:,4) = -dmat(:,4)
dmat(:,5) = -dmat(:,5)
dmat = ele%mat6 - dmat

if (verbosity) then
  print *
  print '(a, 6es12.4, 5x, es12.4)', '1r: ', orb_0r%vec, orb_0r%t
  print '(a, 6es12.4, 5x, es12.4)', '2r: ', orb_1r_sav%vec, orb_1r%t
  print '(a, 6es12.4, 5x, es12.4)', 'Dr: ', dorb_r%vec(1:6), dorb_r%t
  print *
  print '(a, 6es12.4, 5x, es12.4)', '1b: ', orb_0b%vec, orb_0b%t
  print '(a, 6es12.4, 5x, es12.4)', '2b: ', orb_1b_sav%vec, orb_1b%t
  print '(a, 6es12.4, 5x, es12.4)', 'Db: ', dorb_b%vec(1:6), dorb_b%t
  print *
  do n = 1, 6
    print '(6f12.6)', ele%mat6(n,:)
  enddo
  print *
  do n = 1, 6
    print '(6f12.6)', ele2%mat6(n,:) 
  enddo
  print *
  do n = 1, 6
    print '(6f12.6)', dmat(n,:) 
  enddo
  print *
end if

!

str = '"' // trim(ele%name) // '@' // trim(tracking_method_name(ele%tracking_method)) // ':'

write (1, '(a)') trim(line_out(str, 'dorb_r"', maxval(abs(dorb_r%vec))))
write (1, '(a)') trim(line_out(str, 'c*dt_r"', dorb_r%t))
write (1, '(a)') trim(line_out(str, 'dorb_b"', maxval(abs(dorb_b%vec))))
write (1, '(a)') trim(line_out(str, 'c*dt_b"', dorb_b%t))
write (1, '(a)') trim(line_out(str, 'xmat"', maxval(abs(dmat))))

diff_vec_r = maxval([dorb_r%vec, dorb_r%t])
diff_vec_b = maxval([dorb_b%vec, dorb_b%t])
diff_mat = maxval(abs(dmat))

max_diff_vec_r = max(max_diff_vec_r, diff_vec_r)
max_diff_vec_b = max(max_diff_vec_b, diff_vec_b)
max_diff_mat = max(max_diff_mat, diff_mat)

if (verbosity) then
  print '(2a, t50, 3es10.2)', 'Max Diff: ', trim(str), diff_vec_r, diff_vec_b, diff_mat
endif

end subroutine test_this

!-------------------------------------------------------------------------
! contains

function line_out(str1, str2, val) result (str_out)

real(rp) val
character(*) str1, str2
character(100) str_out
character(16) tol

!

str_out = trim(str1) // str2

tol = 'REL 0.1'
if (abs(val) < 1e-14) tol = 'ABS 1e-14'
if (index(str_out, 'Runge') /= 0) tol = 'ABS 1e-9'

select case (str_out)
case ('"QUADRUPOLE5@Bmad_Standard:dorb_r"');    tol = 'ABS 1e-9'
case ('"QUADRUPOLE5@Symp_Lie_Bmad:dorb_r"');    tol = 'ABS 1e-9'
case ('"RFCAVITY1@Runge_Kutta:dorb_r"');        tol = 'ABS 3e-9'
case ('"RFCAVITY2@Bmad_Standard:dorb_r"');      tol = 'ABS 1e-6'
case ('"SBEND4@Bmad_Standard:dorb_r"');         tol = 'ABS 2e-8'
case ('"SBEND4@Runge_Kutta:dorb_r"');           tol = 'ABS 1e-7'
case ('"SBEND4@Runge_Kutta:xmat"');             tol = 'ABS 1e-9'
end select

write (str_out, '(a, t47, a, t60, es12.4)') trim(str_out), tol, val

end function line_out

end program


