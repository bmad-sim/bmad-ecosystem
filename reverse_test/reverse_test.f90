!+
! Program reverse_test
!
! The basic idea is to: 
!   0) Start with some given particle coordinates.
!   1) Track the particle through an element, 
!   2) Reverse the particle's momentum.
!   3) Track the particle in reverse through the element.
!   4) Reverse the particle's momentum.
! Given the right conditions, the particle's final coordinates and spin should equal the particle's 
! initial coordinates. This gives a check as to whether the reverse tracking is correct.
!
! The needed conditions for guaranteeing final = initial are:
!   A) For static magnectic fields: The particle charge is also reversed before reverse tracking.
!   B) For static electric fields: The particle charge is the same in reverse tracking.
!
! For rfcaities, life is more complicated since the fields are time dependent.
! In this case, longitudinal mirror symmetry must be used.
!-

program reverse_test

use bmad
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

real(rp) max_diff_vec_r, max_diff_vec_b, max_diff_mat, max_diff_spin
integer nargs, ie, ib, im
logical :: verbosity = .false.

character(40) :: lat_file  = 'reverse.bmad', tracking_method = ''
character(100) :: str

!

nargs = cesr_iargc()

if (nargs > 2) then
  print *, 'Only one command line arg permitted.'
  call err_exit
endif

if (nargs > 0)then
 call cesr_getarg(1, lat_file)
 print *, 'Using ', trim(lat_file)
 verbosity = .true.
endif

if (nargs == 2) then
  call cesr_getarg(2, tracking_method)
endif

bmad_com%spin_tracking_on = .true.

! Init

open (1, file = 'output.now')

call bmad_parser (lat_file, lat)

max_diff_vec_r = 0
max_diff_vec_b = 0
max_diff_mat = 0
max_diff_spin = 0

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do ie = 1, lat%n_ele_max - 1 ! Do not test end marker
    ele => branch%ele(ie)

    do im = 1, n_methods$
      if (.not. valid_tracking_method(ele, branch%param%particle, im)) cycle
      if (im == mad$  .or. im == symp_map$ .or. im == custom$) cycle
      if ((ele%key == rfcavity$ .or. ele%key == lcavity$) .and. im == taylor$) cycle
      if (tracking_method /= '' .and. upcase(tracking_method_name(im)) /= upcase(tracking_method)) cycle

      ele%tracking_method = im

      if (ele%tracking_method == symp_lie_ptc$) then
        ele%spin_tracking_method = symp_lie_ptc$
      else
        ele%spin_tracking_method = tracking$
      endif

      write (1, *)
      call test_this (ele)

    enddo
  enddo
enddo

! And close

close (1)
if (verbosity) then
  print *
  print *, '!-------------------------------------------------'
  print '(a, 4es12.3)', 'Largest Max Diff: ', max_diff_vec_r, max_diff_vec_b, max_diff_mat, max_diff_spin
endif

!--------------------------------------------------------------------------
contains

subroutine test_this (ele)

type (ele_struct) ele, ele_r, ele_b
type (coord_struct) orb_0f, orb_1f, orb_0r_orient, orb_1r_orient, dorb_r_orient, orb_0b_track, orb_1b_track, dorb_b_track
type (coord_struct) orb_1r_orient_sav, orb_1b_track_sav

real(rp) dmat_r(6,6), dmat_b(6,6), m(6,6), vec1(6), dspin_r_orient(3), dspin_b_track(3)
real(rp) diff_vec_r, diff_vec_b, diff_mat, diff_spin

integer n
logical :: err_flag

!-------------------------------------------------------------------
! Forward tracking

call init_coord (orb_0f, lat%beam_start, ele, upstream_end$)

ele%mat6_calc_method = tracking$

select case (ele%tracking_method)
case (symp_lie_ptc$, symp_lie_bmad$, bmad_standard$)
  ele%mat6_calc_method = ele%tracking_method
end select

call track1 (orb_0f, ele, ele%branch%param, orb_1f)
call make_mat6(ele, ele%branch%param, orb_0f)

str = trim(ele%name) // '@' // tracking_method_name(ele%tracking_method)

if (verbosity) then
  print *, '!--------------------------------'
  print *, str
  print '(a, 6es12.4, 5x, es12.4)', '0: ', orb_0f%vec, orb_0f%t
  print '(a, 6es12.4, 5x, es12.4)', '1: ', orb_1f%vec, orb_1f%t
  print '(a, 3f12.6, 4x, 3f12.6)', 'Spin:', orb_0f%spin, orb_1f%spin - orb_0f%spin
end if

!-------------------------------------------------------------------
! Tracking with element reversed orientation (in the forward direction)

orb_0r_orient         = orb_1f
orb_0r_orient%vec(2)  = -orb_1f%vec(2)
orb_0r_orient%vec(4)  = -orb_1f%vec(4)  
orb_0r_orient%vec(5)  = -orb_1f%vec(5)  
orb_0r_orient%species = antiparticle(orb_0r_orient%species)
orb_0r_orient%t       = -orb_1f%t

ele_r = ele
if (ele_r%key == elseparator$) then
  ele_r%value(hkick$) = -ele%value(hkick$)
  ele_r%value(vkick$) = -ele%value(vkick$)
elseif (ele_r%key == rfcavity$) then
  ele_r%value(phi0$)     = -ele%value(phi0$)
elseif (ele_r%key == lcavity$) then
  ele_r%value(phi0$)     = -ele%value(phi0$)
  ele_r%value(phi0_err$) = -ele%value(phi0_err$)
elseif (ele_r%key == patch$) then
   ele_r%value(upstream_ele_dir$) = -1
   ele_r%value(downstream_ele_dir$) = -1
endif

if (associated(ele_r%a_pole_elec)) then
  ele_r%a_pole_elec = -ele_r%a_pole_elec
  ele_r%b_pole_elec = -ele_r%b_pole_elec
endif

ele_r%orientation = -1

call track1(orb_0r_orient, ele_r, ele_r%branch%param, orb_1r_orient)
call make_mat6(ele_r, ele%branch%param, orb_0r_orient)

orb_1r_orient_sav = orb_1r_orient

orb_1r_orient%vec(2)    = -orb_1r_orient%vec(2)
orb_1r_orient%vec(4)    = -orb_1r_orient%vec(4)  

dorb_r_orient%vec    = orb_1r_orient%vec - orb_0f%vec
dorb_r_orient%vec(5) = (orb_1r_orient%vec(5) - orb_0r_orient%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5))
dorb_r_orient%t      = (orb_1r_orient%t - orb_0r_orient%t) - (orb_1f%t - orb_0f%t)
dspin_r_orient       = orb_1r_orient%spin - orb_0f%spin

! Matrix

call mat_inverse(ele_r%mat6, dmat_r)
dmat_r(2,:) = -dmat_r(2,:)
dmat_r(4,:) = -dmat_r(4,:)
dmat_r(5,:) = -dmat_r(5,:)
dmat_r(:,2) = -dmat_r(:,2)
dmat_r(:,4) = -dmat_r(:,4)
dmat_r(:,5) = -dmat_r(:,5)
dmat_r = ele%mat6 - dmat_r

!-------------------------------------------------------------------
! Tracking backwards (element is unreversed).

orb_0b_track = orb_0r_orient
orb_0b_track%direction = -1

ele_b = ele_r
ele_b%orientation = 1
if (ele_b%key == patch$) then
   ele_b%value(upstream_ele_dir$) = +1
   ele_b%value(downstream_ele_dir$) = +1
endif

call track1(orb_0b_track, ele_b, ele%branch%param, orb_1b_track)
call make_mat6(ele_b, ele%branch%param, orb_0b_track)

orb_1b_track_sav = orb_1b_track

orb_1b_track%vec(2)    = -orb_1b_track%vec(2)
orb_1b_track%vec(4)    = -orb_1b_track%vec(4)  

dorb_b_track%vec    = orb_1b_track%vec - orb_0f%vec
dorb_b_track%vec(5) = (orb_1b_track%vec(5) - orb_0b_track%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5))
dorb_b_track%t      = (orb_1b_track%t - orb_0b_track%t) - (orb_1f%t - orb_0f%t)
dspin_b_track       = orb_1b_track%spin - orb_0f%spin

! Matrix

call mat_inverse(ele_b%mat6, dmat_b)
dmat_b(2,:) = -dmat_b(2,:)
dmat_b(4,:) = -dmat_b(4,:)
dmat_b(5,:) = -dmat_b(5,:)
dmat_b(:,2) = -dmat_b(:,2)
dmat_b(:,4) = -dmat_b(:,4)
dmat_b(:,5) = -dmat_b(:,5)
dmat_b = ele%mat6 - dmat_b

if (verbosity) then
  print *
  print '(a, 6es12.4, 5x, es12.4)', '0r_orient:     ', orb_0r_orient%vec, orb_0r_orient%t
  print '(a, 6es12.4, 5x, es12.4)', '1r_orient:     ', orb_1r_orient_sav%vec, orb_1r_orient%t
  print '(a, 6es12.4, 5x, es12.4)', 'Dr_orient:     ', dorb_r_orient%vec(1:6), dorb_r_orient%t
  print '(a, 3es12.4)',             'dSpin_r_orient:', dspin_r_orient
  print *
  print '(a, 6es12.4, 5x, es12.4)', '0b_track:     ', orb_0b_track%vec, orb_0b_track%t
  print '(a, 6es12.4, 5x, es12.4)', '1b_track:     ', orb_1b_track_sav%vec, orb_1b_track%t
  print '(a, 6es12.4, 5x, es12.4)', 'Db_track:     ', dorb_b_track%vec(1:6), dorb_b_track%t
  print '(a, 3es12.4)',             'dSpin_b_track:', dspin_b_track
  print *
  do n = 1, 6
    print '(6f12.6)', ele%mat6(n,:)
  enddo
  print *
  do n = 1, 6
    print '(6f12.6, 5x, 6f12.6)', ele_r%mat6(n,:), ele_b%mat6(n,:) 
  enddo
  print *
  do n = 1, 6
    print '(6f12.6, 5x, 6f12.6)', dmat_r(n,:), dmat_b(n,:)
  enddo
  print *
end if

!

str = '"' // trim(ele%name) // '@' // trim(tracking_method_name(ele%tracking_method)) // ':'

write (1, '(a)') trim(line_out(str, 'dorb_r_orient"', [dorb_r_orient%vec, c_light * dorb_r_orient%t]))
write (1, '(a)') trim(line_out(str, 'dorb_b_track" ', [dorb_b_track%vec, c_light * dorb_b_track%t]))
write (1, '(a)') trim(line_out(str, 'xmat_r"', [maxval(abs(dmat_r))]))
write (1, '(a)') trim(line_out(str, 'xmat_b"', [maxval(abs(dmat_b))]))
write (1, '(a)') trim(line_out(str, 'dspin_r_orient"', dspin_r_orient))
write (1, '(a)') trim(line_out(str, 'dspin_b_track"', dspin_b_track))

diff_vec_r = maxval([abs(dorb_r_orient%vec), abs(dorb_r_orient%t)])
diff_vec_b = maxval([abs(dorb_b_track%vec), abs(dorb_b_track%t)])
diff_mat  = max(maxval(abs(dmat_r)), maxval(abs(dmat_b)))
diff_spin = max(maxval(abs(dspin_r_orient)), maxval(abs(dspin_b_track)))

max_diff_vec_r = max(max_diff_vec_r, diff_vec_r)
max_diff_vec_b = max(max_diff_vec_b, diff_vec_b)
max_diff_mat  = max(max_diff_mat, diff_mat)
max_diff_spin = max(max_diff_spin, diff_spin)

if (verbosity) then
  print '(2a, t40, 4es10.2)', 'Max Diff: ', trim(str), diff_vec_r, diff_vec_b, diff_mat, diff_spin
endif

end subroutine test_this

!-------------------------------------------------------------------------
! contains

function line_out(str1, str2, val) result (str_out)

real(rp) val(:)
character(*) str1, str2
character(200) str_out
character(16) tol

!

str_out = trim(str1) // str2

tol = 'REL 1E-5'

select case (str_out)
case ('"QUADRUPOLE5@Bmad_Standard:dorb_r_orient"');    tol = 'ABS 1e-9'
case ('"QUADRUPOLE5@Symp_Lie_Bmad:dorb_r_orient"');    tol = 'ABS 1e-9'
case ('"RFCAVITY1@Runge_Kutta:dorb_r_orient"');        tol = 'ABS 3e-9'
case ('"RFCAVITY2@Bmad_Standard:dorb_r_orient"');      tol = 'ABS 1e-6'
case ('"SBEND4@Bmad_Standard:dorb_b_track"');          tol = 'ABS 2e-8'
case ('"SBEND4@Runge_Kutta:dorb_b_track"');            tol = 'ABS 1e-7'
case ('"SBEND4@Bmad_Standard:dorb_r_orient"');         tol = 'ABS 2e-8'
case ('"SBEND4@Runge_Kutta:dorb_r_orient"');           tol = 'ABS 1e-7'
case ('"SBEND4@Runge_Kutta:xmat"');                    tol = 'ABS 2e-5'
end select

write (str_out, '(a, t47, a, t60, 7es12.4)') trim(str_out), tol, val

end function line_out

end program


