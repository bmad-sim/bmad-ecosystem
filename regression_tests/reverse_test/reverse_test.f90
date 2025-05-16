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

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (ele_struct) ele_default
type (branch_struct), pointer :: branch

real(rp) max_diff_vec_r, max_diff_vec_b, max_diff_mat, max_diff_spin
integer, parameter :: n_methods = ubound(tracking_method_name, 1)
integer nargs, ie, ib, im
logical :: verbosity = .false., backwards = .false.

character(40) :: lat_file  = 'reverse.bmad', tracking_method = ''
character(100) :: str

namelist / params / tracking_method, backwards, bmad_com


!

nargs = command_argument_count()

if (nargs > 2) then
  print *, 'Max two command line args permitted.'
  call err_exit
endif

if (nargs > 0)then
 call get_command_argument(1, lat_file)
 print *, 'Using ', trim(lat_file)
 verbosity = .true.
endif

bmad_com%spin_tracking_on = .true.

if (nargs == 2) then
  open (1, file = 'reverse_test.init', status = 'old')
  read (1, nml = params)
  close (1)
endif

! Init

open (1, file = 'output.now')

call bmad_parser (lat_file, lat)

!

do ib = 0, ubound(lat%branch, 1)
  branch => lat%branch(ib)

  do ie = 1, lat%n_ele_max - 1 ! Do not test end marker
    ele => branch%ele(ie)

    ele_default%key = ele%key
    call set_ele_defaults(ele_default)

    max_diff_vec_r = 0
    max_diff_vec_b = 0
    max_diff_mat = 0
    max_diff_spin = 0

    do im = 1, n_methods
      if (im == fixed_step_runge_kutta$ .or. im == fixed_step_time_runge_kutta$) cycle
      if (.not. valid_tracking_method(ele, branch%param%particle, im)) cycle
      if (im == mad$) cycle
      if (im == custom$) cycle
      if (im == taylor$ .and. ele%key == rfcavity$) cycle
      if (im == taylor$ .and. ele%key == lcavity$) cycle
      if (im == symp_lie_ptc$ .and. ele%key == patch$) cycle
      if (tracking_method /= '' .and. upcase(tracking_method_name(im)) /= upcase(tracking_method)) cycle

      ele%tracking_method = im

      select case (ele%tracking_method)
      case (symp_lie_ptc$, symp_lie_bmad$, bmad_standard$, taylor$)
        ele%mat6_calc_method = ele%tracking_method
      case (linear$)
        ele%mat6_calc_method = ele_default%mat6_calc_method
      case default
        ele%mat6_calc_method = tracking$
      end select

      if (ele%tracking_method == symp_lie_ptc$) then
        ele%spin_tracking_method = symp_lie_ptc$
      else
        ele%spin_tracking_method = tracking$
      endif

      write (1, *)
      call test_this (ele)

    enddo

    write (1, *)
    write (1, '(a)') trim(line_out('"' // trim(ele%name) // '@', 'Largest_Max_Diff', &
                                              [max_diff_vec_r, max_diff_mat, max_diff_spin]))
  enddo
enddo

! And close


close (1)

!--------------------------------------------------------------------------
contains

subroutine test_this (ele)

type (ele_struct) ele, ele_r, ele_b
type (coord_struct) orb_0f, orb_1f, orb_0r_orient, orb_1r_orient, dorb_r_orient, orb_0b_track, orb_1b_track, dorb_b_track
type (coord_struct) orb_1r_orient_sav, orb_1b_track_sav

real(rp) dmat_r(6,6), dmat_b(6,6), m(6,6), vec1(6), dspin_r_orient(3), dspin_b_track(3)
real(rp) diff_vec_r, diff_vec_b, diff_mat, diff_spin, dvec0_r(6), diff_vec0_r

integer n
logical :: err_flag

!-------------------------------------------------------------------
! Forward tracking

call init_coord (orb_0f, lat%particle_start, ele, upstream_end$)

call make_mat6(ele, ele%branch%param, orb_0f)
call track1 (orb_0f, ele, ele%branch%param, orb_1f)

str = trim(ele%name) // '@' // tracking_method_name(ele%tracking_method)

if (verbosity) then
  print *, '!--------------------------------'
  print *, str
  print '(a, 6es12.4, 5x, es12.4)', '0: ', orb_0f%vec, orb_0f%t
  print '(a, 6es12.4, 5x, es12.4)', '1: ', orb_1f%vec, orb_1f%t
  print '(a, 3f13.7, 4x, 3f13.7)', 'Spin:', orb_0f%spin, orb_1f%spin - orb_0f%spin
end if

!-------------------------------------------------------------------
! Tracking with element reversed orientation but in the forward (orb%direction = 1) direction.

orb_0r_orient         = orb_1f
orb_0r_orient%vec(2)  = -orb_1f%vec(2)
orb_0r_orient%vec(4)  = -orb_1f%vec(4)  
orb_0r_orient%vec(5)  = -orb_1f%vec(5)  
orb_0r_orient%species = antiparticle(orb_0r_orient%species)
orb_0r_orient%t       = -orb_1f%t

ele_r = ele
ele_r%value(ref_time_start$) = -ele%ref_time  ! So z, t, and ref_t are consistant
ele_r%ref_time = ele_r%value(ref_time_start$) + ele_r%value(delta_ref_time$)
ele_r%orientation = -1

if (ele_r%key == elseparator$) then
  ele_r%value(hkick$) = -ele%value(hkick$)
  ele_r%value(vkick$) = -ele%value(vkick$)
elseif (ele_r%key == rfcavity$) then
  ele_r%value(phi0$)     = -ele%value(phi0$)
elseif (ele_r%key == lcavity$) then
  ele_r%value(phi0$)     = -ele%value(phi0$)
  ele_r%value(phi0_err$) = -ele%value(phi0_err$)
elseif (ele_r%key == patch$) then
   ele_r%value(upstream_coord_dir$) = -1
   ele_r%value(downstream_coord_dir$) = -1
endif

if (associated(ele_r%a_pole_elec)) then
  ele_r%a_pole_elec = -ele_r%a_pole_elec
  ele_r%b_pole_elec = -ele_r%b_pole_elec
endif

call kill_taylor (ele_r%taylor)
call kill_taylor (ele_r%spin_taylor)
call make_mat6(ele_r, ele%branch%param, orb_0r_orient)
call track1(orb_0r_orient, ele_r, ele_r%branch%param, orb_1r_orient)

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

dvec0_r = matmul(ele_r%mat6, [ele%vec0(1), -ele%vec0(2), ele%vec0(3), -ele%vec0(4), -ele%vec0(5), ele%vec0(6)])
dvec0_r = ele_r%vec0 + dvec0_r

!-------------------------------------------------------------------
! Tracking backwards (element is unreversed).

orb_0b_track   = orb_0r_orient
orb_0b_track%t = -orb_1f%t
orb_0b_track%direction = -1
orb_0b_track%vec(5) = -orb_0b_track%vec(5)

ele_b = ele
ele_b%ref_time    = -ele%ref_time
ele_b%value(ref_time_start$) = ele_b%ref_time - ele_b%value(delta_ref_time$)

if (ele_r%key == elseparator$) then
  ele_r%value(hkick$) = -ele%value(hkick$)
  ele_r%value(vkick$) = -ele%value(vkick$)
elseif (ele_r%key == rfcavity$) then
  ele_r%value(phi0$)     = -ele%value(phi0$)
elseif (ele_r%key == lcavity$) then
  ele_r%value(phi0$)     = -ele%value(phi0$)
  ele_r%value(phi0_err$) = -ele%value(phi0_err$)
endif

if (associated(ele_r%a_pole_elec)) then
  ele_r%a_pole_elec = -ele_r%a_pole_elec
  ele_r%b_pole_elec = -ele_r%b_pole_elec
endif

call track1(orb_0b_track, ele_b, ele%branch%param, orb_1b_track)

orb_1b_track_sav = orb_1b_track

orb_1b_track%vec(2)    = -orb_1b_track%vec(2)
orb_1b_track%vec(4)    = -orb_1b_track%vec(4)  

dorb_b_track%vec    = orb_1b_track%vec - orb_0f%vec
if (orb_1b_track%beta == orb_0b_track%beta) then
  dorb_b_track%vec(5) = (orb_1b_track%vec(5) - orb_0b_track%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5))
  dorb_b_track%vec(5) = dorb_b_track%vec(5) + 2 * orb_1b_track%beta * c_light * ele_b%value(delta_ref_time$)
else
  dorb_b_track%vec(5) = 0
endif
dorb_b_track%t      = (orb_1b_track%t - orb_0b_track%t) - (orb_1f%t - orb_0f%t)
dspin_b_track       = orb_1b_track%spin - orb_0f%spin

! Matrix

if (backwards) then
  call mat_inverse(ele_b%mat6, dmat_b)
  dmat_b(2,:) = -dmat_b(2,:)
  dmat_b(4,:) = -dmat_b(4,:)
  dmat_b(5,:) = -dmat_b(5,:)
  dmat_b(:,2) = -dmat_b(:,2)
  dmat_b(:,4) = -dmat_b(:,4)
  dmat_b(:,5) = -dmat_b(:,5)
  dmat_b = ele%mat6 - dmat_b
endif

!-------------------------------------------------------------------

if (verbosity) then
  print *
  print '(a, 6es12.4, 5x, es12.4)', '0r_orient:     ', orb_0r_orient%vec, orb_0r_orient%t
  print '(a, 6es12.4, 5x, es12.4)', '1r_orient:     ', orb_1r_orient_sav%vec, orb_1r_orient%t
  print '(a, 6es12.4, 5x, es12.4)', 'Dr_orient:     ', dorb_r_orient%vec(1:6), dorb_r_orient%t
  print '(a, 3es12.4)',             'dSpin_r_orient:', dspin_r_orient
  print *
  if (backwards) then
    print '(a, 6es12.4, 5x, es12.4)', '0b_track:     ', orb_0b_track%vec, orb_0b_track%t
    print '(a, 6es12.4, 5x, es12.4)', '1b_track:     ', orb_1b_track_sav%vec, orb_1b_track%t
    print '(a, 6es12.4, 5x, es12.4)', 'Db_track:     ', dorb_b_track%vec(1:6), dorb_b_track%t
    print '(a, 3es12.4)',             'dSpin_b_track:', dspin_b_track
    print *
  endif
  do n = 1, 6
    print '(6f13.7, 5x, f13.7)', ele%mat6(n,:), ele%vec0(n)
  enddo
  print *
  do n = 1, 6
    print '(6f13.7, 5x, f13.7)', ele_r%mat6(n,:) , ele_r%vec0(n)
  enddo
  print *
  do n = 1, 6
    print '(6f13.7, 5x, f13.7)', dmat_r(n,:), dvec0_r(n)
  enddo
  print *
end if

!

str = '"' // trim(ele%name) // '@' // trim(tracking_method_name(ele%tracking_method)) // ':'

write (1, '(a)') trim(line_out(str, 'dorb_r_orient"', [real(rp):: dorb_r_orient%vec, c_light * dorb_r_orient%t]))
if (backwards) write (1, '(a)') trim(line_out(str, 'dorb_b_track" ', [real(rp):: dorb_b_track%vec, c_light * dorb_b_track%t]))
write (1, '(a)') trim(line_out(str, 'xmat_r"', [maxval(abs(dmat_r))]))
if (backwards) write (1, '(a)') trim(line_out(str, 'xmat_b"', [maxval(abs(dmat_b))]))
write (1, '(a)') trim(line_out(str, 'dspin_r_orient"', dspin_r_orient))
if (backwards) write (1, '(a)') trim(line_out(str, 'dspin_b_track"', dspin_b_track))

diff_vec_r = maxval([real(rp):: abs(dorb_r_orient%vec), abs(dorb_r_orient%t)])
if (backwards) diff_vec_b = maxval([real(rp):: abs(dorb_b_track%vec), abs(dorb_b_track%t)])
diff_mat  = max(maxval(abs(dmat_r)), maxval(abs(dvec0_r)))
if (backwards) diff_spin = max(maxval(abs(dspin_r_orient)), maxval(abs(dspin_b_track)))
diff_spin = maxval(abs(dspin_r_orient))

max_diff_vec_r = max(max_diff_vec_r, diff_vec_r)
if (backwards) max_diff_vec_b = max(max_diff_vec_b, diff_vec_b)
max_diff_mat  = max(max_diff_mat, diff_mat)
max_diff_spin = max(max_diff_spin, diff_spin)

if (verbosity) then
if (backwards) print '(2a, t45, 4es10.2)', 'Max Diff: ', trim(str), diff_vec_r, diff_vec_b, diff_mat, diff_spin
  print '(2a, t45, 4es10.2)', 'Max Diff: ', trim(str), diff_vec_r, diff_mat, diff_spin
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

tol = 'ABS 1E-12'

select case (str_out)
end select

write (str_out, '(a, t52, a, t62, 7es12.4)') trim(str_out), tol, val

end function line_out

end program


