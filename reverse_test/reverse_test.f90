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
! The particle's final coords tracking forward and reversed should be the same modulo
! a constant time and z offset
!-

program reverse_test

use bmad
use write_lat_file_mod

implicit none

type (lat_struct), target :: lat
character(40) :: lat_file  = 'reverse.bmad'
type (ele_struct), pointer :: ele_f, ele_r
type (coord_struct) orb_0f, orb_1f, orb_0r, orb_1r, dorb

real(rp) mat_f(6,6), m(6,6), vec1(6)
real(rp) beta_ref, dt_ref
logical :: err_flag
logical :: verbosity = .false.
integer nargs, i, length
character(20) :: str

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

DO i = 1, lat%n_ele_max - 1

  orb_0f = lat%beam_start
  ele_f => lat%branch(0)%ele(i)
  call init_coord (orb_0f, orb_0f%vec, ele_f, upstream_end$)

  call track1 (orb_0f, ele_f, ele_f%branch%param, orb_1f)
  call make_mat6 (ele_f, ele_f%branch%param, orb_0f)

  str = trim(ele_f%name)

  if (verbosity) then
    print *, str
    print '(a, 6es12.4, 5x, es12.4)', '0: ', orb_0f%vec, orb_0f%t
    print '(a, 6es12.4, 5x, es12.4)', '1: ', orb_1f%vec, orb_1f%t
  end if

  orb_0r        = orb_1f
  orb_0r%vec(2) = -orb_1f%vec(2)
  orb_0r%vec(4) = -orb_1f%vec(4)

  ele_r => lat%branch(1)%ele(i)

  if (ele_f%key == elseparator$) then
    ele_r%branch%param%rel_tracking_charge = ele_f%branch%param%rel_tracking_charge
  elseif (ele_f%key == rfcavity$) then
    ele_r%branch%param%rel_tracking_charge = ele_f%branch%param%rel_tracking_charge
    ele_r%value(x_pitch$) = -ele_f%value(x_pitch$)
    ele_r%value(y_pitch$) = -ele_f%value(y_pitch$)
    orb_0r = orb_0f
    beta_ref = ele_f%value(p0c$) / ele_f%value(e_tot$)
    dt_ref = ele_f%value(l$) / (c_light * beta_ref)
    orb_0r%vec(5) = orb_0r%vec(5) - dt_ref * orb_0r%beta
    orb_0r%t = orb_0r%t + dt_ref
  else
    ele_r%branch%param%rel_tracking_charge = -ele_f%branch%param%rel_tracking_charge
  endif

  call track1 (orb_0r, ele_r, ele_r%branch%param, orb_1r)
  call make_mat6 (ele_r, ele_r%branch%param, orb_0r)

  if (ele_f%key == rfcavity$) then
    dorb%vec = orb_1r%vec - orb_1f%vec
    dorb%vec(5) = (orb_1r%vec(5) - orb_0r%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5)) 
    dorb%t = c_light * ((orb_1r%t - orb_0r%t) - (orb_1f%t - orb_0f%t))

  else
    orb_1r%vec(2) = -orb_1r%vec(2)
    orb_1r%vec(4) = -orb_1r%vec(4)
    dorb%vec = orb_1r%vec - orb_0f%vec
    dorb%vec(5) = (orb_1r%vec(5) - orb_0r%vec(5)) - (orb_1f%vec(5) - orb_0f%vec(5)) 
    dorb%t = c_light * ((orb_1r%t - orb_0r%t) - (orb_1f%t - orb_0f%t))
  endif

  if (verbosity) then
     print '(a, 6es12.4, 5x, es12.4)', '1: ', orb_0r%vec, orb_0r%t
     print '(a, 6es12.4, 5x, es12.4)', '2: ', orb_1r%vec, orb_1r%t
     print *
     print '(a, 6es12.4, 5x, es12.4)', 'D: ', dorb%vec(1:6), dorb%t
  end if

  !

  if (ele_f%key == rfcavity$) then
    mat_f = ele_f%mat6
  else
    m(1, 1:6) = [1,  0,  0,  0,  0,  0]
    m(2, 1:6) = [0, -1,  0,  0,  0,  0]
    m(3, 1:6) = [0,  0,  1,  0,  0,  0]
    m(4, 1:6) = [0,  0,  0, -1,  0,  0]
    m(5, 1:6) = [0,  0,  0,  0, -1,  0]
    m(6, 1:6) = [0,  0,  0,  0,  0,  1]
    call mat_inverse (ele_f%mat6, mat_f)
    mat_f = matmul(matmul(m, mat_f), m)
  endif

  !

  str = '"' // trim(ele_f%name) // ':'
  length = len(trim(ele_f%name))

  if (i > 1) write (1, *)
  write (1, '(a, a, a, es16.8)') str(1:length+2), 'dorb(1)"', tolerance(trim(ele_f%name),'dorb(1)'), dorb%vec(1)
  write (1, '(a, a, a, es16.8)') str(1:length+2), 'dorb(2)"', tolerance(trim(ele_f%name),'dorb(2)'), dorb%vec(2)
  write (1, '(a, a, a, es16.8)') str(1:length+2), 'dorb(3)"', tolerance(trim(ele_f%name),'dorb(3)'), dorb%vec(3)
  write (1, '(a, a, a, es16.8)') str(1:length+2), 'dorb(4)"', tolerance(trim(ele_f%name),'dorb(4)'), dorb%vec(4)
  write (1, '(a, a, a, es16.8)') str(1:length+2), 'dorb(5)"', tolerance(trim(ele_f%name),'dorb(5)'), dorb%vec(5)
  write (1, '(a, a, a, es16.8)') str(1:length+2), 'dorb(6)"', tolerance(trim(ele_f%name),'dorb(6)'), dorb%vec(6)
  write (1, '(a, a, a, es16.8)') str(1:length+2), 'c*dt"   ', tolerance(trim(ele_f%name),'c*dt'), dorb%t

  if (verbosity) then
     write (1, *)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_f(1,:)" ABS 2e-14 ', ele_f%mat6(1,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_f(2,:)" ABS 2e-14 ', ele_f%mat6(2,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_f(3,:)" ABS 2e-14 ', ele_f%mat6(3,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_f(4,:)" ABS 2e-14 ', ele_f%mat6(4,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_f(5,:)" ABS 2e-14 ', ele_f%mat6(5,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_f(6,:)" ABS 2e-14 ', ele_f%mat6(6,:)

     write (1, *)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_r(1,:)" ABS 2e-14 ', ele_r%mat6(1,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_r(2,:)" ABS 2e-14 ', ele_r%mat6(2,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_r(3,:)" ABS 2e-14 ', ele_r%mat6(3,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_r(4,:)" ABS 2e-14 ', ele_r%mat6(4,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_r(5,:)" ABS 2e-14 ', ele_r%mat6(5,:)
     write (1, '(a, a, 6es16.8)') str(1:length+2), 'mat_r(6,:)" ABS 2e-14 ', ele_r%mat6(6,:)
  end if

  write (1, *)
  write (1, '(a, a, a, 6es16.8)') str(1:length+2), 'dmat(1,:)"', tolerance(trim(ele_f%name),'dmat(1,:)'), ele_r%mat6(1,:) - mat_f(1,:)
  write (1, '(a, a, a, 6es16.8)') str(1:length+2), 'dmat(2,:)"', tolerance(trim(ele_f%name),'dmat(2,:)'), ele_r%mat6(2,:) - mat_f(2,:)
  write (1, '(a, a, a, 6es16.8)') str(1:length+2), 'dmat(3,:)"', tolerance(trim(ele_f%name),'dmat(3,:)'), ele_r%mat6(3,:) - mat_f(3,:)
  write (1, '(a, a, a, 6es16.8)') str(1:length+2), 'dmat(4,:)"', tolerance(trim(ele_f%name),'dmat(4,:)'), ele_r%mat6(4,:) - mat_f(4,:)
  write (1, '(a, a, a, 6es16.8)') str(1:length+2), 'dmat(5,:)"', tolerance(trim(ele_f%name),'dmat(5,:)'), ele_r%mat6(5,:) - mat_f(5,:)
  write (1, '(a, a, a, 6es16.8)') str(1:length+2), 'dmat(6,:)"', tolerance(trim(ele_f%name),'dmat(6,:)'), ele_r%mat6(6,:) - mat_f(6,:)

  if (verbosity) then
     write (1, *)
     write (1, '(a, a,  3es16.8)') str(1:length+2), 'max(dvec, dmat)" ', maxval(abs([dorb%vec, dorb%t])), &
           maxval(abs(ele_r%mat6 - mat_f))
  end if

END DO

! And close

close (1)

!-------------------------------------------------------------------------
contains

character(11) function tolerance(instr1,instr2)
character(*) :: instr1
character(*) :: instr2

tolerance = ' ABS 2e-14 '

select case (instr1)
  case('SOL_QUAD2')
    select case (instr2)
    case('dorb(1)') ;   tolerance = ' ABS 1e-13 '
    case('dorb(4)') ;   tolerance = ' ABS 1e-13 '
    case('dmat(1,:)') ; tolerance = ' ABS 1e-11 '
    case('dmat(2,:)') ; tolerance = ' ABS 1e-11 '
    case('dmat(3,:)') ; tolerance = ' ABS 1e-11 '
    case('dmat(4,:)') ; tolerance = ' ABS 1e-11 '
    case('dmat(5,:)') ; tolerance = ' ABS 1e-12 '
    end select

  case('SBEND4')
    select case (instr2)
    case('dorb(1)') ;   tolerance = ' ABS 2e-12 '
    case('dorb(3)') ;   tolerance = ' ABS 2e-11 '
    case('dorb(5)') ;   tolerance = ' ABS 2e-08 '
    case('c*dt') ;      tolerance = ' ABS 2e-08 '
    case('dmat(1,:)') ; tolerance = ' ABS 2e-09 '
    case('dmat(3,:)') ; tolerance = ' ABS 2e-08 '
    case('dmat(5,:)') ; tolerance = ' ABS 2e-08 '
    end select
end select

end function tolerance

end program
