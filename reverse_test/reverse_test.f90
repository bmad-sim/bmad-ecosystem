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

do i = 1, lat%n_ele_max - 1

  orb_0f = lat%beam_start
  ele_f => lat%branch(0)%ele(i)
  call init_coord (orb_0f, orb_0f%vec, ele_f, upstream_end$)

  call track1 (orb_0f, ele_f, ele_f%branch%param, orb_1f)

  str = trim(ele_f%name)

  if (verbosity) then
    print *, str
    print '(a, 6es12.4, 5x, es12.4)', '0: ', orb_0f%vec, orb_0f%t
    print '(a, 6es12.4, 5x, es12.4)', '1: ', orb_1f%vec, orb_1f%t
  end if

  orb_0r           = orb_1f
  if (ele_f%key /= elseparator$) then
    orb_0r%species   = flip_species_charge(orb_0r%species)
  endif
  orb_0r%direction = -1

  call track1(orb_0r, ele_f, ele_f%branch%param, orb_1r)

  dorb%vec = orb_1r%vec - orb_0f%vec
  dorb%vec = orb_1r%vec - orb_0f%vec
  dorb%t   = orb_1r%t + orb_0f%t - 2 * orb_1f%t

  if (verbosity) then
     print '(a, 6es12.4, 5x, es12.4)', '1: ', orb_0r%vec, orb_0r%t
     print '(a, 6es12.4, 5x, es12.4)', '2: ', orb_1r%vec, orb_1r%t
     print *
     print '(a, 6es12.4, 5x, es12.4)', 'D: ', dorb%vec(1:6), dorb%t
  end if

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

end do

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
