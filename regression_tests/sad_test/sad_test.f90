program sad_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orbit(:)

integer nargs
logical exist, existpy3, debug_mode
character(100) str, pythontouse

pythontouse = 'python'

!

debug_mode = .false.
nargs = command_argument_count()

if (nargs > 0) then
  call get_command_argument(1, str)
  debug_mode = .true.
endif

!
existpy3 = .false.

inquire (file = 'python3', exist = exist)
if (exist) then
  existpy3 = .true.
endif

!if (existpy3) then
!   print *,'exists 1'
!endif

inquire (file = '/usr/bin/python3', exist = exist)
if (exist) then
  existpy3 = .true.
endif

!if (existpy3) then
!   print *,'exists 2'
!endif

inquire (file = '/usr/common/python/python3', exist = exist)
if (exist) then
  existpy3 = .true.
endif

!if (existpy3) then
!   print *,'exists 3'
!endif

if (existpy3) then
  pythontouse = 'python3'
endif

if (.not. debug_mode) then
  inquire (file = '../../util_programs/sad_to_bmad/sad_to_bmad.py', exist = exist)
  if (exist) then
!    call system_command ('python ../../util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_ptc.params')
    call system_command (pythontouse//' ../../util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_ptc.params')
  else
!    call system_command ('python $ACC_ROOT_DIR/util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_ptc.params')
    call system_command (pythontouse//' $ACC_ROOT_DIR/util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_ptc.params')
  endif
endif

call bmad_parser('sler_1689.bmad', lat)
call twiss_and_track(lat, orbit)

ele => lat%ele(0)

open (1, file = 'output.now')
write (1, '(a, 6es16.8)') '"Orb0-PTC"  ABS 1E-10', orbit(0)%vec
write (1, '(a, 2es16.8)') '"Betas-PTC"  REL 1E-6', ele%a%beta, ele%b%beta
write (1, '(a, 2es16.8)') '"Alphas-PTC" ABS 1E-8', ele%a%alpha, ele%b%alpha

!

if (.not. debug_mode) then
  inquire (file = '../../util_programs/sad_to_bmad/sad_to_bmad.py', exist = exist)
  if (exist) then
!    call system_command ('python ../../util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_bmad.params')
    call system_command (pythontouse//' ../../util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_bmad.params')
  else
!    call system_command ('python $ACC_ROOT_DIR/util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_bmad.params')
    call system_command (pythontouse//' $ACC_ROOT_DIR/util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_bmad.params')
  endif
endif

call bmad_parser('sler_1689.bmad', lat)
call twiss_and_track(lat, orbit)

ele => lat%ele(0)

open (1, file = 'output.now')
write (1, '(a, 6es16.8)') '"Orb0-Bmad"  ABS 1E-10', orbit(0)%vec
write (1, '(a, 2es16.8)') '"Betas-Bmad"  REL 1E-6', ele%a%beta, ele%b%beta
write (1, '(a, 2es16.8)') '"Alphas-Bmad" ABS 1E-8', ele%a%alpha, ele%b%alpha

close (1)

end program
