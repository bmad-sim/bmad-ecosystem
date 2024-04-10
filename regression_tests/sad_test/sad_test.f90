program sad_test

use bmad

implicit none

type (lat_struct), target :: lat
type (ele_struct), pointer :: ele
type (coord_struct), allocatable :: orbit(:)

integer nargs
logical exist, debug_mode
character(100) str

!

debug_mode = .false.
nargs = command_argument_count()

if (nargs > 0) then
  call get_command_argument(1, str)
  debug_mode = .true.
endif

!

if (.not. debug_mode) then
  inquire (file = '../../util_programs/sad_to_bmad/sad_to_bmad.py', exist = exist)
  if (exist) then
    call system_command ('python3 ../../util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_ptc.params')
  else
    call system_command ('python3 $ACC_ROOT_DIR/util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_ptc.params')
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
    call system_command ('python3 ../../util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_bmad.params')
  else
    call system_command ('python3 $ACC_ROOT_DIR/util_programs/sad_to_bmad/sad_to_bmad.py sad_to_bmad_bmad.params')
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
