!+
! Program e_cooling
!
! Program to simulate coherent electron cooling
!-

program e_cooling

use e_cooling_mod
use beam_mod
use beam_utils

implicit none

type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (coord_struct), pointer :: p
type (branch_struct), pointer :: branch
type (coord_struct), allocatable :: closed_orb(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (beam_init_struct) beam_init
procedure (track1_bunch_hook_def) :: track1_bunch_hook

integer i, ie, ib, n_loc, status
logical err_flag
character(200) master_input_file

!

namelist / params / ec, bmad_com, beam_init
namelist / wake_and_diffusion / ec

!

track1_bunch_hook_ptr => track1_bunch_hook

! Read master_input_file

master_input_file = ''
call get_command_argument(1, master_input_file)
if (master_input_file == '') master_input_file = 'e_cooling.init'
print '(2a)', 'Initialization file: ', trim(master_input_file)

open (1, file = master_input_file, status = 'old', action = 'read')
read (1, nml = params)

ec%beam_init = beam_init

! Read wake_and_diffusion_file. If not set, this file is same as the master_input_file

if (ec%wake_and_diffusion_file /= '') then
  close(1)
  open (1, file = ec%wake_and_diffusion_file, status = 'old', action = 'read')
endif

read (1, nml = wake_and_diffusion)
close(1)

do i = 1, size(ec%wd%xm)
  if (ec%wd%xm(i)%s == real_garbage$) exit
enddo
ec%wd%n_size_array = i - 1

! Parse lattice

call bmad_parser (ec%lat_file, ec_com%lat)
ec_com%branch => ec_com%lat%branch(0)
branch => ec_com%branch           ! Typing shortcut
call twiss_and_track(ec_com%lat, closed_orb, status, branch%ix_branch)

! Find cooler element and input, output elements.

call lat_ele_locator('FEEDBACK::*', ec_com%lat, eles, n_loc, err_flag)
if (n_loc == 0) then
  print *, 'CANNOT FIND FEEDBACK ELEMENT!'
  stop
endif
if (n_loc > 1) then
  print *, 'MULTIPLE FEEDBACK ELEMENTS! I DO NOT KNOW WHAT TO DO!'
  stop
endif

ec_com%cool_ele   => eles(1)%ele
ec_com%input_ele  => pointer_to_slave(ec_com%cool_ele, 1)
ec_com%output_ele => pointer_to_slave(ec_com%cool_ele, 2)

! Test track through lattice

call init_beam_distribution(branch%ele(0), ec_com%lat%param, ec%beam_init, beam, err_flag)
bunch => beam%bunch(1)  ! Only track 1 bunch

do ie = 1, branch%n_ele_track
  call track1_bunch(beam%bunch(1), branch%ele(ie), err_flag)
enddo

! Print end distribution

print *
print *, 'Final beam distribution:'
do ib = 1, size(bunch%particle)
  p => bunch%particle(ib)
  print '(i6, 6es12.4)', ib, p%vec
enddo

end program
