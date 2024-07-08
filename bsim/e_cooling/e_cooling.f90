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

integer i, ie, ib, n_loc, n, n_turns, status
logical err_flag
character(200) master_input_file

!

namelist / params / ec, bmad_com, beam_init, n_turns
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

do i = 1, size(ec%wd%xm%s)
  if (ec%wd%xm%s(i) == real_garbage$) exit
enddo
ec%wd%n_size_array = i - 1

! TODO: make this a subroutine to be cleaner
! Make splines
! First, fill in the information in the splines
do i = 1, ec%wd%n_size_array
  ec%wd%xm%A_spl(i)%x0 = ec%wd%xm%s(i)
  ec%wd%xm%k_spl(i)%x0 = ec%wd%xm%s(i)
  ec%wd%xm%lambda_spl(i)%x0 = ec%wd%xm%s(i)

  ec%wd%xm%A_spl(i)%y0 = ec%wd%xm%A(i)
  ec%wd%xm%k_spl(i)%y0 = ec%wd%xm%k(i)
  ec%wd%xm%lambda_spl(i)%y0 = ec%wd%xm%lambda(i)


  ec%wd%ym%A_spl(i)%x0 = ec%wd%ym%s(i)
  ec%wd%ym%k_spl(i)%x0 = ec%wd%ym%s(i)
  ec%wd%ym%lambda_spl(i)%x0 = ec%wd%ym%s(i)

  ec%wd%ym%A_spl(i)%y0 = ec%wd%ym%A(i)
  ec%wd%ym%k_spl(i)%y0 = ec%wd%ym%k(i)
  ec%wd%ym%lambda_spl(i)%y0 = ec%wd%ym%lambda(i)


  ec%wd%xk%A_spl(i)%x0 = ec%wd%xk%s(i)
  ec%wd%xk%k_spl(i)%x0 = ec%wd%xk%s(i)
  ec%wd%xk%lambda_spl(i)%x0 = ec%wd%xk%s(i)
  ec%wd%xk%D_h_spl(i)%x0 = ec%wd%xk%s(i)
  ec%wd%xk%D_e11_spl(i)%x0 = ec%wd%xk%s(i)
  ec%wd%xk%D_e12_spl(i)%x0 = ec%wd%xk%s(i)
  ec%wd%xk%D_e22_spl(i)%x0 = ec%wd%xk%s(i)

  ec%wd%xk%A_spl(i)%y0 = ec%wd%xk%A(i)
  ec%wd%xk%k_spl(i)%y0 = ec%wd%xk%k(i)
  ec%wd%xk%lambda_spl(i)%y0 = ec%wd%xk%lambda(i)
  ec%wd%xk%D_h_spl(i)%y0 = ec%wd%xk%D_h(i)
  ec%wd%xk%D_e11_spl(i)%y0 = ec%wd%xk%D_e11(i)
  ec%wd%xk%D_e12_spl(i)%y0 = ec%wd%xk%D_e12(i)
  ec%wd%xk%D_e22_spl(i)%y0 = ec%wd%xk%D_e22(i)


  ec%wd%yk%A_spl(i)%x0 = ec%wd%yk%s(i)
  ec%wd%yk%k_spl(i)%x0 = ec%wd%yk%s(i)
  ec%wd%yk%lambda_spl(i)%x0 = ec%wd%yk%s(i)
  ec%wd%yk%D_h_spl(i)%x0 = ec%wd%yk%s(i)
  ec%wd%yk%D_e11_spl(i)%x0 = ec%wd%yk%s(i)
  ec%wd%yk%D_e12_spl(i)%x0 = ec%wd%yk%s(i)
  ec%wd%yk%D_e22_spl(i)%x0 = ec%wd%yk%s(i)

  ec%wd%yk%A_spl(i)%y0 = ec%wd%yk%A(i)
  ec%wd%yk%k_spl(i)%y0 = ec%wd%yk%k(i)
  ec%wd%yk%lambda_spl(i)%y0 = ec%wd%yk%lambda(i)
  ec%wd%yk%D_h_spl(i)%y0 = ec%wd%yk%D_h(i)
  ec%wd%yk%D_e11_spl(i)%y0 = ec%wd%yk%D_e11(i)
  ec%wd%yk%D_e12_spl(i)%y0 = ec%wd%yk%D_e12(i)
  ec%wd%yk%D_e22_spl(i)%y0 = ec%wd%yk%D_e22(i)
enddo

! Fill rest with zeros, and make x0 locations equally spaced.
do i = ec%wd%n_size_array, size(ec%wd%xm%s)
  ec%wd%xm%A_spl(i)%x0 = 2*ec%wd%xm%A_spl(i-1)%x0 - ec%wd%xm%A_spl(i-2)%x0
  ec%wd%xm%k_spl(i)%x0 = 2*ec%wd%xm%k_spl(i-1)%x0 - ec%wd%xm%k_spl(i-2)%x0
  ec%wd%xm%lambda_spl(i)%x0 = 2*ec%wd%xm%lambda_spl(i-1)%x0 - ec%wd%xm%lambda_spl(i-2)%x0

  ec%wd%xm%A_spl(i)%y0 = 0
  ec%wd%xm%k_spl(i)%y0 = 0
  ec%wd%xm%lambda_spl(i)%y0 = 0


  ec%wd%ym%A_spl(i)%x0 = 2*ec%wd%ym%A_spl(i-1)%x0 - ec%wd%ym%A_spl(i-2)%x0
  ec%wd%ym%k_spl(i)%x0 = 2*ec%wd%ym%k_spl(i-1)%x0 - ec%wd%ym%k_spl(i-2)%x0
  ec%wd%ym%lambda_spl(i)%x0 = 2*ec%wd%ym%lambda_spl(i-1)%x0 - ec%wd%ym%lambda_spl(i-2)%x0

  ec%wd%ym%A_spl(i)%y0 = 0
  ec%wd%ym%k_spl(i)%y0 = 0
  ec%wd%ym%lambda_spl(i)%y0 = 0


  ec%wd%xk%A_spl(i)%x0 = 2*ec%wd%xk%A_spl(i-1)%x0 - ec%wd%xk%A_spl(i-2)%x0
  ec%wd%xk%k_spl(i)%x0 = 2*ec%wd%xk%k_spl(i-1)%x0 - ec%wd%xk%k_spl(i-2)%x0
  ec%wd%xk%lambda_spl(i)%x0 = 2*ec%wd%xk%lambda_spl(i-1)%x0 - ec%wd%xk%lambda_spl(i-2)%x0
  ec%wd%xk%D_h_spl(i)%x0 = 2*ec%wd%xk%D_h_spl(i-1)%x0 - ec%wd%xk%D_h_spl(i-2)%x0
  ec%wd%xk%D_e11_spl(i)%x0 = 2*ec%wd%xk%D_e11_spl(i-1)%x0 - ec%wd%xk%D_e11_spl(i-2)%x0
  ec%wd%xk%D_e12_spl(i)%x0 = 2*ec%wd%xk%D_e12_spl(i-1)%x0 - ec%wd%xk%D_e12_spl(i-2)%x0
  ec%wd%xk%D_e22_spl(i)%x0 = 2*ec%wd%xk%D_e22_spl(i-1)%x0 - ec%wd%xk%D_e22_spl(i-2)%x0

  ec%wd%xk%A_spl(i)%y0 = 0
  ec%wd%xk%k_spl(i)%y0 = 0
  ec%wd%xk%lambda_spl(i)%y0 = 0
  ec%wd%xk%D_h_spl(i)%y0 = 0
  ec%wd%xk%D_e11_spl(i)%y0 = 0
  ec%wd%xk%D_e12_spl(i)%y0 = 0
  ec%wd%xk%D_e22_spl(i)%y0 = 0


  ec%wd%yk%A_spl(i)%x0 = 2*ec%wd%yk%A_spl(i-1)%x0 - ec%wd%yk%A_spl(i-2)%x0
  ec%wd%yk%k_spl(i)%x0 = 2*ec%wd%yk%k_spl(i-1)%x0 - ec%wd%yk%k_spl(i-2)%x0
  ec%wd%yk%lambda_spl(i)%x0 = 2*ec%wd%yk%lambda_spl(i-1)%x0 - ec%wd%yk%lambda_spl(i-2)%x0
  ec%wd%yk%D_h_spl(i)%x0 = 2*ec%wd%yk%D_h_spl(i-1)%x0 - ec%wd%yk%D_h_spl(i-2)%x0
  ec%wd%yk%D_e11_spl(i)%x0 = 2*ec%wd%yk%D_e11_spl(i-1)%x0 - ec%wd%yk%D_e11_spl(i-2)%x0
  ec%wd%yk%D_e12_spl(i)%x0 = 2*ec%wd%yk%D_e12_spl(i-1)%x0 - ec%wd%yk%D_e12_spl(i-2)%x0
  ec%wd%yk%D_e22_spl(i)%x0 = 2*ec%wd%yk%D_e22_spl(i-1)%x0 - ec%wd%yk%D_e22_spl(i-2)%x0

  ec%wd%yk%A_spl(i)%y0 = 0
  ec%wd%yk%k_spl(i)%y0 = 0
  ec%wd%yk%lambda_spl(i)%y0 = 0
  ec%wd%yk%D_h_spl(i)%y0 = 0
  ec%wd%yk%D_e11_spl(i)%y0 = 0
  ec%wd%yk%D_e12_spl(i)%y0 = 0
  ec%wd%yk%D_e22_spl(i)%y0 = 0
enddo

! Now actually do the interpolation.
call spline_akima(ec%wd%xm%A_spl, err_flag)
call spline_akima(ec%wd%xm%k_spl, err_flag)
call spline_akima(ec%wd%xm%lambda_spl, err_flag)

call spline_akima(ec%wd%ym%A_spl, err_flag)
call spline_akima(ec%wd%ym%k_spl, err_flag)
call spline_akima(ec%wd%ym%lambda_spl, err_flag)

call spline_akima(ec%wd%xk%A_spl, err_flag)
call spline_akima(ec%wd%xk%k_spl, err_flag)
call spline_akima(ec%wd%xk%lambda_spl, err_flag)
call spline_akima(ec%wd%xk%D_h_spl, err_flag)
call spline_akima(ec%wd%xk%D_e11_spl, err_flag)
call spline_akima(ec%wd%xk%D_e12_spl, err_flag)
call spline_akima(ec%wd%xk%D_e22_spl, err_flag)

call spline_akima(ec%wd%yk%A_spl, err_flag)
call spline_akima(ec%wd%yk%k_spl, err_flag)
call spline_akima(ec%wd%yk%lambda_spl, err_flag)
call spline_akima(ec%wd%yk%D_h_spl, err_flag)
call spline_akima(ec%wd%yk%D_e11_spl, err_flag)
call spline_akima(ec%wd%yk%D_e12_spl, err_flag)
call spline_akima(ec%wd%yk%D_e22_spl, err_flag)

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

do n = 1, n_turns
  do ie = 1, branch%n_ele_track
    call track1_bunch(beam%bunch(1), branch%ele(ie), err_flag)
  enddo
enddo

! Print end distribution

print *
print *, 'Final (partial) beam distribution:'
do ib = 1, size(bunch%particle)
  p => bunch%particle(ib)
  if(ib < 10) then
    print '(i6, 6es16.8)', ib, p%vec
  endif
enddo

end program
