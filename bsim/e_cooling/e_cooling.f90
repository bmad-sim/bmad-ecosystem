!+
! Program e_cooling
!
! Program to simulate coherent electron cooling
!-

program e_cooling

use bmad

implicit none

! Structs for coherent electron cooling

type e_cooling_wake_struct
  real(rp) :: s = real_garbage$
  real(rp) :: A = real_garbage$
  real(rp) :: k = real_garbage$
  real(rp) :: lambda = real_garbage$
end type 

type e_cooling_diffusion_struct
  real(rp) :: s = real_garbage$
  real(rp) :: A = real_garbage$
  real(rp) :: k = real_garbage$
  real(rp) :: lambda = real_garbage$
  real(rp) :: D_h = real_garbage$
  real(rp) :: D_e11 = real_garbage$
  real(rp) :: D_e12 = real_garbage$
  real(rp) :: D_e22 = real_garbage$
end type

type e_cooling_wd_struct
  integer :: n_size = 0
  real(rp) :: Ie_peak = 0             ! Peak electron current.
  real(rp) :: sigma_ze = 0            ! Electron bunch length.
  real(rp) :: supergaussian_order = 0 ! Order of supergaussian used to model the longitudinal 
                                      !   electron bunch distribution
  real(rp) :: off_E_reduction = 0     ! How much to reduce the wake for off-energy protons
  real(rp) :: phi_avg = 0             ! Average phase advance of electrons through a single amplifier straight
  type (e_cooling_wake_struct) :: xm(100), ym(100)
  type (e_cooling_diffusion_struct) :: xk(100), yk(100)
end type


!

type e_cooling_param_struct
  character(200) master_input_file, lat_file, wake_and_diffusion_file
  type (e_cooling_wd_struct) :: wd
end type

type e_cooling_common_struct
  type (lat_struct) lat
end type

character(200) master_input_file, lat_file, wake_and_diffusion_file

type (e_cooling_param_struct), target :: ec
type (e_cooling_common_struct), target :: ec_com

integer i

!

namelist / params / ec, bmad_com
namelist / wake_and_diffusion / ec

! Read master_input_file

ec%master_input_file = ''
call get_command_argument(1, ec%master_input_file)
if (ec%master_input_file == '') ec%master_input_file = 'e_cooling.init'
print '(2a)', 'Initialization file: ', trim(ec%master_input_file)

open (1, file = ec%master_input_file, status = 'old', action = 'read')
read (1, nml = params)
close(1)

! Read wake_and_diffusion_file

if (ec%wake_and_diffusion_file /= '') then
  open (1, file = ec%wake_and_diffusion_file, status = 'old', action = 'read')
  read (1, nml = wake_and_diffusion)
  close(1)
endif

do i = 1, size(ec%wd%xm)
  if (ec%wd%xm(i)%s == real_garbage$) exit
enddo
ec%wd%n_size = i - 1

!

call bmad_parser (ec%lat_file, ec_com%lat)

end program
