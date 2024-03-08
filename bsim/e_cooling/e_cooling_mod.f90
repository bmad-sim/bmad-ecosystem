module e_cooling_mod

use bmad

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
  integer :: n_size_array = 0         ! Calculated by program
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
  character(200) lat_file, wake_and_diffusion_file
  type (e_cooling_wd_struct) :: wd
  type (beam_init_struct) beam_init
end type

type e_cooling_common_struct
  type (lat_struct) lat
  type (branch_struct), pointer :: branch   ! Pointer to particular branch bunch is tracked through.
  type (ele_struct), pointer :: cool_ele    ! Feedback element.
  type (ele_struct), pointer :: input_ele   ! Pickup element hooked to feedback ele.
  type (ele_struct), pointer :: output_ele  ! Kicker element hooked to feedback ele.
end type

! Common 

type (e_cooling_param_struct), target :: ec            ! Parameters set by User.
type (e_cooling_common_struct), target :: ec_com       ! Not set by User.

end module
