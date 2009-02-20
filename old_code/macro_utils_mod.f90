module macro_utils_mod

use bmad
use macroparticle_mod

!type macro_bunch_lat_param_struct
!  real(rp) beta, alpha, gamma
!  real(rp) sigma, p_sigma
!  real(rp) emit ! normalized emittance
!endtype

!type macro_bunch_params_struct
!  type (macro_bunch_lat_param_struct) :: x, y
!  real(rp) dpz_dz, z_sigma, p_z_sigma
!  type (coord_struct) :: centroid
!  real(rp) :: charge ! bunch charge NOT including lost particles
!endtype

!type macro_beam_params_struct
!  type (macro_bunch_params_struct), pointer :: bunch(:) => null()
!endtype

type macro_bunch_lat_param_struct
  real(rp) beta, alpha, gamma
  real(rp) eta, etap
  real(rp) sigma, p_sigma
  real(rp) dpx_dx ! x x' correlation
  real(rp) norm_emit ! normalized emittance
end type

type macro_bunch_params_struct
  type (macro_bunch_lat_param_struct) :: x, y, z, a, b
  type (coord_struct) :: centroid !Lab frame
  real(rp) :: charge ! bunch charge NOT including lost particles
end type

type macro_beam_params_struct
  type (macro_bunch_params_struct), pointer :: bunch(:) => null()
end type

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_macro_bunch_params (bunch, ele, params)
!
! Finds all macroparticle bunch parameters defined in macro_bunch_params_struct, 
! both normal-mode and projected parameters found.
!
! Modules needed:
!  use beam_mod
!
! Input:
!  bunch     -- Macro_unch_struct
!  ele       -- ele_struct: element to find parameters at
!
! Output     -- Macro_bunch_params_struct
!                 %a%alpha; %a%beta; %a%gamma
!                 %a%sigma; %a%p_sigma
!                 %a%emit;  %a%dpx_dx
!                 %b%alpha; %b%beta; %b%gamma
!                 %b%sigma; %b%p_sigma
!                 %b%emit;  %b%dpx_dx
!                 %z%alpha; %z%beta; %z%gamma
!                 %z%sigma; %z%p_sigma
!                 %z%emit;  %z%dpx_dx
!                 %a%alpha; %a%beta; %a%gamma
!                 %a%sigma; %a%p_sigma
!                 %a%emit;  %a%dpx_dx
!                 %b%alpha; %b%beta; %b%gamma
!                 %b%sigma; %b%p_sigma
!                 %b%emit;  %b%dpx_dx
!                 %centroid
!                 %n_part
!
!-

subroutine calc_macro_bunch_params (bunch, ele, params)

  implicit none

  type (macro_bunch_struct), intent(in) :: bunch
  type (macro_beam_struct), save :: a_beam ! converted to normal modes
  type (ele_struct) :: ele
  type (macro_bunch_params_struct) params
  type (macro_struct), pointer :: macro(:)

  real(rp) exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
  real(rp) avg_energy 
  real(rp) avg_delta, exp_delta2, var1, var2
  real(rp) v_mat(4,4), v_inv_mat(4,4), sigma_a(4,4)

  integer i, ii, j
  
  ! centroid, charge, average energy, delta spread and delta center
  params%charge = 0.0
  params%centroid%vec = 0.0
  avg_energy = 0.0
  avg_delta = 0.0
  do i = 1, size(bunch%slice)
    macro => bunch%slice(i)%macro
    do j = 1, 6
      params%centroid%vec(j) = params%centroid%vec(j) + sum (macro%charge*macro%r%vec(j), &
                                                       mask = .not. macro%lost)
    enddo
    avg_energy = avg_energy + sum((1+macro%r%vec(6))*macro%charge, mask = .not. macro%lost)
    avg_delta = avg_delta + sum(macro%r%vec(6)*macro%charge, mask = .not. macro%lost)
    params%charge = params%charge + sum(macro%charge, mask = .not. macro%lost)
  enddo

  params%charge = params%charge
  params%centroid%vec = params%centroid%vec / params%charge
  avg_energy = avg_energy * ele%value(E_TOT$) / params%charge
  avg_delta = avg_delta  / params%charge
    
  if (params%charge < e_charge) then
    params%centroid%vec = 0.0
    call zero_plane (params%a)
    call zero_plane (params%b)
    call zero_plane (params%z)
    call zero_plane (params%a)
    call zero_plane (params%b)
    return
  endif

  ! delta spread
  ! <delta**2>_tot = sqrt ( s_66**2 + <delta**2>**2)
  exp_delta2 = 0.0
  var1 = 0.0
  var2 = 0.0
  do i = 1, size(bunch%slice)
    macro => bunch%slice(i)%macro
    var1 = var1 + sum (macro%sigma(s66$)*macro%charge, mask = .not. macro%lost)
    var2 = var2 + sum (macro%charge*(macro%r%vec(6) - params%centroid%vec(6))**2, &
                                                            mask = .not. macro%lost)
  enddo
  exp_delta2  = sqrt((var1/params%charge)**2 + (var2/params%charge)**2)
  
  ! Projected Parameters
  ! X
  call find_expectations (bunch, 'x', exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  call param_stuffit (params%a, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
     
  ! Y
  call find_expectations (bunch, 'y', exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  call param_stuffit (params%b, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
  
  ! Z
  call find_expectations (bunch, 'z', exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  call param_stuffit (params%z, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
  
  !***
  ! Normal-Mode Parameters
  
  ! take out coupling
  call reallocate_macro_beam (a_beam, 1, size(bunch%slice), size(bunch%slice(1)%macro))
  a_beam%bunch(1) = bunch
  call make_v_mats (ele, v_mat, v_inv_mat)
  do i = 1, size(a_beam%bunch(1)%slice)
    macro => a_beam%bunch(1)%slice(i)%macro
    do ii = 1, size(macro)
      macro(ii)%r%vec(1:4) = matmul(v_inv_mat, macro(ii)%r%vec(1:4))
      call mp_sigma_to_mat (macro(ii)%sigma, sigma_a)
      sigma_a = matmul (transpose(v_mat), matmul (sigma_a, v_mat))
      call mat_to_mp_sigma (sigma_a, macro(ii)%sigma)
    enddo
  enddo 
  
  ! A
  call find_expectations (a_beam%bunch(1), 'a', exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  call param_stuffit (params%a, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
     
  ! B
  call find_expectations (a_beam%bunch(1), 'b', exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  call param_stuffit (params%b, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
     
  
contains
!----------------------------------------------------------------------
subroutine zero_plane (param)

  implicit none

  type (macro_bunch_lat_param_struct), intent(out) :: param

  param%beta       = 0.0
  param%alpha      = 0.0
  param%gamma      = 0.0
  param%eta        = 0.0
  param%etap       = 0.0
  param%sigma      = 0.0
  param%p_sigma    = 0.0
  param%dpx_dx     = 0.0
  param%norm_emit  = 0.0

end subroutine zero_plane
  
!----------------------------------------------------------------------
subroutine find_expectations (bunch, plane, exp_x2, exp_p_x2, exp_x_p_x, &
                              exp_x_d, exp_px_d)

  implicit none

  character(1), intent(in) :: plane
  type (macro_bunch_struct), intent(in) :: bunch
  type (macro_struct), pointer :: macro(:)
  real(rp), intent(out) ::  exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
  real(rp) ::  var_x2, var_p_x2, var_x_p_x, var_x_d, var_px_d
  real(rp) avg_x, avg_p_x, eta, etap
  real(rp), allocatable, save :: x(:), p_x(:), d(:)
  integer p11, p22, p12, p16, p26
  
  integer i, index

  logical normal_mode_flag

  if (.not. allocated(x)) allocate (x(size(bunch%slice(1)%macro)), &
                  p_x(size(bunch%slice(1)%macro)), d(size(bunch%slice(1)%macro)))
  if (size(x) .ne. size(bunch%slice(1)%macro)) then
    deallocate(x, p_x, d)
    allocate (x(size(bunch%slice(1)%macro)), &
                  p_x(size(bunch%slice(1)%macro)), d(size(bunch%slice(1)%macro)))

  endif
  
  if (plane .eq. 'x') then
    index = 1
    normal_mode_flag = .false.
    p11 = s11$
    p22 = s22$
    p12 = s12$
    p16 = s16$
    p26 = s26$
  elseif (plane .eq. 'y') then
    index = 3
    normal_mode_flag = .false.
    p11 = s33$
    p22 = s44$
    p12 = s34$
    p16 = s36$
    p26 = s46$
  elseif (plane .eq. 'z') then
    index = 5
    normal_mode_flag = .false.
    p11 = s55$
    p22 = s66$
    p12 = s56$
    p16 = s56$
    p26 = s66$
  elseif (plane .eq. 'a') then
    index = 1
    normal_mode_flag = .true.
    p11 = s11$
    p22 = s22$
    p12 = s12$
    p16 = s16$
    p26 = s26$
  elseif (plane .eq. 'b') then
    index = 3
    normal_mode_flag = .true.
    p11 = s33$
    p22 = s44$
    p12 = s34$
    p16 = s36$
    p26 = s46$
  endif
  
  avg_x = params%centroid%vec(index)
  avg_p_x = params%centroid%vec(index+1)

  exp_x_d   = 0.0
  exp_px_d  = 0.0
  exp_x2    = 0.0
  exp_p_x2  = 0.0
  exp_x_p_x = 0.0
  var_x_d   = 0.0
  var_px_d  = 0.0
  var_x2    = 0.0
  var_p_x2  = 0.0
  var_x_p_x = 0.0

  ! <xy>_tot = sqrt ( s_xy**2 + <xy>**2)
  
  ! centroid part
  do i = 1, size(bunch%slice)
    macro => bunch%slice(i)%macro

    x   = macro%r%vec(index) - avg_x 
    p_x = macro%r%vec(index+1) - avg_p_x
    d   = macro%r%vec(6) - avg_delta

    var_x_d   = var_x_d + sum(x*d*macro%charge, mask = .not. macro%lost)
    var_px_d  = var_px_d + sum(p_x*d*macro%charge, mask = .not. macro%lost)
    
    var_x2    = var_x2 + sum(x**2*macro%charge, mask = .not. macro%lost)
    var_p_x2  = var_p_x2 + sum(p_x**2*macro%charge, mask = .not. macro%lost)
    var_x_p_x = var_x_p_x + sum(x*p_x*macro%charge, mask = .not. macro%lost)
  enddo

  
  ! sigma part
  do i = 1, size(bunch%slice)
    macro => bunch%slice(i)%macro

    exp_x_d   = exp_x_d + sum(macro%sigma(p16)*macro%charge, mask = .not. macro%lost)
    exp_px_d  = exp_px_d + sum(macro%sigma(p26)*macro%charge, mask = .not. macro%lost)

    exp_x2    = exp_x2 + sum(macro%sigma(p11)*macro%charge, mask = .not. macro%lost)
    exp_p_x2  = exp_p_x2 + sum(macro%sigma(p22)*macro%charge, mask = .not. macro%lost)
    exp_x_p_x = exp_x_p_x + sum(macro%sigma(p12)*macro%charge, mask = .not. macro%lost)
  enddo

  ! add sigma and centroid parts in quadrature
  exp_x_d   = sqrt((exp_x_d   /params%charge)**2 + (var_x_d   / params%charge)**2 ) 
  exp_px_d  = sqrt((exp_px_d  /params%charge)**2 + (var_px_d  / params%charge)**2 ) 
  exp_x2    = sqrt((exp_x2    /params%charge)**2 + (var_x2    / params%charge)**2 )  
  exp_p_x2  = sqrt((exp_p_x2  /params%charge)**2 + (var_p_x2  / params%charge)**2 ) 
  exp_x_p_x = sqrt((exp_x_p_x /params%charge)**2 + (var_x_p_x / params%charge)**2 ) 
 
  if (normal_mode_flag) then
    eta   = exp_x_d / exp_delta2
    etap  = exp_px_d / exp_delta2

    exp_x2    = exp_x2 - 2*eta*exp_x_d + (eta**2)*exp_delta2
    exp_p_x2  = exp_p_x2 - 2*etap*exp_px_d + (etap**2)*exp_delta2
    exp_x_p_x = exp_x_p_x - etap*exp_x_d - eta*exp_px_d + eta*etap*exp_delta2
  endif

end subroutine find_expectations

! contains
!----------------------------------------------------------------------
subroutine param_stuffit (param, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  implicit none

  type (macro_bunch_lat_param_struct), intent(out) :: param
  real(rp), intent(in) :: exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
  real(rp) emit

  if (exp_x2*exp_p_x2 < exp_x_p_x**2) then
    emit = 0.0
    param%alpha = 0.0
    param%beta  = 0.0
    param%gamma = 0.0
    param%eta   = 0.0
    param%etap  = 0.0
    param%norm_emit = 0.0
    param%sigma = 0.0
    param%p_sigma = 0.0
    param%dpx_dx = 0.0
    return
  else
    emit = SQRT(exp_x2*exp_p_x2 - exp_x_p_x**2)
  endif

  param%alpha = exp_x_p_x / emit
  param%beta  = exp_x2 / emit
  param%gamma = exp_p_x2 / emit
  
  param%eta   = exp_x_d / exp_delta2
  param%etap  = exp_px_d / exp_delta2

  param%norm_emit = (avg_energy/m_electron) * emit

  param%sigma = SQRT(exp_x2)
  param%p_sigma = SQRT(exp_p_x2)

  param%dpx_dx = exp_x_p_x / exp_x2

end subroutine param_stuffit

end subroutine calc_macro_bunch_params
  
 end module
