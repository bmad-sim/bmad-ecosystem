module macro_utils_mod

use bmad
use macroparticle_mod

!type macro_bunch_param_struct
!  real(rp) beta, alpha, gamma
!  real(rp) sigma, p_sigma
!  real(rp) emitt ! normalized emittance
!endtype

!type macro_bunch_params_struct
!  type (macro_bunch_param_struct) :: x, y
!  real(rp) dpz_dz, z_sigma, p_z_sigma
!  type (coord_struct) :: centroid
!  real(rp) :: charge ! bunch charge NOT including lost particles
!endtype

!type macro_beam_params_struct
!  type (macro_bunch_params_struct), pointer :: bunch(:) => null()
!endtype

type macro_bunch_param_struct
  real(rp) beta, alpha, gamma
  real(rp) eta, etap
  real(rp) sigma, p_sigma
  real(rp) dpx_dx ! x x' correlation
  real(rp) norm_emitt ! normalized emittance
endtype

type macro_bunch_params_struct
  type (macro_bunch_param_struct) :: x, y, z, a, b
  type (coord_struct) :: centroid !Lab frame
  real(rp) :: charge ! bunch charge NOT including lost particles
endtype

type macro_beam_params_struct
  type (macro_bunch_params_struct), pointer :: bunch(:) => null()
endtype

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
!                 %x%alpha; %x%beta; %x%gamma
!                 %x%sigma; %x%p_sigma
!                 %x%emitt; %x%dpx_dx
!                 %y%alpha; %y%beta; %y%gamma
!                 %y%sigma; %y%p_sigma
!                 %y%emitt; %y%dpx_dx
!                 %z%alpha; %z%beta; %z%gamma
!                 %z%sigma; %z%p_sigma
!                 %z%emitt; %z%dpx_dx
!                 %a%alpha; %a%beta; %a%gamma
!                 %a%sigma; %a%p_sigma
!                 %a%emitt; %a%dpx_dx
!                 %b%alpha; %b%beta; %b%gamma
!                 %b%sigma; %b%p_sigma
!                 %b%emitt; %b%dpx_dx
!                 %centroid
!                 %n_part
!
!-

subroutine calc_macro_bunch_params (bunch, ele, params)

  implicit none

  type (macro_bunch_struct), intent(in) :: bunch
  type (macro_beam_struct) :: a_beam ! converted to normal modes
  type (ele_struct), intent(in) :: ele
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
  avg_energy = avg_energy * ele%value(beam_energy$) / params%charge
  avg_delta = avg_delta  / params%charge
    
  if (params%charge .lt. e_charge) then
    params%centroid%vec = 0.0
    call zero_plane (params%x)
    call zero_plane (params%y)
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
    var1 = sum (macro%sigma(s66$)*macro%charge, mask = .not. macro%lost)
    var2 = sum (macro%charge*(macro%r%vec(6) - params%centroid%vec(6))**2, &
                                                            mask = .not. macro%lost)
  enddo
  exp_delta2  = sqrt((var1/params%charge)**2 + (var2/params%charge)**2)
  
  ! Projected Parameters
  ! X
  call find_expectations (bunch, 'x', exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  call param_stuffit (params%x, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
     
  ! Y
  call find_expectations (bunch, 'y', exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)

  call param_stuffit (params%y, exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d)
  
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

  type (macro_bunch_param_struct), intent(out) :: param

  param%beta       = 0.0
  param%alpha      = 0.0
  param%gamma      = 0.0
  param%eta        = 0.0
  param%etap       = 0.0
  param%sigma      = 0.0
  param%p_sigma    = 0.0
  param%dpx_dx     = 0.0
  param%norm_emitt = 0.0

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
  real(rp), avg_x, avg_p_x, eta, etap
  real(rp), automatic, dimension(size(bunch%slice(1)%macro)) :: x, p_x, d
  integer p11, p22, p12, p16, p26
  
  integer i, index

  logical normal_mode_flag

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

  type (macro_bunch_param_struct), intent(out) :: param
  real(rp), intent(in) :: exp_x2, exp_p_x2, exp_x_p_x, exp_x_d, exp_px_d
  real(rp) emitt

  emitt = SQRT(exp_x2*exp_p_x2 - exp_x_p_x**2)

  param%alpha = exp_x_p_x / emitt
  param%beta  = exp_x2 / emitt
  param%gamma = exp_p_x2 / emitt
  
  param%eta   = exp_x_d / exp_delta2
  param%etap  = exp_px_d / exp_delta2

  param%norm_emitt = (avg_energy/m_electron) * emitt

  param%sigma = SQRT(exp_x2)
  param%p_sigma = SQRT(exp_p_x2)

  param%dpx_dx = exp_x_p_x / exp_x2

end subroutine param_stuffit

end subroutine calc_macro_bunch_params
  
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_macro_bunch_params (bunch, ele, params)
!
! Finds all the bunch parameters defined in bunch_params_struct
!
! Input:
!  bunch      -- macro_bunch_struct
!  ele        -- ele_struct: element to find parameters at
!
! Output:
!  params     -- macro_bunch_params_struct
!                 %x%alpha; %x%beta; %x%gamma
!                 %x%sigma; %x%p_sigma
!                 %x%emitt
!                 %y%alpha; %y%beta; %y%gamma
!                 %y%sigma; %y%p_sigma
!                 %y%emitt
!                 %dpz_dz; %z_sigma; %p_z_sigma
!                 %centroid; %charge
!
!-
!
!Subroutine calc_macro_bunch_params (bunch, ele, params)
!
!implicit none
!
!type (macro_bunch_struct) bunch
!type (ele_struct) ele
!type (macro_bunch_params_struct) params
!
!call calc_macro_bunch_centroid (bunch, params)
!call calc_macro_bunch_twiss_and_emittance (bunch, ele, params, .false.)
!call calc_macro_bunch_sigma (bunch, params, .false.)
!call calc_macro_bunch_dpz_dz (bunch, params, .false.)
!
!end subroutine
!
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!+
!! Subroutine calc_macro_bunch_twiss_and_emittance (bunch, ele, params, calc_centroid)
!! 
!! Subroutine to calculate normalized x and y emittance and twiss
!! parameters of a bunch at element ele. If calc_centroid if false then the
!! centroid will not be calculated. If calc_centroid is true then the bunch
!! centroid and charge will calculated with calc_macro_bunch_centroid.
!! The default is calc_centroid = .true.
!!
!! Modules needed:
!!   use macro_utils_mod
!!
!! Input:
!!   bunch         -- macro_bunch_struct: bunch for which emittance should be calculated.
!!   ele           -- ele_struct: Element at which to calculate emittance.
!!   calc_centroid -- (optional) Logical: calculate the bunch centroid and charge?
!!
!! Output:   
!!   params  -- macro_bunch_params_struct
!!               %x%alpha; %x%beta; %x%gamma
!!               %x%emitt
!!               %y%alpha; %y%beta; %y%gamma
!!               %y%emitt
!!               %centroid (optional)
!!               %charge (optional)
!!-
!!
!!subroutine calc_macro_bunch_twiss_and_emittance (bunch, ele, params, calc_centroid)
!!
!implicit none
!
!type (macro_bunch_struct), target :: bunch
!type(ele_struct) :: ele
!type (macro_bunch_params_struct) :: params
!type (macro_struct), pointer :: macro(:), mp
!logical, optional :: calc_centroid
!
!integer i, j
!real(rp) x_norm_emit, bunch_x, bunch_xp, avg_x2, avg_xxp, avg_xp2
!real(rp) y_norm_emit, bunch_y, bunch_yp, avg_y2, avg_yyp, avg_yp2
!real(rp) charge, avg_energy
!real(rp) x, xp, s11, s12, s22, s66, x_eta, x_etap
!real(rp) y, yp, s33, s34, s44, y_eta, y_etap
!real(rp) zp, s16, s26, s36, s46
!real(rp) v_mat(4,4), v_inv_mat(4,4), orb(4), sigma_a(4,4)
!
!! initialize 
!
!  bunch_x      = 0
!  bunch_xp     = 0
!  bunch_y      = 0
!  bunch_yp     = 0
!  avg_energy   = 0
!
!  x_eta        = ele%x%eta
!  x_etap       = ele%x%etap
!  y_eta        = ele%y%eta
!  y_etap       = ele%y%etap
!  
!  ! testing....
!  x_eta        = ele%x%eta_lab
!  x_etap       = ele%x%etap_lab
!  y_eta        = ele%y%eta_lab
!  y_etap       = ele%y%etap_lab
!
!! x_eta        = 0.0
!! x_etap       = 0.0
!! y_eta        = 0.0
!! y_etap       = 0.0
!!  ! testing...
!!      
!!  avg_x2       = 0
!!  avg_xxp      = 0
!!  avg_xp2      = 0
!  avg_y2       = 0
!  avg_yyp      = 0
!  avg_yp2      = 0
!  
!  ! first calculate average x/x'and y/y' of bunch weighted by charge
!  
!  if (present(calc_centroid)) then
!    if (calc_centroid) call calc_macro_bunch_centroid (bunch, params)
!  else
!    call calc_macro_bunch_centroid (bunch, params)
!  endif
!
!  if (params%charge == 0) then
!    params%x%beta  = 0.0
!    params%x%gamma = 0.0
!    params%x%alpha = 0.0
!    params%y%beta  = 0.0
!    params%y%gamma = 0.0
!    params%y%alpha = 0.0
!    params%x%emitt = 0.0
!    params%y%emitt = 0.0
!    return
!  endif
!
!  bunch_x  = params%centroid%vec(1)
!  bunch_xp = params%centroid%vec(2)
!  bunch_y  = params%centroid%vec(3)
!  bunch_yp = params%centroid%vec(4)
!  
!
!  ! now calculate emittance
!  
!  do i=1, size(bunch%slice)
!  
!    do j = 1, size(bunch%slice(i)%macro)
!
!      mp => bunch%slice(i)%macro(j)
!      if (mp%lost) cycle
!  
!      zp = mp%r%vec(6)
!      charge = mp%charge
!
!      s11 = mp%sigma(s11$)
!      s12 = mp%sigma(s12$)
!      s22 = mp%sigma(s22$)
!      s33 = mp%sigma(s33$)
!      s34 = mp%sigma(s34$)
!      s44 = mp%sigma(s44$)
!      s66 = mp%sigma(s66$)
!      s16 = mp%sigma(s16$)
!      s26 = mp%sigma(s26$)
!      s36 = mp%sigma(s36$)
!      s46 = mp%sigma(s46$)
!  
!      ! testing...
!!     call make_v_mats (ele, v_mat, v_inv_mat)
!!     
!!     call mp_sigma_to_mat (mp%sigma, sigma_a)
!!     sigma_a = matmul (transpose(v_mat), matmul (sigma_a, v_mat))
!
!!     s11 = sigma_a(1,1)
!!     s12 = sigma_a(1,2) 
!!     s22 = sigma_a(2,2) 
!!     s33 = sigma_a(3,3) 
!!     s34 = sigma_a(3,4) 
!!     s44 = sigma_a(4,4) 
!!     s66 = mp%sigma(s66$)
!!     s16 = x_eta * s66
!!     s26 = x_etap* s66
!!     s36 = y_eta * s66
!!     s46 = y_etap* s66
!
!      ! testing...
!      
!      x  = mp%r%vec(1) - bunch_x  - x_eta *  zp
!      xp = mp%r%vec(2) - bunch_xp - x_etap * zp
!      y  = mp%r%vec(3) - bunch_y  - y_eta *  zp
!      yp = mp%r%vec(4) - bunch_yp - y_etap * zp
!  
!      
!      ! testing....
!!     orb = mp%r%vec(1:4) - (/ bunch_x, bunch_xp, bunch_y, bunch_yp /)
!!     orb = matmul(v_inv_mat,orb)
!!     x  = orb(1) - zp*x_eta
!!     xp = orb(2) - zp*x_etap
!!     y  = orb(3) - zp*y_eta
!!     yp = orb(4) - zp*y_etap
!!      ! testing....
!!
!!      
!!      avg_x2  = avg_x2 + charge * &
!               (s11 + x_eta**2 * s66 + x**2 - 2 * s16*x_eta)
!      avg_xp2 = avg_xp2+ charge * &
!             (s22 + x_etap**2 * s66 + xp**2 - 2 * s26*x_etap)
!      avg_xxp = avg_xxp+ charge * &
!               (s12 + x_eta*x_etap*s66 + x*xp - s16*x_etap - s26*x_eta)
!      avg_y2  = avg_y2 + charge * &
!               (s33 + y_eta**2 * s66 + y**2 - 2 * s36*y_eta)
!      avg_yp2 = avg_yp2 + charge * &
!             (s44 + y_etap**2 * s66 + yp**2 - 2 * s46*y_etap)
!      avg_yyp = avg_yyp + charge * &
!               (s34 + y_eta*y_etap*s66 + y*yp - s36*y_etap - s46*y_eta)
!      avg_energy = avg_energy + charge * (1+mp%r%vec(6)) * ele%value(beam_energy$)
!  
!    end do
!  end do
!  
!  avg_x2  = avg_x2  / params%charge
!  avg_xp2 = avg_xp2 / params%charge
!  avg_xxp = avg_xxp / params%charge
!  avg_y2  = avg_y2  / params%charge
!  avg_yp2 = avg_yp2 / params%charge
!  avg_yyp = avg_yyp / params%charge
!  avg_energy = avg_energy / params%charge
!  
!  params%x%emitt = sqrt(avg_x2*avg_xp2 - avg_xxp**2)
!  params%y%emitt = sqrt(avg_y2*avg_yp2 - avg_yyp**2)
!
!  params%x%beta  = avg_x2   / params%x%emitt
!  params%x%gamma = avg_xp2  / params%x%emitt
!  params%x%alpha = -avg_xxp / params%x%emitt
!  
!  params%y%beta  = avg_y2   / params%y%emitt
!  params%y%gamma = avg_yp2  / params%y%emitt
!  params%y%alpha = -avg_yyp / params%y%emitt
!    
!  ! now normalize the emittance
!  params%x%emitt = (avg_energy/m_electron) * params%x%emitt
!  params%y%emitt = (avg_energy/m_electron) * params%y%emitt
!
!end subroutine
!
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!+
!! Subroutine calc_macro_bunch_centroid (bunch, params)
!!
!! Finds the bunch centroid and charge, keeping track of lost macroparticles
!!
!! Input:
!!  bunch         -- macro_bunch_struct
!!
!! Output:
!!  params      -- macro_bunch_params_struct
!!                  %centorid  -- bunch centroid
!!                  %charge    -- bunch charge excluding lost macroparticles
!!  
!!-
!
!subroutine calc_macro_bunch_centroid (bunch, params)
!
!implicit none
!
!type (macro_bunch_struct), target :: bunch
!type (macro_bunch_params_struct) :: params
!type (macro_struct), pointer :: macro(:)
!type (coord_struct) :: centroid
!
!real(rp) tot_charge
!
!integer i, j
!
!  centroid%vec = 0.0
!  tot_charge = 0.0
!
!  do i = 1, size(bunch%slice)
!    macro => bunch%slice(i)%macro
!    do j = 1, 6
!      centroid%vec(j) = centroid%vec(j) + sum (macro%charge*macro%r%vec(j), &
!                                                       mask = .not. macro%lost)
!    enddo
!    tot_charge = tot_charge + sum(macro%charge, mask = .not. macro%lost)
!  enddo
!
!  params%charge = tot_charge
!    
!  if (tot_charge .eq. 0) then
!    params%centroid%vec = 0.0
!  else
!    params%centroid%vec = centroid%vec / tot_charge
!  endif
!
!end subroutine calc_macro_bunch_centroid
!  
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!+
!! Subroutine calc_macro_bunch_sigma (bunch, params, calc_centroid)
!!
!! Finds the bunch sigmas. If calc_cenroid if false then the
!! centroid will not be calculated. If calc_centorid is true then the bunch
!! centroid and charge will calculated with calc_macro_bunch_centroid.
!! The default is calc_centroid = .true.
!
!!
!! Input:
!!  bunch        -- macro_bunch_struct
!!
!! Output:
!!  params       -- macro_bunch_params_struct
!!                  %x%sigma; %x%p_sigma
!!                  %y%sigma; %y%p_sigma
!!                  %z%sigma; %z%p_sigma
!!                  %centroid (optional)
!!                  %charge (optional)
!!
!!-
!
!subroutine calc_macro_bunch_sigma (bunch, params, calc_centroid)
!
!implicit none
!
!type (macro_bunch_struct), target :: bunch
!type (macro_bunch_params_struct) :: params
!logical, optional :: calc_centroid
!
!type (macro_struct), pointer :: macro(:)
!type (coord_struct) :: centroid
!real(rp) var1(6), var2(6)
!
!integer i, j, ix
!integer, parameter :: s_plane(6) = (/ s11$, s22$, s33$, s44$, s55$, s66$ /)
!
!  var1 = 0
!  var2 = 0
! 
!  if (present(calc_centroid)) then
!    if (calc_centroid) call calc_macro_bunch_centroid (bunch, params)
!  else
!    call calc_macro_bunch_centroid (bunch, params)
!  endif
!
!  if (params%charge .eq. 0) then
!    params%x%sigma   = 0.0
!    params%x%p_sigma = 0.0
!    params%y%sigma   = 0.0
!    params%y%p_sigma = 0.0
!    params%z_sigma   = 0.0
!    params%p_z_sigma = 0.0
!    return
!  endif
!
!  ! This is the macroparticle centroid part (computes variance (sigma**2))
!  
!  do i = 1, size(bunch%slice)
!    macro => bunch%slice(i)%macro
!    do j = 1, 6
!      var1(j) = var1(j) + sum(macro%charge*(macro%r%vec(j) - params%centroid%vec(j))**2, &
!                                               mask = .not. macro%lost)
!    enddo
!  enddo
! 
!  ! This is the macroparticle sigma part (again, the variance)
! 
!  do i = 1, size(bunch%slice)
!    macro => bunch%slice(i)%macro
!    do j = 1, 6
!      var2(j) = var2(j) + sum(macro%charge*macro%sigma(s_plane(j)), mask = .not. macro%lost)
!    enddo
!  enddo
!
!  ! Now take root of sum of variances to get net sigma, that is,
!  !     (sigma**2 = sigma1**2 + sigma2**2
!  params%x%sigma   = sqrt((var1(1) + var2(1))/params%charge)
!  params%x%p_sigma = sqrt((var1(2) + var2(2))/params%charge)
!  params%y%sigma   = sqrt((var1(3) + var2(3))/params%charge)
!  params%y%p_sigma = sqrt((var1(4) + var2(4))/params%charge)
!  params%z_sigma   = sqrt((var1(5) + var2(5))/params%charge)
!  params%p_z_sigma = sqrt((var1(6) + var2(6))/params%charge)
!  
!end subroutine calc_macro_bunch_sigma
! 
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!--------------------------------------------------------------------------
!!+
!! Subroutine calc_macro_bunch_dpz_dz (bunch, params, calc_centroid)
!!
!! Finds the energy-Z correlation for a bunch. If calc_centorid is true then
!! the bunch centroid and charge will calculated with calc_macro_bunch_centroid.
!! The default is calc_centroid = .true.
!
!!
!! Input:
!!  bunch        -- macro_bunch_struct
!!
!! Output        -- macro_bunch_params_struct
!!
!!-
!
!subroutine calc_macro_bunch_dpz_dz (bunch, params, calc_centroid)
!
!implicit none
!  
!type (macro_bunch_struct), target :: bunch
!type (macro_bunch_params_struct) :: params
!logical, optional :: calc_centroid
!
!type (macro_struct), pointer :: macro(:)
!
!real(rp) z, zp, ave_zz, ave_zzp
!integer*4 i
!
!  if (present(calc_centroid)) then
!    if (calc_centroid) call calc_macro_bunch_centroid (bunch, params)
!  else
!    call calc_macro_bunch_centroid (bunch, params)
!  endif
!
!  if (params%charge .eq. 0) then
!    params%dpz_dz = 0.0
!    return
!  endif
!
!  ave_zzp = 0.0
!  ave_zz = 0.0
!  
!  do i = 1,size(bunch%slice)
!    macro => bunch%slice(i)%macro
!
!    ave_zzp = ave_zzp + sum(macro%charge*(macro%sigma(s56$) + &
!         	            (macro%r%vec(5) - params%centroid%vec(5)) * &
!                            (macro%r%vec(6) - params%centroid%vec(6))), &
!    	                    mask = .not. macro%lost)
!         
!    ave_zz = ave_zz + sum(macro%charge*(macro%sigma(s55$) + &
!            	          (macro%r%vec(5) - params%centroid%vec(5)) * &
!                          (macro%r%vec(5) - params%centroid%vec(5))), &
!		          mask = .not. macro%lost)
!	     
!  enddo
!  
!  ave_zzp = ave_zzp / bunch%charge
!  ave_zz = ave_zz / bunch%charge
!  params%dpz_dz = ave_zzp / ave_zz
!  
!end subroutine calc_macro_bunch_dpz_dz
!
 end module
