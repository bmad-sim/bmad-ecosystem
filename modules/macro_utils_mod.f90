module macro_utils_mod

use bmad
use macroparticle_mod

type bunch_param_struct
  real(rp) beta, alpha, gamma
  real(rp) sigma, p_sigma
  real(rp) emitt ! normalized emittance
endtype

type bunch_params_struct
  type (bunch_param_struct) :: x, y
  real(rp) dpz_dz, z_sigma, p_z_sigma
  type (coord_struct) :: centroid
  real(rp) :: charge ! bunch charge NOT including lost particles
endtype

type beam_params_struct
  type (bunch_params_struct), pointer :: bunch(:) => null()
endtype

contains

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_params (bunch, ele, params)
!
! Finds all the bunch parameters defined in bunch_params_struct
!
! Input:
!  bunch      -- bunch_struct
!  ele        -- ele_struct: element to find parameters at
!
! Output:
!  params     -- bunch_params_struct
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

Subroutine calc_bunch_params (bunch, ele, params)

implicit none

type (bunch_struct) bunch
type (ele_struct) ele
type (bunch_params_struct) params

call calc_bunch_centroid (bunch, params)
call calc_bunch_twiss_and_emittance (bunch, ele, params, .false.)
call calc_bunch_sigma (bunch, params, .false.)
call calc_bunch_dpz_dz (bunch, params, .false.)

end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_twiss_and_emittance (bunch, ele, params, calc_centroid)
! 
! Subroutine to calculate normalized x and y emittance and twiss
! parameters of a bunch at element ele. If calc_centroid if false then the
! centroid will not be calculated. If calc_centroid is true then the bunch
! centroid and charge will calculated with calc_bunch_centroid.
! The default is calc_centroid = .true.
!
! Modules needed:
!   use macro_utils_mod
!
! Input:
!   bunch         -- bunch_struct: bunch for which emittance should be calculated.
!   ele           -- ele_struct: Element at which to calculate emittance.
!   calc_centroid -- (optional) Logical: calculate the bunch centroid and charge?
!
! Output:   
!   params  -- bunch_params_struct
!               %x%alpha; %x%beta; %x%gamma
!               %x%emitt
!               %y%alpha; %y%beta; %y%gamma
!               %y%emitt
!               %centroid (optional)
!               %charge (optional)
!-

subroutine calc_bunch_twiss_and_emittance (bunch, ele, params, calc_centroid)

implicit none

type (bunch_struct), target :: bunch
type(ele_struct) :: ele
type (bunch_params_struct) :: params
type (macro_struct), pointer :: macro(:), mp
logical, optional :: calc_centroid

integer i, j
real(rp) x_norm_emit, bunch_x, bunch_xp, avg_x2, avg_xxp, avg_xp2
real(rp) y_norm_emit, bunch_y, bunch_yp, avg_y2, avg_yyp, avg_yp2
real(rp) charge, avg_energy
real(rp) x, xp, s11, s12, s22, s66, x_eta, x_etap
real(rp) y, yp, s33, s34, s44, y_eta, y_etap
real(rp) zp, s16, s26, s36, s46

! initialize 

  bunch_x      = 0
  bunch_xp     = 0
  bunch_y      = 0
  bunch_yp     = 0
  avg_energy   = 0

  x_eta        = ele%x%eta
  x_etap       = ele%x%etap
  y_eta        = ele%y%eta
  y_etap       = ele%y%etap
  
  avg_x2       = 0
  avg_xxp      = 0
  avg_xp2      = 0
  avg_y2       = 0
  avg_yyp      = 0
  avg_yp2      = 0
  
  ! first calculate average x/x'and y/y' of bunch weighted by charge
  
  if (present(calc_centroid)) then
    if (calc_centroid) call calc_bunch_centroid (bunch, params)
  else
    call calc_bunch_centroid (bunch, params)
  endif

  if (params%charge == 0) then
    params%x%beta  = 0.0
    params%x%gamma = 0.0
    params%x%alpha = 0.0
    params%y%beta  = 0.0
    params%y%gamma = 0.0
    params%y%alpha = 0.0
    params%x%emitt = 0.0
    params%y%emitt = 0.0
    return
  endif

  bunch_x  = params%centroid%vec(1)
  bunch_xp = params%centroid%vec(2)
  bunch_y  = params%centroid%vec(3)
  bunch_yp = params%centroid%vec(4)
  

  ! now calculate emittance
  
  do i=1, size(bunch%slice)
  
    do j = 1, size(bunch%slice(i)%macro)

      mp => bunch%slice(i)%macro(j)
      if (mp%lost) cycle
  
      zp = mp%r%vec(6)
      charge = mp%charge

      s11 = mp%sigma(s11$)
      s12 = mp%sigma(s12$)
      s22 = mp%sigma(s22$)
      s33 = mp%sigma(s33$)
      s34 = mp%sigma(s34$)
      s44 = mp%sigma(s44$)
      s66 = mp%sigma(s66$)
      s16 = mp%sigma(s16$)
      s26 = mp%sigma(s26$)
      s36 = mp%sigma(s36$)
      s46 = mp%sigma(s46$)
  
      x  = mp%r%vec(1) - bunch_x  - x_eta *  zp
      xp = mp%r%vec(2) - bunch_xp - x_etap * zp
      y  = mp%r%vec(3) - bunch_y  - y_eta *  zp
      yp = mp%r%vec(4) - bunch_yp - y_etap * zp
  
      avg_x2  = avg_x2 + charge * &
               (s11 + x_eta**2 * s66 + x**2 - 2 * s16*x_eta)
      avg_xp2 = avg_xp2+ charge * &
             (s22 + x_etap**2 * s66 + xp**2 - 2 * s26*x_etap)
      avg_xxp = avg_xxp+ charge * &
               (s12 + x_eta*x_etap*s66 + x*xp - s16*x_etap - s26*x_eta)
      avg_y2  = avg_y2 + charge * &
               (s33 + y_eta**2 * s66 + y**2 - 2 * s36*y_eta)
      avg_yp2 = avg_yp2 + charge * &
             (s44 + y_etap**2 * s66 + yp**2 - 2 * s46*y_etap)
      avg_yyp = avg_yyp + charge * &
               (s34 + y_eta*y_etap*s66 + y*yp - s36*y_etap - s46*y_eta)
      avg_energy = avg_energy + charge * (1+mp%r%vec(6)) * ele%value(beam_energy$)
  
    end do
  end do
  
  avg_x2  = avg_x2  / params%charge
  avg_xp2 = avg_xp2 / params%charge
  avg_xxp = avg_xxp / params%charge
  avg_y2  = avg_y2  / params%charge
  avg_yp2 = avg_yp2 / params%charge
  avg_yyp = avg_yyp / params%charge
  avg_energy = avg_energy / params%charge
  
  params%x%emitt = (avg_energy/m_electron) * sqrt(avg_x2*avg_xp2 - avg_xxp**2)
  params%y%emitt = (avg_energy/m_electron) * sqrt(avg_y2*avg_yp2 - avg_yyp**2)

  params%x%beta  = avg_x2   / params%x%emitt
  params%x%gamma = avg_xp2  / params%x%emitt
  params%x%alpha = -avg_xxp / params%x%emitt
  
  params%y%beta  = avg_y2   / params%y%emitt
  params%y%gamma = avg_yp2  / params%y%emitt
  params%y%alpha = -avg_yyp / params%y%emitt
    
end subroutine

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_centroid (bunch, params)
!
! Finds the bunch centroid and charge, keeping track of lost macroparticles
!
! Input:
!  bunch         -- bunch_struct
!
! Output:
!  params      -- bunch_params_struct
!                  %centorid  -- bunch centroid
!                  %charge    -- bunch charge excluding lost macroparticles
!  
!-

subroutine calc_bunch_centroid (bunch, params)

implicit none

type (bunch_struct), target :: bunch
type (bunch_params_struct) :: params
type (macro_struct), pointer :: macro(:)
type (coord_struct) :: centroid

real(rp) tot_charge

integer i, j

  centroid%vec = 0.0
  tot_charge = 0.0

  do i = 1, size(bunch%slice)
    macro => bunch%slice(i)%macro
    do j = 1, 6
      centroid%vec(j) = centroid%vec(j) + sum (macro%charge*macro%r%vec(j), &
                                                       mask = .not. macro%lost)
    enddo
    tot_charge = tot_charge + sum(macro%charge, mask = .not. macro%lost)
  enddo

  params%charge = tot_charge
    
  if (tot_charge .eq. 0) then
    params%centroid%vec = 0.0
  else
    params%centroid%vec = centroid%vec / tot_charge
  endif

end subroutine calc_bunch_centroid
  
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_sigma (bunch, params, calc_centroid)
!
! Finds the bunch sigmas. If calc_centroid if false then the
! centroid will not be calculated. If calc_centorid is true then the bunch
! centroid and charge will calculated with calc_bunch_centroid.
! The default is calc_centroid = .true.

!
! Input:
!  bunch        -- bunch_struct
!
! Output:
!  params       -- bunch_params_struct
!                  %x%sigma; %x%p_sigma
!                  %y%sigma; %y%p_sigma
!                  %z%sigma; %z%p_sigma
!                  %centroid (optional)
!                  %charge (optional)
!
!-

subroutine calc_bunch_sigma (bunch, params, calc_centroid)

implicit none

type (bunch_struct), target :: bunch
type (bunch_params_struct) :: params
logical, optional :: calc_centroid

type (macro_struct), pointer :: macro(:)
type (coord_struct) :: centroid
real(rp) var1(6), var2(6)

integer i, j, ix
integer, parameter :: s_plane(6) = (/ s11$, s22$, s33$, s44$, s55$, s66$ /)

  var1 = 0
  var2 = 0
 
  if (present(calc_centroid)) then
    if (calc_centroid) call calc_bunch_centroid (bunch, params)
  else
    call calc_bunch_centroid (bunch, params)
  endif

  if (params%charge .eq. 0) then
    params%x%sigma   = 0.0
    params%x%p_sigma = 0.0
    params%y%sigma   = 0.0
    params%y%p_sigma = 0.0
    params%z_sigma   = 0.0
    params%p_z_sigma = 0.0
    return
  endif

  ! This is the macroparticle centroid part (computes variance (sigma**2))
  
  do i = 1, size(bunch%slice)
    macro => bunch%slice(i)%macro
    do j = 1, 6
      var1(j) = var1(j) + sum(macro%charge*(macro%r%vec(j) - params%centroid%vec(j))**2, &
                                               mask = .not. macro%lost)
    enddo
  enddo
 
  ! This is the macroparticle sigma part (again, the variance)
 
  do i = 1, size(bunch%slice)
    macro => bunch%slice(i)%macro
    do j = 1, 6
      var2(j) = var2(j) + sum(macro%charge*macro%sigma(s_plane(j)), mask = .not. macro%lost)
    enddo
  enddo

  ! Now take root of sum of variances to get net sigma, that is,
  !     (sigma**2 = sigma1**2 + sigma2**2
  params%x%sigma   = sqrt((var1(1) + var2(1))/params%charge)
  params%x%p_sigma = sqrt((var1(2) + var2(2))/params%charge)
  params%y%sigma   = sqrt((var1(3) + var2(3))/params%charge)
  params%y%p_sigma = sqrt((var1(4) + var2(4))/params%charge)
  params%z_sigma   = sqrt((var1(5) + var2(5))/params%charge)
  params%p_z_sigma = sqrt((var1(6) + var2(6))/params%charge)
  
end subroutine calc_bunch_sigma
 
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine calc_bunch_dpz_dz (bunch, params, calc_centroid)
!
! Finds the energy-Z correlation for a bunch. If calc_centorid is true then
! the bunch centroid and charge will calculated with calc_bunch_centroid.
! The default is calc_centroid = .true.

!
! Input:
!  bunch        -- bunch_struct
!
! Output        -- bunch_params_struct
!
!-

subroutine calc_bunch_dpz_dz (bunch, params, calc_centroid)

implicit none
  
type (bunch_struct), target :: bunch
type (bunch_params_struct) :: params
logical, optional :: calc_centroid

type (macro_struct), pointer :: macro(:)

real(rp) z, zp, ave_zz, ave_zzp
integer*4 i

  if (present(calc_centroid)) then
    if (calc_centroid) call calc_bunch_centroid (bunch, params)
  else
    call calc_bunch_centroid (bunch, params)
  endif

  if (params%charge .eq. 0) then
    params%dpz_dz = 0.0
    return
  endif

  ave_zzp = 0.0
  ave_zz = 0.0
  
  do i = 1,size(bunch%slice)
    macro => bunch%slice(i)%macro
!    do k = 1,size(bunch%slice(1)%macro)
!       z =  bunch%slice(j)%macro(k)%r%vec(5) - params%centroid%vec(5)
!       zp =  bunch%slice(j)%macro(k)%r%vec(6) - params%centroid%vec(6)
!       ave_zzp = ave_zzp + bunch%slice(j)%macro(k)%charge* &
!            (bunch%slice(j)%macro(k)%sigma(s56$) + z * zp)
!       ave_zz = ave_zz + bunch%slice(j)%macro(k)%charge* &
!            (bunch%slice(j)%macro(k)%sigma(s55$) + z * z)

    ave_zzp = ave_zzp + sum(macro%charge*(macro%sigma(s56$) + &
         	            (macro%r%vec(5) - params%centroid%vec(5)) * &
                            (macro%r%vec(6) - params%centroid%vec(6))), &
    	                    mask = .not. macro%lost)
         
    ave_zz = ave_zz + sum(macro%charge*(macro%sigma(s55$) + &
            	          (macro%r%vec(5) - params%centroid%vec(5)) * &
                          (macro%r%vec(5) - params%centroid%vec(5))), &
		          mask = .not. macro%lost)
	     
!    end do
  enddo
  
  ave_zzp = ave_zzp / bunch%charge
  ave_zz = ave_zz / bunch%charge
  params%dpz_dz = ave_zzp / ave_zz
  
end subroutine calc_bunch_dpz_dz

end module
