#include "CESR_platform.inc"

MODULE touschek_mod

USE bmad_struct
USE bmad_interface
USE nr

TYPE touschek_common_struct
  REAL(rp) B1
  REAL(rp) B2
  REAL(rp) tau_m
END TYPE

! This common block is needed in order to pass parameters to
! integrand.  We use NR's qromb function, which requires that the
! function it integrates be one variable only.

TYPE(touschek_common_struct), private :: touschek_com

PUBLIC touschek_lifetime ! subroutine called by programs referencing this module
PRIVATE integrand        ! local function called by touschek lifetime
PRIVATE exp_bessi0       ! local function called by integrand

CONTAINS

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!+
! Subroutine touschek_lifetime(mode, Tl, ring, orb)
!
! Calculates the touschek lifetime for a lattice.
! This calculation is based on Piwinski 1998 "The Touschek Effect In
! Strong Focusing Storage Rings".  This is the most general case, equation
! 42.
!
! This function assumes that the twiss parameters and closed orbit have
! been calculated, and that mode has been populated.
!
! Modules needed:
!   use touschek_mod
!
! Input:
!   mode            -- normal_modes_struct: beam properties
!     %pz_aperture    -- Real: momentum apperature of lattice
!     %sigE_E         -- Real: deltaE/E
!     %a%emittance    -- Real: horizontal emittance
!     %b%emittance    -- Real: vertical emittance, should be set to 
!                          a%emittance * coupling_constant
!     %sig_z          -- Real: bunch length
!   ring            -- lat_struct: Lattice for tracking
!     %param$n_part   -- Real: Number of particles in bunch, must be set 
!                          prior to calling touschek_lifetime
!   orb             -- Coord_struct array: closed orbit
! Output:
!   Tl              -- Real: Touschek lifetime in seconds
!
!-

SUBROUTINE touschek_lifetime(mode, Tl, ring, orb)

IMPLICIT NONE

TYPE(normal_modes_struct), INTENT(IN) :: mode
REAL(rp), INTENT(OUT) :: Tl
TYPE(lat_struct), INTENT(IN), target :: ring
TYPE(coord_struct), INTENT(IN) :: orb(0:)

TYPE(ele_struct), pointer :: ele

REAL(rp) E_tot, gamma, KE, beta, g2, beta2
REAL(rp) sigma_p2, sigma_z
REAL(rp) sigma_x2, sigma_y2, sigma_x_beta2, sigma_y_beta2
REAL(rp) NB,alpha_a,alpha_b,beta_a2,beta_b2
REAL(rp) Dx,Dy,Dxp,Dyp,Dxt,Dyt
REAL(rp) TL_inv,sum_Tl_inv,norm_sum_Tl_inv
REAL(rp) k_m, pi_2, integral, F
REAL(rp) emit_x, emit_y
REAL(rp) sigma_h2
REAL(rp) sigma_x_t2

INTEGER i

!

pi_2 = pi/2

sigma_p2 = mode%sigE_E**2
sigma_z = mode%sig_z
emit_x = mode%a%emittance
! this should be set to something sane like C*emit_x, by calling program
emit_y = mode%b%emittance 

E_tot = ring%ele(0)%value(E_TOT$)
CALL convert_total_energy_to(E_tot, ring%param%particle, gamma, KE, beta)
g2 = gamma**2
beta2 = beta**2

NB = ring%param%n_part

touschek_com%tau_m = beta2 * mode%pz_aperture**2

sum_Tl_inv = 0.0

DO i=1,ring%n_ele_track
  ele => ring%ele(i)

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a2 = ele%a%beta**2
  beta_b2 = ele%b%beta**2
  Dx = ele%a%eta
  Dxp = ele%a%etap
  Dy = ele%b%eta
  Dyp = ele%b%etap
  Dxt = alpha_a*Dx + SQRT(beta_a2)*Dxp
  Dyt = alpha_b*Dy + SQRT(beta_b2)*Dyp
  sigma_x_beta2 = SQRT(beta_a2) * emit_x
  sigma_y_beta2 = SQRT(beta_b2) * emit_y
  sigma_x2 = sigma_x_beta2 + (Dx**2)*sigma_p2
  sigma_y2 = sigma_y_beta2 + (Dy**2)*sigma_p2
  sigma_x_t2 = sigma_x2 + sigma_p2*(Dxt**2)

  sigma_h2 = 1/( 1/sigma_p2 &
    + ( (Dx**2)+(Dxt**2))/sigma_x_beta2 + ((Dy**2)+(Dyt**2))/sigma_y_beta2 )

  touschek_com%B1 = &
    beta_a2/2.0/beta2/g2/sigma_x_beta2*(1.0 - sigma_h2*(Dxt**2)/sigma_x_beta2) &
    + beta_b2/2.0/beta2/g2/sigma_y_beta2*(1.0 - sigma_h2*(Dyt**2)/sigma_y_beta2) 


  touschek_com%B2 = &
    SQRT( touschek_com%B1**2 - beta_a2*beta_b2*sigma_h2/(beta2**2)/(g2**2) &
    / (sigma_x_beta2**2)/(sigma_y_beta2**2)/(sigma_p2) &
    * (sigma_x2*sigma_y2 - (sigma_p2**2)*(Dx**2)*(Dy**2)) )

  k_m = ATAN(SQRT(touschek_com%tau_m))
  integral = qromb(integrand,k_m,pi_2)

  F = 2.0*SQRT( pi*(touschek_com%B1**2 - touschek_com%B2**2) ) * &
          touschek_com%tau_m * integral

  Tl_inv = (r_e**2)*(c_light)*(NB)/(8*pi)/(g2)/(sigma_z)/touschek_com%tau_m &
    / SQRT( sigma_x2*sigma_y2 - (sigma_p2**2)*(Dx**2)*(Dy**2) ) * F

  sum_Tl_inv = sum_Tl_inv + Tl_inv &
    * ( ring%ele(i)%value(l$)/2 + ring%ele(i+1)%value(l$)/2 )

ENDDO

norm_sum_Tl_inv = sum_Tl_inv / ring%param%total_length
Tl = 1 / norm_sum_Tl_inv

END SUBROUTINE touschek_lifetime

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------

FUNCTION integrand(k)
!this is the tan version of the integral in equation 42
!from Piwinski 1998
REAL(rp), DIMENSION(:), INTENT(IN) :: k
REAL(rp), DIMENSION(size(k)) :: integrand, t
REAL(rp) tau_m, B1, B2

INTEGER i

tau_m = touschek_com%tau_m
B1 = touschek_com%B1
B2 = touschek_com%B2

DO i=1, size(k)
  t(i) = TAN(k(i))**2
  integrand(i) = ( ((2.0*t(i)+1.0)**2)*(t(i)/tau_m/(1.0+t(i))-1.0)/t(i) &
    + t(i) - SQRT(t(i)*tau_m*(1.0+t(i))) &
    - (2.0 + 1.0/2.0/t(i))*LOG(t(i)/tau_m/(1.0+t(i))) )*SQRT(1.0+t(i)) &
    * exp_bessi0(t(i), B1, B2)

ENDDO

END FUNCTION integrand

!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
!--------------------------------------------------------------------------
! this is Numerial Recipes bessi0_s with exp(-B1*t) multiplied to the result

FUNCTION exp_bessi0(t, B1, B2)

USE nrtype; USE nrutil, ONLY : poly
IMPLICIT NONE
REAL(dp), INTENT(IN) :: t, B1, B2
REAL(dp) :: b2t
REAL(dp) :: exp_bessi0
REAL(DP), DIMENSION(7) :: p = (/1.0_dp,3.5156229_dp,&
  3.0899424_dp,1.2067492_dp,0.2659732_dp,0.360768e-1_dp,&
  0.45813e-2_dp/)
REAL(DP), DIMENSION(9) :: q = (/0.39894228_dp,0.1328592e-1_dp,&
  0.225319e-2_dp,-0.157565e-2_dp,0.916281e-2_dp,&
  -0.2057706e-1_dp,0.2635537e-1_dp,-0.1647633e-1_dp,&
  0.392377e-2_dp/)
b2t = B2 * t
IF (b2t < 3.75) THEN
  exp_bessi0=POLY(REAL((b2t/3.75_dp)**2,dp),p) * EXP(-B1*t)
ELSE
  exp_bessi0=(EXP((B2-B1)*t)/SQRT(b2t))*POLY(REAL(3.75_dp/(b2t),dp),q)
ENDIF

END FUNCTION exp_bessi0

END MODULE touschek_mod
