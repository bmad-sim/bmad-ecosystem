!+
! Module touschek_mod
!
! This module calculates the Touschek scattering rate found in Eq. 31 of
! Piwinski's paper "The Touschek Effect in Strong Focusing Storage Rings".
! This rate is useful for calculating beam lifetimes (in the case of storage
! rings), and background radiation from particle loss (in the case of linacs).
!
! Piwinski's Touschek formula is very general in that it uses the relativistic
! Molller scattering cross-section, and so is useful in cases where particles
! in the COM frame are relativistic, and also in that it takes into account
! dispersion, which makes the formula useful in cases where the beam pipe
! aperture is significant to the momentum aperture.
!-

MODULE touschek_mod

USE bmad

! This common block is needed in order to pass parameters to
! integrand.  We use NR's qromb function, which requires that the
! function it integrates be one variable only.
TYPE touschek_common_struct
  REAL(rp) B1    !argument to the exponential in the integrand.  Piwinski Eqn. 33
  REAL(rp) B2    !argument to the Bessel function in the integrand.  Piwinski Eqn. 34
  REAL(rp) tau_m !momentum aperture
END TYPE touschek_common_struct
TYPE(touschek_common_struct), PRIVATE :: touschek_com

!This structure contains the positive and negative momentum aperture for locations s.
TYPE momentum_aperture_struct
  REAL(rp) s
  REAL(rp) pos
  REAL(rp) neg
END TYPE momentum_aperture_struct

! Functions and Subroutines
PUBLIC touschek_lifetime       ! Calculates Touschek lifetime by averaging over elements
PUBLIC touschek_lifetime_with_aperture     ! Calculates Touschek lifetime by averaging over set of locations
PUBLIC touschek_rate1          ! Returns touschek rate for either a location s or an element ix
PRIVATE integrand_base         ! Integral from equation 42 in Piwinski 1998
PRIVATE integrand_base_cov     ! Change of variables from t to exp(y) for integrand base
PRIVATE exp_bessi0             ! Exponential times Bessel I0, needed for integrand_base

! Data Types
PUBLIC touschek_common_struct
PUBLIC momentum_aperture_struct

CONTAINS

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine touschek_lifetime(mode, Tl, lat)
!
! Calculates the touschek lifetime for a lattice by calling touschek_rate1
! for each element.
! The loss rate at each element is averaged over one turn to obtain the lifetime.
!
! This function assumes that the twiss parameters and closed orbit have
! been calculated, and that mode has been populated.
!
! This subroutine assumes a fixed momentum aperture.  The loss rate at each element
! uses the same momentum aperture, mode%pz_aperture.
!
! A common way to call this function is to first populate mode using
! radiation integrals.  If an ideal lattice is used, the vertical
! emittance must also be set to a reasonable value.  If the vertical
! emittance is due only to quantum excitation, then it will likely be
! several orders of magnitude smaller than any real physical situation, in which
! case the integral in this function will have problems converging.
!
! In addition to setting mode, also set lat%param%n_part to the number of particles
! per bunch.
!
! Modules needed:
!   use touschek_mod
!
! Input:
!   mode             -- TYPE(normal_modes_struct), INTENT(INOUT): beam properties
!       %pz_aperture -- Real(rp): momentum aperture
!   lat              -- TYPE(lat_struct): Accelerator Lattice
!      %param%n_part -- Real(rp): number particles per bunch
!
! Output:
!   Tl               -- Real(rp): Touschek lifetime in seconds
!
!-

SUBROUTINE touschek_lifetime(mode, Tl, lat)

IMPLICIT NONE

TYPE(normal_modes_struct), INTENT(INOUT) :: mode
REAL(rp), INTENT(OUT) :: Tl
TYPE(lat_struct) :: lat

INTEGER i
REAL(rp) NB
REAL(rp) sum_Tl_inv, rate, norm_sum_Tl_inv,Tl_inv

NB = lat%param%n_part

sum_Tl_inv = 0.0

DO i=1,lat%n_ele_track
  CALL touschek_rate1(mode, rate, lat, ix=i)
  Tl_inv = rate / NB

  sum_Tl_inv = sum_Tl_inv + Tl_inv &
    * ( lat%ele(i)%value(l$)/2.0 + lat%ele(i+1)%value(l$)/2.0)
ENDDO

norm_sum_Tl_inv = sum_Tl_inv / lat%param%total_length
Tl = 1 / norm_sum_Tl_inv

END SUBROUTINE touschek_lifetime

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine touschek_lifetime_detailed_aperture(mode, Tl, lat, momentum_aperture)
!
! Calculates the touschek lifetime for a lattice by calling touschek_rate1
! for each s-coordinate in the momentum_aperture array of momentum_aperture_structs.
! This calculation is based on Piwinski 1998 "The Touschek Effect In
! Strong Focusing Storage Rings".  This is the most general case, equation 31.
! 42.
!
! A common way to call this function is to first populate mode using
! radiation integrals.  If an ideal lattice is used, the vertical
! emittance must also be set to a reasonable value.  If the vertical
! emittance is due only to quantum excitation, then it will likely be
! several orders of magnitude smaller than any real physical situation, in which
! case the integral in this function will have problems converging.
!
! In addition to setting mode, also set lat%param%n_part to the number of particles
! per bunch.
!
! This function assumes that the twiss parameters 
! been calculated, and that mode has been populated with emittance and bunch length.
!
! This function assumes that momentum_aperture(0)%s==0 and momentum_aperture(last)%s==lat%param%total_length.
!
! Modules needed:
!   use touschek_mod
!
! Input:
!   mode                     -- TYPE(normal_modes_struct), INTENT(INOUT): beam properties
!   lat                      -- TYPE(lat_struct), INTENT(IN): Lattice
!   momentum_aperture(:)     -- TYPE(momentum_aperture_struct), INTENT(IN): loc-by-loc unsymmatric apertures
!                       %s   -- Real(rp): location in lattice
!                       %pos -- Real(rp): positive momentum aperture
!                       %neg -- Real(rp): negative momentum aperture
!
! Output:
!   Tl                       -- Real(rp): Touschek lifetime in seconds
!-

SUBROUTINE touschek_lifetime_with_aperture(mode, Tl, lat, momentum_aperture)

IMPLICIT NONE

TYPE(normal_modes_struct), INTENT(INOUT) :: mode
REAL(rp), INTENT(OUT) :: Tl
TYPE(lat_struct), INTENT(IN) :: lat
TYPE(momentum_aperture_struct), INTENT(IN) :: momentum_aperture(:)

INTEGER i
INTEGER nlocs
REAL(rp) pos_rate, neg_rate
REAL(rp) sum_Tl_inv, rate, norm_sum_Tl_inv,Tl_inv

nlocs = SIZE(momentum_aperture)
sum_Tl_inv = 0.0

DO i=1,nlocs-1  !skip last entry, since that is either the end of a linac,
                !or is taken into account as the start of the storage ring
  mode%pz_aperture = momentum_aperture(i)%pos
  CALL touschek_rate1(mode, rate, lat, s=momentum_aperture(i)%s)
  pos_rate = rate / 2.0 !divide by two because Piwinski's formula assumes that two particles
                        !are lost per scattering event

  mode%pz_aperture = momentum_aperture(i)%neg
  CALL touschek_rate1(mode, rate, lat, s=momentum_aperture(i)%s)
  neg_rate = rate / 2.0 

  Tl_inv = (pos_rate + neg_rate) / lat%param%n_part

  sum_Tl_inv = sum_Tl_inv + Tl_inv &
    * ( momentum_aperture(i+1)%s - momentum_aperture(i)%s )
ENDDO

norm_sum_Tl_inv = sum_Tl_inv / momentum_aperture(nlocs)%s
Tl = 1 / norm_sum_Tl_inv

END SUBROUTINE touschek_lifetime_with_aperture

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine touschek_rate1(mode, rate, lat, s=s)
! Subroutine touschek_rate1(mode, rate, lat, ix=ix)
!
! Calculates the touschek rate at the location specified by s or ix
! This calculation is based on Piwinski 1998 "The Touschek Effect In
! Strong Focusing Storage Rings".  This is the most general case, equation
! 31.
!
! This function uses twiss_and_track_at_s to determine the Twiss parameters
! at the location s or element index ix.
!
! A common way to call this function is to first populate mode using
! radiation integrals.  If an ideal lattice is used, the vertical
! emittance must also be set to a reasonable value.  If the vertical
! emittance is due only to quantum excitation, then it will likely be
! several orders of magnitude smaller than any real physical situation, in which
! case the integral in this function will have problems converging.
! Additionally, mode%pz_aperture needs to be set to the momentum aperture.
!
! In addition to setting mode, also set lat%param%n_part to the number of particles
! per bunch.
!
! IMPORTANT NOTE: If the lattice type is a circular lattice, then 
!                 mode%a%emittance and mode%b%emittance are assumed to
!                 contain the normalized emittences.  If lattice type is
!                 linear lattice, then the emittances are assumed to be
!                 unnormalized.
!
! IMPORTANT NOTE: The output of this subroutine is the loss rate assuming
!                 that two particles are lost per collision, one with too
!                 much energy, and one with too little energy.  This agrees
!                 with Piwinski's original derivation, which assumes that the 
!                 positive energy aperture is equal in magnitude to the
!                 negative energy aperture.  If you are studying an
!                 accelerator with a non-symmetric energy aperture, then 
!                 this subroutine should be called twice, once with the positive
!                 aperture, and once with the negative aperture, and rate from 
!                 each call should be halved and summed.
!
! Modules needed:
!   use touschek_mod
!
! Input:
!   mode            -- TYPE(normal_modes_struct), INTENT(IN): beam properties
!   lat             -- TYPE(lat_struct): Lattice
!   s=s             -- REAL(rp), OPTIONAL, INTENT(IN): location in meters (either s or ix must be specified) 
!   ix=ix           -- INTEGER, OPTIONAL, INTENT(IN): element index (either s or ix must be specified)
! Output:
!   rate            -- REAL(rp), INTENT(OUT): Touschek rate, in units particle per second, assuming
!                                two particles per event.
!-

SUBROUTINE touschek_rate1(mode, rate, lat, ix, s)

  use nr, only: qtrap

  IMPLICIT NONE

  TYPE(normal_modes_struct), INTENT(IN) :: mode
  REAL(rp), INTENT(OUT) :: rate
  TYPE(lat_struct), INTENT(IN) :: lat
  REAL(rp), INTENT(IN), OPTIONAL :: s
  INTEGER, INTENT(IN), OPTIONAL :: ix

  REAL(rp) E_tot, gamma, KE, beta, g2, beta2
  REAL(rp) sigma_p2, sigma_z
  REAL(rp) sigma_x2, sigma_y2, sigma_x_beta2, sigma_y_beta2
  REAL(rp) NB,alpha_a,alpha_b,beta_a2,beta_b2
  REAL(rp) Dx,Dy,Dxp,Dyp,Dxt,Dyt
  REAL(rp) pi_2, integral
  REAL(rp) emit_x, emit_y
  REAL(rp) sigma_h2
  REAL(rp) sigma_x_t2
  REAL(rp) tau_min, tau_max
  TYPE(ele_struct) ele

  pi_2 = pi/2._rp

  IF(PRESENT(s) .and. PRESENT(IX)) THEN
    WRITE(*,*) "ERROR: ix and s cannot be specified at the same time."
    STOP
  ELSEIF(PRESENT(s)) THEN
    CALL twiss_and_track_at_s(lat, s, ele)
  ELSEIF(PRESENT(ix)) THEN
    ele = lat%ele(ix)
  ELSE
    WRITE(*,*) "ERROR: Either ix or s must be specified when calling touschek_rate1."
    STOP
  ENDIF

  IF (lat%param%lattice_type == circular_lattice$) THEN
    sigma_p2 = mode%sigE_E**2
    sigma_z = mode%sig_z
    E_tot = lat%ele(0)%value(E_TOT$)
    CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)
    g2 = gamma**2
    beta2 = beta**2
    !Emittance is assumed to be normalized here.
    emit_x = mode%a%emittance
    emit_y = mode%b%emittance 
    touschek_com%tau_m = beta2 * mode%pz_aperture**2
  ELSEIF (lat%param%lattice_type == linear_lattice$) THEN
    sigma_p2 = ele%z%sigma_p**2
    sigma_z = ele%z%sigma
    E_tot = ele%value(E_TOT$)
    CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)
    g2 = gamma**2
    beta2 = beta**2
    !Emittance is assumes to be UNNORMALIZED here
    emit_x = mode%a%emittance / gamma
    emit_y = mode%b%emittance / gamma
    touschek_com%tau_m = beta2 * mode%pz_aperture**2
  ELSE
    WRITE(*,*) "ERROR: lattice_type unknown. Halting."
    STOP
  ENDIF

  NB = lat%param%n_part

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a2 = ele%a%beta**2
  beta_b2 = ele%b%beta**2
  Dx = ele%a%eta
  Dxp = ele%a%etap
  Dy = ele%b%eta
  Dyp = ele%b%etap
  Dxt = alpha_a*Dx + ele%a%beta*Dxp
  Dyt = alpha_b*Dy + ele%b%beta*Dyp
  sigma_x_beta2 = ele%a%beta * emit_x
  sigma_y_beta2 = ele%b%beta * emit_y
  sigma_x2 = sigma_x_beta2 + (Dx**2)*sigma_p2
  sigma_y2 = sigma_y_beta2 + (Dy**2)*sigma_p2
  sigma_x_t2 = sigma_x2 + sigma_p2*(Dxt**2)

  sigma_h2 = 1._rp/( 1._rp/sigma_p2 &
    + ( (Dx**2)+(Dxt**2))/sigma_x_beta2 + ((Dy**2)+(Dyt**2))/sigma_y_beta2 )

  touschek_com%B1 = &
    beta_a2/2.0_rp/beta2/g2/sigma_x_beta2*(1.0_rp - sigma_h2*(Dxt**2)/sigma_x_beta2) &
    + beta_b2/2.0_rp/beta2/g2/sigma_y_beta2*(1.0_rp - sigma_h2*(Dyt**2)/sigma_y_beta2) 

  touschek_com%B2 = SQRT( (0.25_rp)/(beta2**2)/(g2**2) &
                    *( beta_a2/sigma_x_beta2*(1._rp-sigma_h2*(Dxt**2)/sigma_x_beta2) - &
                      beta_b2/sigma_y_beta2*(1._rp-sigma_h2*(Dyt**2)/sigma_y_beta2) )**2 &
                    + ( (sigma_h2**2)*beta_a2*beta_b2*(Dxt**2)*(Dyt**2)/(beta2**2)/(g2**2) &
                        /(sigma_x_beta2**2)/(sigma_x_beta2**2) ) )

  tau_min = touschek_com%tau_m
  tau_max = 1.0_rp

  integral = qtrap(integrand_base_cov, LOG(tau_min), LOG(tau_max))

  rate = (r_e**2)*c_light*(NB**2)/8.0_rp/SQRT(pi)/ &
         (g2*g2)/beta2/sigma_z/SQRT(sigma_p2/sigma_h2)/emit_x/emit_y* &
         integral

  IF(rate .lt. 0.0_rp) rate = 0.0_rp

END SUBROUTINE touschek_rate1

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Function integrand_base(t)
!
! This vectorized private function is the integrand in equation 31 of Piwinski's paper.
!
! This intetegrand has a sharp exponential decay, and so a change of variables from t to y where t=exp(y)
! is applied.  This COV makes the integrand more evenly distributed over the domain of integration,
! which makes it easier for qtrap to integrate.
!
! The change of variables is done using integrand_base_cov, which is then integrated
! using qtrap.
!
! Input:
!   t(:)           -- REAL(rp), INTENT(IN): Array of reals over which to evaluate the integrand.
!
! Output:
!   <return value> -- REAL(rp): Array of reals containing values of integrand at t(:).
!-

FUNCTION integrand_base(t)

!this is the the integral in equation 42 from Piwinski 1998

implicit none

REAL(rp), DIMENSION(:), INTENT(IN) :: t
REAL(rp), DIMENSION(size(t)) :: integrand_base
REAL(rp) tm, B1, B2

INTEGER i

tm = touschek_com%tau_m
B1 = touschek_com%B1
B2 = touschek_com%B2

DO i=1, size(t)
  integrand_base(i) = ( (2._rp+1._rp/t(i))**2*(t(i)/tm/(1.+t(i))-1._rp)+ 1._rp - SQRT(1._rp+t(i))/SQRT(t(i)/tm) &
                        - 0.5_rp/t(i)*(4._rp+1._rp/t(i))*LOG(t(i)/tm/(1._rp+t(i))) ) * SQRT(t(i)) / SQRT(1._rp+t(i)) &
                      * exp_bessi0(t(i), B1, B2)
  integrand_base(i) = MAX(integrand_base(i),0._rp)
ENDDO

END FUNCTION integrand_base

!+
! Function integrand_base_cov
!
! This vectorized private function performs a change of variables on integrand_base.
! The cov is from t to y where t=exp(y).  Accordingly, qtrap needs to be called with
! lower bound log(tlower) and upper bound log(tupper).
!
! Input:
!   y(:)           -- Real(rp), INTENT(IN): Array of reals over which to evaluate the integrand.
!
! Output:
!   <return value> -- Real(rp): Array of reals containing the values of the integrand evaluated at t(:)
!-

FUNCTION integrand_base_cov(y)

implicit none

REAL(rp), DIMENSION(:), INTENT(IN) :: y
REAL(rp), DIMENSION(size(y)) :: integrand_base_cov

integrand_base_cov=integrand_base(EXP(y))*EXP(y)

END FUNCTION integrand_base_cov

!+
! Function exp_bessi0(t, B1, B2)
!
! This is essentially the Numercal Recipes bessi0 function multiplied by exp(-B1*t).
!
! This overcomes an issue where exp(B2*t) may be huge and exp(-B1*t) may be small.
! Evaluating exp(B2*t) may result in overflow, but exp((B2-B1)*t) has a moderate value.
! Simplifying the algebra of B2-B1 suggests that is should always have a moderate magnitude.
!
! Input:
!   t              -- Real(rp), INTENT(IN): Scalar agrument to evaluate function at.
!   B1             -- Real(rp), INTENT(IN): Scalar value.  Eq. 33 from Piwinski's paper.
!   B2             -- Real(rp), INTENT(IN): Scalar value.  Eq. 34 from Piwinski's paper.
!
! Output:
!   <return value> -- Real(rp): Scalar return value.
!-
FUNCTION exp_bessi0(t, B1, B2)

USE nrtype; USE nrutil, ONLY : poly

IMPLICIT NONE

REAL(rp), INTENT(IN) :: t, B1, B2
REAL(rp) :: b2t
REAL(rp) :: exp_bessi0

REAL(rp), DIMENSION(7) :: p = [1.0_dp, 3.5156229_dp, 3.0899424_dp, 1.2067492_dp, 0.2659732_dp, &
                              0.360768e-1_dp, 0.45813e-2_dp]
REAL(rp), DIMENSION(9) :: q = [0.39894228_dp, 0.1328592e-1_dp, 0.225319e-2_dp, -0.157565e-2_dp, &
                   0.916281e-2_dp, -0.2057706e-1_dp, 0.2635537e-1_dp, -0.1647633e-1_dp, 0.392377e-2_dp]
b2t = B2 * t

IF (b2t < 3.75) THEN
  exp_bessi0=POLY(REAL((b2t/3.75_dp)**2,dp),p) * EXP(-B1*t)
ELSE
  exp_bessi0=(EXP((B2-B1)*t)/SQRT(b2t))*POLY(REAL(3.75_dp/(b2t),dp),q)
ENDIF

END FUNCTION exp_bessi0

END MODULE touschek_mod
