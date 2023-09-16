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

module touschek_mod

use twiss_and_track_mod
use fgsl
use, intrinsic :: iso_c_binding

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
PUBLIC integrand_base         ! Integral from equation 42 in Piwinski 1998
PUBLIC integrand_base_cov     ! Change of variables from t to exp(y) for integrand base
PUBLIC exp_bessi0             ! Exponential times Bessel I0, needed for integrand_base

! Data Types
PUBLIC momentum_aperture_struct

REAL(fgsl_double), PARAMETER :: eps7 = 1.0d-7
INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t

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
! Subroutine touschek_lifetime_ele_by_ele(mode, Tl, lat, momentum_aperture)
!
! Calculates the touschek lifetime for a lattice by calling touschek_rate1
! for each element the momentum_aperture array of momentum_aperture_structs.
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
! Input:
!   mode                     -- TYPE(normal_modes_struct), INTENT(INOUT): beam properties
!   lat                      -- TYPE(lat_struct), INTENT(IN): Lattice
!   momentum_aperture(:)     -- TYPE(momentum_aperture_struct), INTENT(IN): ele-by-ele unsymmatric apertures
!                       %s   -- Real(rp): ignored: in this subroutine, momentum_aperture is indexed by element index
!                       %pos -- Real(rp): positive momentum aperture
!                       %neg -- Real(rp): negative momentum aperture
!
! Output:
!   Tl                       -- Real(rp): Touschek lifetime in seconds
!-

subroutine touschek_lifetime_ele_by_ele(mode, Tl, lat, momentum_aperture)

implicit none

type(normal_modes_struct), intent(inout) :: mode
real(rp), intent(out) :: Tl
type(lat_struct), intent(in) :: lat
type(momentum_aperture_struct), intent(in) :: momentum_aperture(:)

integer i
real(rp) pos_rate, neg_rate
real(rp) sum_Tl_inv, rate, norm_sum_Tl_inv,Tl_inv

sum_Tl_inv = 0.0

do i=1,lat%n_ele_track
  if( lat%ele(i)%value(l$) .gt. 0.0001 ) then
    mode%pz_aperture = momentum_aperture(i)%pos
    call touschek_rate1(mode, rate, lat, ix=i)
    !call touschek_rate1_zap(mode, rate, lat, ix=i)
    pos_rate = rate / 2.0 !divide by two because Piwinski's formula assumes that two particles
                          !are lost per scattering event

    mode%pz_aperture = momentum_aperture(i)%neg
    call touschek_rate1(mode, rate, lat, ix=i)
    !call touschek_rate1_zap(mode, rate, lat, ix=i)
    neg_rate = rate / 2.0 

    Tl_inv = (pos_rate + neg_rate) / lat%param%n_part

    sum_Tl_inv = sum_Tl_inv + Tl_inv * lat%ele(i)%value(l$)
  endif
enddo

norm_sum_Tl_inv = sum_Tl_inv / lat%param%total_length
Tl = 1 / norm_sum_Tl_inv

end subroutine touschek_lifetime_ele_by_ele

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

subroutine touschek_lifetime_with_aperture(mode, Tl, lat, momentum_aperture)

implicit none

type(normal_modes_struct), intent(inout) :: mode
real(rp), intent(out) :: Tl
type(lat_struct), intent(in) :: lat
type(momentum_aperture_struct), intent(in) :: momentum_aperture(:)

integer i
integer nlocs
real(rp) pos_rate, neg_rate
real(rp) sum_Tl_inv, rate, norm_sum_Tl_inv,Tl_inv

nlocs = size(momentum_aperture)
sum_Tl_inv = 0.0

do i=2,nlocs
  mode%pz_aperture = momentum_aperture(i)%pos
  call touschek_rate1(mode, rate, lat, s=momentum_aperture(i)%s)
  !call touschek_rate1_zap(mode, rate, lat, s=momentum_aperture(i)%s)
  pos_rate = rate / 2.0 !divide by two because Piwinski's formula assumes that two particles
                        !are lost per scattering event

  mode%pz_aperture = momentum_aperture(i)%neg
  call touschek_rate1(mode, rate, lat, s=momentum_aperture(i)%s)
  !call touschek_rate1_zap(mode, rate, lat, s=momentum_aperture(i)%s)
  neg_rate = rate / 2.0 

  Tl_inv = (pos_rate + neg_rate) / lat%param%n_part

  sum_Tl_inv = sum_Tl_inv + Tl_inv * ( momentum_aperture(i)%s - momentum_aperture(i-1)%s )
enddo

norm_sum_Tl_inv = sum_Tl_inv / momentum_aperture(nlocs)%s
Tl = 1 / norm_sum_Tl_inv

end subroutine touschek_lifetime_with_aperture

!-
!-
! Touschek ZAP
!-
subroutine touschek_rate1_zap(mode, rate, lat, ix, s)

  implicit none

  type(normal_modes_struct) mode
  real(rp) rate
  type(lat_struct) lat
  real(rp), optional :: s
  integer, optional :: ix

  real(rp) E_tot, gamma, KE, beta, g2, beta2
  real(rp) sigma_p2, sigma_z
  real(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
  real(rp) sigma_xp
  real(rp) NB
  real(rp) Da,Db,Dap
  real(rp) integral
  real(rp) emit_a, emit_b
  real(rp) tau_m,  tau_param
  real(rp) Ha
  type(ele_struct) ele

  real(fgsl_double), target :: args(1:1)
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  if(present(s) .and. present(ix)) then
    write(*,*) "ERROR: ix and s cannot be specified at the same time."
    stop
  elseif(present(s)) then
    call twiss_and_track_at_s(lat, s, ele)
  elseif(present(ix)) then
    ele = lat%ele(ix)
  else
    write(*,*) "error: either ix or s must be specified when calling touschek_rate1."
    stop
  endif

  sigma_p2 = mode%sigE_E**2
  sigma_z = mode%sig_z
  E_tot = lat%ele(0)%value(E_TOT$)
  calL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)

  !Emittance is assumed to be normalized here.
  emit_a = mode%a%emittance
  emit_b = mode%b%emittance 
  tau_m = mode%pz_aperture

  NB = lat%param%n_part

  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap

  sigma_x_beta = sqrt(ele%a%beta * emit_a)
  sigma_y_beta = sqrt(ele%b%beta * emit_b)
  sigma_x = sqrt(sigma_x_beta**2 + (Da**2)*sigma_p2)
  sigma_y = sqrt(sigma_y_beta**2 + (Db**2)*sigma_p2)
  Ha = ele%a%gamma*Da*Da + 2.0d0*ele%a%alpha*Da*Dap + ele%a%beta*Dap*Dap
  sigma_xp = emit_a/sigma_x*sqrt(1.0d0 + Ha*sigma_p2/emit_a)

  tau_param = ( tau_m / gamma / sigma_xp ) ** 2

  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(integrand_zap, ptr)
  args = (/tau_param/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  integral = integration_result
  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)

  rate = classical_radius(lat%param%particle)**2 * c_light * NB**2 * integral / &
                (8.0_rp * pi * (gamma**3) * sigma_z * sigma_x * sigma_y * sigma_xp * tau_m * tau_m)

  if(rate .lt. 0.0_rp) rate = 0.0_rp

end subroutine touschek_rate1_zap

function integrand_zap(u,args) bind(c)

implicit none

real(c_double), value :: u
type(c_ptr), value :: args
real(c_double) :: integrand_zap
real(c_double), pointer :: local_args(:)
real(rp) tau_param
real(rp) oou

call c_f_pointer(args,local_args,[1])

tau_param = local_args(1)
oou = 1.0d0/u

integrand_zap = ( oou - 0.5d0*log(oou) - 1.0d0 ) * exp(-tau_param*oou)

!integrand_zap = MAX(integrand_base,0._rp)

END FUNCTION integrand_zap

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
!+
! Subroutine touschek_rate1(mode, rate, lat, ix, s)
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
!                 contain the normalized emittences.  If lattice geometry is
!                 open, the emittances are assumed to be
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
! Input:
!   mode            -- TYPE(normal_modes_struct), INTENT(IN): beam properties
!   lat             -- TYPE(lat_struct): Lattice
!   ix              -- INTEGER, OPTIONAL, INTENT(IN): element index (either s or ix must be specified)
!   s               -- REAL(rp), OPTIONAL, INTENT(IN): location in meters (either s or ix must be specified) 
! Output:
!   rate            -- REAL(rp), INTENT(OUT): Touschek rate, in units particle per second, assuming
!                                two particles per event.
!-

SUBROUTINE touschek_rate1(mode, rate, lat, ix, s)

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
  REAL(rp) Da,Db,Dap,Dbp,Dat,Dbt
  REAL(rp) pi_2, integral
  REAL(rp) emit_a, emit_b
  REAL(rp) sigma_h2
  REAL(rp) sigma_x_t2
  REAL(rp) tau_min, tau_max
  REAL(rp) tau_m
  REAL(rp) B1
  REAL(rp) B2, B2alt
  TYPE(ele_struct) ele

  REAL(fgsl_double), TARGET :: args(1:3)
  TYPE(fgsl_function) :: integrand_ready
  REAL(fgsl_double) :: integration_result
  REAL(fgsl_double) :: abserr
  TYPE(c_ptr) :: ptr
  TYPE(fgsl_integration_workspace) :: integ_wk
  INTEGER(fgsl_int) :: fgsl_status

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

  IF (lat%param%geometry == closed$) THEN
    sigma_p2 = mode%sigE_E**2
    sigma_z = mode%sig_z
    E_tot = lat%ele(0)%value(E_TOT$)
    CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)
    g2 = gamma**2
    beta2 = beta**2
    !Emittance is assumed to be normalized here.
    emit_a = mode%a%emittance
    emit_b = mode%b%emittance 
    tau_m = beta2 * mode%pz_aperture**2
  ELSEIF (lat%param%geometry == open$) THEN
    sigma_p2 = ele%z%sigma_p**2
    sigma_z = ele%z%sigma
    E_tot = ele%value(E_TOT$)
    CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)
    g2 = gamma**2
    beta2 = beta**2
    !Emittance is assumes to be UNNORMALIZED here
    emit_a = mode%a%emittance / gamma
    emit_b = mode%b%emittance / gamma
    tau_m = beta2 * mode%pz_aperture**2
  ELSE
    WRITE(*,*) "ERROR: geometry unknown. Halting."
    STOP
  ENDIF

  NB = lat%param%n_part

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a2 = ele%a%beta**2
  beta_b2 = ele%b%beta**2
  Da = ele%a%eta
  Dap = ele%a%etap
  Db = ele%b%eta
  Dbp = ele%b%etap
  Dat = alpha_a*Da + ele%a%beta*Dap
  Dbt = alpha_b*Db + ele%b%beta*Dbp
  sigma_x_beta2 = ele%a%beta * emit_a
  sigma_y_beta2 = ele%b%beta * emit_b
  sigma_x2 = sigma_x_beta2 + (Da**2)*sigma_p2
  sigma_y2 = sigma_y_beta2 + (Db**2)*sigma_p2
  sigma_x_t2 = sigma_x2 + sigma_p2*(Dat**2)

  sigma_h2 = 1._rp/( 1._rp/sigma_p2 &
    + ( (Da**2)+(Dat**2))/sigma_x_beta2 + ((Db**2)+(Dbt**2))/sigma_y_beta2 )

  B1 = beta_a2/2.0_rp/beta2/g2/sigma_x_beta2*(1.0_rp - sigma_h2*(Dat**2)/sigma_x_beta2) &
       + beta_b2/2.0_rp/beta2/g2/sigma_y_beta2*(1.0_rp - sigma_h2*(Dbt**2)/sigma_y_beta2) 

  ! B2 = SQRT( (0.25_rp)/(beta2**2)/(g2**2) &
  !                   *( beta_a2/sigma_x_beta2*(1._rp-sigma_h2*(Dat**2)/sigma_x_beta2) - &
  !                      beta_b2/sigma_y_beta2*(1._rp-sigma_h2*(Dbt**2)/sigma_y_beta2) )**2 &
  !                   + ( (sigma_h2**2)*beta_a2*beta_b2*(Dat**2)*(Dbt**2)/(beta2**2)/(g2**2) &
  !                     /(sigma_x_beta2**2)/(sigma_y_beta2**2) ) )

  B2 = B1*B1 - beta_a2*beta_b2*sigma_h2/beta2/beta2/g2/g2/sigma_x_beta2/sigma_x_beta2/sigma_y_beta2/sigma_y_beta2/sigma_p2 * &
                  (sigma_x2*sigma_y2-sigma_p2*sigma_p2*Da*Da*Db*Db)
  B2 = sqrt(B2)

  tau_min = tau_m
  tau_max = 1.0_rp

  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(integrand_base_cov, ptr)
  args = (/tau_m, B1, B2/)
  fgsl_status = fgsl_integration_qag(integrand_ready, LOG(tau_min), LOG(tau_max), eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  integral = integration_result
  CALL fgsl_integration_workspace_free(integ_wk)
  CALL fgsl_function_free(integrand_ready)

  !integral = qtrap(integrand_base_cov, LOG(tau_min), LOG(tau_max))

  rate = (classical_radius(lat%param%particle)**2)*c_light*(NB**2)/8.0_rp/SQRT(pi)/ &
                           (g2*g2)/beta2/sigma_z/SQRT(sigma_p2/sigma_h2)/emit_a/emit_b* integral

  !Simplified (a small bit) from thesis
  !rate = (r_e**2)*c_light*(NB**2)/8.0_rp/SQRT(pi)/ &
  !       (g2*g2)/sigma_z/SQRT(sigma_p2/sigma_h2)/emit_a/emit_b* &
  !       integral


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

FUNCTION integrand_base(t,args)

!this is the the integral in equation 42 from Piwinski 1998

implicit none

REAL(rp) :: t
REAL(rp) :: integrand_base
REAL(rp) args(:)
REAL(rp) tm, B1, B2

tm = args(1)
B1 = args(2)
B2 = args(3)

! Classic Piwinski
integrand_base = ( (2._rp+1._rp/t)**2*(t/tm/(1.+t)-1._rp)+ 1._rp - SQRT(1._rp+t)/SQRT(t/tm) &
                        - 0.5_rp/t*(4._rp+1._rp/t)*LOG(t/tm/(1._rp+t)) ) * SQRT(t) / SQRT(1._rp+t) &
                      * exp_bessi0(t, B1, B2)

! Simplified from thesis
! Also set tm=deltaE,max**2
!integrand_base =   ( (1.0d0/t+2.0d0)**2*(t/tm/(1.0d0+t) - 1.0d0) + 1.0d0 - SQRT(tm*(1.0d0+t)/t) &
!                     + (1.0d0/t/t+4.0d0/t)*LOG(SQRT(tm*(1.0d0+t)/t)) ) * SQRT(t) / SQRT(1.0d0+t) &
!                  * exp_bessi0(t, B1, B2)
                
integrand_base = MAX(integrand_base,0._rp)

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

FUNCTION integrand_base_cov(y,params) BIND(c)

implicit none

REAL(c_double), VALUE :: y
TYPE(c_ptr), VALUE :: params
REAL(c_double) :: integrand_base_cov
REAL(c_double), POINTER :: args(:)
REAL(c_double) :: tau_m, B1, B2

CALL c_f_pointer(params,args,[3])

integrand_base_cov=integrand_base(EXP(y),args)*EXP(y)

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

use super_recipes_mod, only: super_poly

implicit none

REAL(rp), INTENT(IN) :: t, B1, B2
REAL(rp) :: b2t
REAL(rp) :: exp_bessi0

REAL(rp), parameter :: p(7) = [1.0_rp, 3.5156229_rp, 3.0899424_rp, 1.2067492_rp, 0.2659732_rp, &
                                                                          0.360768e-1_rp, 0.45813e-2_rp]
REAL(rp), parameter :: q(9) = [0.39894228_rp, 0.1328592e-1_rp, 0.225319e-2_rp, -0.157565e-2_rp, &
                   0.916281e-2_rp, -0.2057706e-1_rp, 0.2635537e-1_rp, -0.1647633e-1_rp, 0.392377e-2_rp]
b2t = B2 * t

IF (b2t < 3.75) THEN
  exp_bessi0=super_POLY(REAL((b2t/3.75_rp)**2,rp),p) * EXP(-B1*t)
ELSE
  exp_bessi0=(EXP((B2-B1)*t)/SQRT(b2t))*super_POLY(REAL(3.75_rp/(b2t),rp),q)
ENDIF

END FUNCTION exp_bessi0

END MODULE touschek_mod
