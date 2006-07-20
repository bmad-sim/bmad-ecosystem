#include "CESR_platform.inc"

!+
! Module ibs_mod
!
! Module for interbeam scattering calculations.
!
! Module developer: Michael Ehrlichman
!-

MODULE ibs_mod

USE bmad_struct
USE bmad_interface
USE nr

TYPE ibs_struct
  REAL(rp) inv_Tx
  REAL(rp) inv_Ty
  REAL(rp) inv_Tz
END TYPE

PUBLIC ibs_equilibrium
PUBLIC ibs_rates
PRIVATE cimp
PRIVATE bane
PRIVATE bjmt

CONTAINS

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  Subroutine ibs_equilibrium(ring,inmode,ibsmode,formula,coupling)
!
!  Computes equilibrium beam sizes taking into account
!  radiation damping, IBS growth rates, and coupling.
!
!  Modules needed:
!    use ibs_mod
!
!  Input:
!    ring             -- ring_struct: lattice for tracking
!      %param$n_part  -- Real: number of particles in bunch
!    inmode           -- modes_struct: natural beam parameters 
!    formula          -- character(4): IBS formulation to use (see ibs_rates)
!    coupling         -- real: horizontal to vertical emittanc coupling
!
!  Output:
!    ibsmode          -- modes_struct: beam parameters after IBS effects
!
!  See ibs_rates subroutine for available IBS rate formulas.
!-

SUBROUTINE ibs_equilibrium(ring,inmode,ibsmode,formula,coupling)

IMPLICIT NONE

TYPE(ring_struct), INTENT(IN), target :: ring
TYPE(modes_struct), INTENT(IN) :: inmode
TYPE(modes_struct), INTENT(OUT) :: ibsmode
CHARACTER(*), INTENT(IN) :: formula
REAL(rp), INTENT(IN), OPTIONAL :: coupling
TYPE(ele_struct), pointer :: ele
TYPE(ibs_struct) rates
TYPE(modes_struct) :: naturalmode

REAL(rp) time_for_one_turn
REAL(rp) tau_x, tau_y, tau_z
REAL(rp) Tx, Ty, Tz
REAL(rp) emit_x, emit_y, emit_z
REAL(rp) emit_x0, emit_y0, emit_z0
REAL(rp) dxdt, dydt, dzdt, dT
REAL(rp) threshold
REAL(rp) sigE_E0, sigma_z0, L_ratio
REAL(rp) ka, ka_one, ka_small
LOGICAL converged
INTEGER counter

! natural mode here means emittances before IBS effects

naturalmode = inmode
IF(present(coupling)) THEN
  ka = coupling
  naturalmode%b%emittance = coupling * inmode%a%emittance
ELSE
  ka = inmode%b%emittance / inmode%a%emittance
ENDIF

! compute the SR damping times
time_for_one_turn = ring%param%total_length / c_light
tau_x = time_for_one_turn / naturalmode%a%alpha_damp
tau_y = time_for_one_turn / naturalmode%b%alpha_damp
tau_z = time_for_one_turn / naturalmode%z%alpha_damp 

sigma_z0 = naturalmode%sig_z
sigE_E0 = naturalmode%sigE_E
emit_x0 = naturalmode%a%emittance
emit_y0 = naturalmode%b%emittance
emit_z0 = sigma_z0 * sigE_E0
L_ratio = sigma_z0 / sigE_E0
threshold = .00001   ! fractional changes in emittance smaller than this
                     ! indicate convergence
converged = .false.
dT = tau_x / 40.0       ! Time to advance per iteration
ka_one = 1. / (1.+ka)   ! Used in determining de_dt
ka_small = ka / (1.+ka) ! Used in determining de_dt

counter = 0
ibsmode = naturalmode
DO WHILE(.not.converged)
  CALL ibs_rates(ring,ibsmode,rates,formula)
  counter = counter + 1
  Tx = 1.0/rates%inv_Tx
  Ty = 1.0/rates%inv_Ty
  Tz = 1.0/rates%inv_Tz
  emit_x = ibsmode%a%emittance
  emit_y = ibsmode%b%emittance
  emit_z = ibsmode%sigE_E*ibsmode%sig_z

  ! Compute change in emittance per time for x,y,z dimensions, taking
  ! into account radiation damping, IBS blow-up, and x-y coupling
  dxdt = -(emit_x-ka_one*emit_x0)/tau_x + emit_x/Tx
  dydt = -(emit_y-ka_small*emit_x)/tau_y + emit_y/Ty
  dzdt = -(emit_z - emit_z0)/tau_z + emit_z/Tz

  IF( (dxdt*dT)/emit_x .lt. threshold ) THEN
    IF( (dydt*dT)/emit_y .lt. threshold ) THEN
      IF( (dzdt*dT)/emit_z .lt. threshold ) THEN
        converged = .true.
      ENDIF
    ENDIF
  ENDIF

  IF(.not.converged) THEN
    ibsmode%a%emittance = emit_x + dxdt*dT
    ibsmode%b%emittance = emit_y + dydt*dT
    emit_z = emit_z + dzdt*dT
    ibsmode%sig_z =  SQRT(emit_z*L_ratio)
    ibsmode%sigE_E = emit_z / ibsmode%sig_z
  ENDIF
ENDDO

END SUBROUTINE ibs_equilibrium

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  Subroutine ibs_rates(ring, mode, rates, formula)
!
!  Calculates IBS risetimes for given ring and mode.
!  This is basically a front-end for the various formulas 
!  available in this module of calculating IBS rates.
!
!  Available IBS formulas:
!    cimp - Modified Piwinski
!    bane - Bane approximation of Bjorken-Mtingwa formulation
!    bjmt - General Bjorken-Mtingwa formulation for bunched beams (slow)
!
!  Input:
!    ring             -- ring_struct: lattice for tracking
!      %param$n_part  -- Real: number of particles in bunch
!    mode             -- modes_struct: beam parameters 
!    formula          -- character(4): IBS formulation to use
!
!  Output:
!    rates          -- ibs_struct: ibs rates in x,y,z 
!-

SUBROUTINE ibs_rates(ring, mode, rates, formula)

IMPLICIT NONE

TYPE(modes_struct), INTENT(IN) :: mode
TYPE(ring_struct), INTENT(IN), target :: ring
TYPE(ibs_struct), INTENT(OUT) :: rates
CHARACTER(*), INTENT(IN) ::  formula

IF(formula == 'cimp') THEN
  CALL cimp(ring, mode, rates)
ELSEIF(formula == 'bjmt') THEN
  CALL bjmt(ring, mode, rates)
ELSEIF(formula == 'bane') THEN
  CALL bane(ring, mode, rates)
ELSE
  WRITE(*,*) "Invalid IBS formula selected ... returning zero"
  rates%inv_Tx = 0.0
  rates%inv_Ty = 0.0
  rates%inv_Tz = 0.0
ENDIF

END SUBROUTINE ibs_rates

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  subroutine bjmt(ring, mode, rates)
!
!  This is a private subroutine.  To access this subroutine, call
!  ibs_rates.
!
!  This is an implementation of equations 1-9 from "Intrabeam
!  scattering formulas for high energy beams" Kubo, Mtingwa, and Wolski.
!  It is the most general form (for bunched beams) of the 
!  Bjorken-Mtingwa IBS formulation.
!
!  This formulation takes a very long time to evaluate.
!-

SUBROUTINE bjmt(ring, mode, rates)

IMPLICIT NONE

TYPE bjmt_struct
  REAL(rp) Lp(3,3)
  REAL(rp) Lh(3,3)
  REAL(rp) Lv(3,3)
  REAL(rp) L_sum(3,3)
END TYPE

! bjmt_struct is a bjmt common struct for passing matrices
! to the integrands this subroutine contains.

TYPE(bjmt_struct) :: bjmt_com

TYPE(modes_struct), INTENT(IN) :: mode
TYPE(ring_struct), INTENT(IN), target :: ring
TYPE(ibs_struct), INTENT(OUT) :: rates
TYPE(ele_struct), pointer :: ele

REAL(rp) sigma_p, emit_x, emit_y, sigma_z, E_tot
REAL(rp) gamma, KE, beta, beta_x, beta_y
REAL(rp) sigma_y
REAL(rp) Dx, Dy, Dxp, Dyp
REAL(rp) alpha_x, alpha_y, coulomb_log
REAL(rp) NB, big_A
REAL(rp) sigma_H, Hx, Hy
REAL(rp) sum_inv_Tz, sum_inv_Tx, sum_inv_Ty
REAL(rp) inv_Tz, inv_Tx, inv_Ty
REAL(rp) length_multiplier
INTEGER i,j

REAL(rp) phi_h, phi_v
REAL(rp) piece,last_piece,sum_x,sum_y,sum_z
REAL(rp) endpoint
LOGICAL far_enough

NB = ring%param%n_part
sigma_p = mode%sigE_E
emit_x = mode%a%emittance
emit_y = mode%b%emittance
sigma_z = mode%sig_z

E_tot = ring%ele_(0)%value(beam_energy$)
CALL convert_total_energy_to(E_tot, ring%param%particle, gamma, KE, beta)

sum_inv_Tz = 0.0
sum_inv_Tx = 0.0
sum_inv_Ty = 0.0
big_A = (r_e**2)*c_light*NB/64.0/(pi**2)/(beta**3)/(gamma**4)/ &
                                                emit_x/emit_y/sigma_z/sigma_p
DO i=1,ring%n_ele_use
  ele => ring%ele_(i)

  alpha_x = ele%x%alpha
  alpha_y = ele%y%alpha
  beta_x = ele%x%beta
  beta_y = ele%y%beta
  Dx = ele%x%eta
  Dy = ele%y%eta
  Dxp = ele%x%etap
  Dyp = ele%y%etap
  sigma_y = SQRT(beta_y*emit_y + (Dy*sigma_p)**2)

  coulomb_log = LOG( (gamma**2)*sigma_y*emit_x/r_e/beta_x )

  Hx = ( Dx**2 + (beta_x*Dxp + alpha_x*Dx)**2 ) / beta_x
  Hy = ( Dy**2 + (beta_y*Dyp + alpha_y*Dy)**2 ) / beta_y

  phi_h = Dxp + alpha_x*Dx/beta_x
  phi_v = Dyp + alpha_y*Dy/beta_y

  bjmt_com%Lp = 0.0
  bjmt_com%Lp(2,2) = 1.0
  bjmt_com%Lp = (gamma**2)/(sigma_p**2)*bjmt_com%Lp

  bjmt_com%Lh = 0.0
  bjmt_com%Lh(1,1) = 1.0
  bjmt_com%Lh(1,2) = -1.0*gamma*phi_h
  bjmt_com%Lh(2,1) = -1.0*gamma*phi_h
  bjmt_com%Lh(2,2) = (gamma**2)*Hx/beta_x
  bjmt_com%Lh = beta_x/emit_x*bjmt_com%Lh

  bjmt_com%Lv = 0.0
  bjmt_com%Lv(2,2) = (gamma**2)*Hy/beta_y
  bjmt_com%Lv(2,3) = -1.0*gamma*phi_v
  bjmt_com%Lv(3,2) = -1.0*gamma*phi_v
  bjmt_com%Lv(3,3) = 1.0
  bjmt_com%Lv = beta_y/emit_y*bjmt_com%Lv

  bjmt_com%L_sum = bjmt_com%Lp + bjmt_com%Lh + bjmt_com%Lv

  ! The bjmt integrals are slow to converge.  To help out qromo, the
  ! integrals are evaluated in segments, each segment 10 times longer
  ! than its predecessor.  When the integral of the segment contributes
  ! less than 1% to the total thus far, we call it converged.  This typically
  ! results in the integrals being evaluated out to 10^19 and higher.

  far_enough = .false.
  endpoint = 1.0E5
  piece = qromo(bjmt_int_p, 0.0_rp, endpoint, midpnt)
  sum_z = piece
  last_piece = piece
  DO WHILE(.not. far_enough) 
    piece = qromo(bjmt_int_p, endpoint, endpoint*1.0E1, midpnt)
    IF( (ABS(piece) .le. ABS(last_piece)) .and. (ABS(piece/sum_z) .le. .01) ) THEN
      far_enough = .true.
    ELSE
      sum_z = sum_z + piece
      last_piece = piece
      endpoint = endpoint*1.0E1
    ENDIF
  ENDDO
  inv_Tz=4.0*pi*big_A*coulomb_log*sum_z

  far_enough = .false.
  endpoint = 1.0E5
  piece = qromo(bjmt_int_h, 0.0_rp, endpoint, midpnt)
  sum_x = piece
  last_piece = piece
  DO WHILE(.not. far_enough) 
    piece = qromo(bjmt_int_h, endpoint, endpoint*1.0E1, midpnt)
    IF( (ABS(piece) .le. ABS(last_piece)) .and. (ABS(piece/sum_x) .le. .01) ) THEN
      far_enough = .true.
    ELSE
      sum_x = sum_x + piece
      last_piece = piece
      endpoint = endpoint*1.0E1
    ENDIF
  ENDDO
  inv_Tx=4.0*pi*big_A*coulomb_log*sum_x

  far_enough = .false.
  endpoint = 1.0E5
  piece = qromo(bjmt_int_v, 0.0_rp, endpoint, midpnt)
  sum_y = piece
  last_piece = piece
  DO WHILE(.not. far_enough) 
    piece = qromo(bjmt_int_v, endpoint, endpoint*1.0E1, midpnt)
    IF( (ABS(piece) .le. ABS(last_piece)) .and. (ABS(piece/sum_y) .le. .01) ) THEN
      far_enough = .true.
    ELSE
      sum_y = sum_y + piece
      last_piece = piece
      endpoint = endpoint*1.0E1
    ENDIF
  ENDDO
  inv_Ty=4.0*pi*big_A*coulomb_log*sum_y

  length_multiplier = ring%ele_(i)%value(l$)/2.0 + ring%ele_(i+1)%value(l$)/2.0
  sum_inv_Tz = sum_inv_Tz + inv_Tz * length_multiplier
  sum_inv_Tx = sum_inv_Tx + inv_Tx * length_multiplier
  sum_inv_Ty = sum_inv_Ty + inv_Ty * length_multiplier
ENDDO

rates%inv_Tz = sum_inv_Tz / ring%param%total_length
rates%inv_Tx = sum_inv_Tx / ring%param%total_length
rates%inv_Ty = sum_inv_Ty / ring%param%total_length

!------------------------------------------------------------------------
CONTAINS

FUNCTION bjmt_int_p(u)

real(rp), DIMENSION(:), INTENT(IN) :: u
real(rp), DIMENSION(size(u)) :: bjmt_int_p
    
real(rp) det, TrLi, TrInv, TrMult
real(rp) inv(3,3), ident(3,3), mult(3,3)
integer i

ident = 0.0
ident(1,1) = 1.0
ident(2,2) = 1.0
ident(3,3) = 1.0

TrLi = bjmt_com%Lp(1,1)+bjmt_com%Lp(2,2)+bjmt_com%Lp(3,3)

DO i=1, size(u)
  CALL mat_det(bjmt_com%L_sum + u(i)*ident, det)
  CALL mat_inverse(bjmt_com%L_sum + u(i)*ident, inv)
  TrInv = inv(1,1)+inv(2,2)+inv(3,3)
  mult = matmul(bjmt_com%Lp,inv)
  TrMult = mult(1,1)+mult(2,2)+mult(3,3)

  bjmt_int_p(i)=(u(i)**.5)/(det**.5)*(TrLi*TrInv-3.0*TrMult)
ENDDO

END FUNCTION bjmt_int_p

!------------------------------------------------------------------------
! CONTAINS

FUNCTION bjmt_int_h(u)

real(rp), DIMENSION(:), INTENT(IN) :: u
real(rp), DIMENSION(size(u)) :: bjmt_int_h

real(rp) det, TrLi, TrInv, TrMult
real(rp) inv(3,3), ident(3,3), mult(3,3)
integer i

ident = 0.0
ident(1,1) = 1.0
ident(2,2) = 1.0
ident(3,3) = 1.0

TrLi = bjmt_com%Lh(1,1)+bjmt_com%Lh(2,2)+bjmt_com%Lh(3,3)

DO i=1, size(u)
  CALL mat_det(bjmt_com%L_sum + u(i)*ident, det)
  CALL mat_inverse(bjmt_com%L_sum + u(i)*ident, inv)
  TrInv = inv(1,1)+inv(2,2)+inv(3,3)
  mult = matmul(bjmt_com%Lh,inv)
  TrMult = mult(1,1)+mult(2,2)+mult(3,3)

  bjmt_int_h(i)=(u(i)**.5)/(det**.5)*(TrLi*TrInv-3.0*TrMult)
ENDDO


END FUNCTION bjmt_int_h

!------------------------------------------------------------------------
! CONTAINS

FUNCTION bjmt_int_v(u)

real(rp), DIMENSION(:), INTENT(IN) :: u
real(rp), DIMENSION(size(u)) :: bjmt_int_v

real(rp) det, TrLi, TrInv, TrMult
real(rp) inv(3,3), ident(3,3), mult(3,3)
integer i

ident = 0.0
ident(1,1) = 1.0
ident(2,2) = 1.0
ident(3,3) = 1.0

TrLi = bjmt_com%Lv(1,1)+bjmt_com%Lv(2,2)+bjmt_com%Lv(3,3)

DO i=1, size(u)
  CALL mat_det(bjmt_com%L_sum + u(i)*ident, det)
  CALL mat_inverse(bjmt_com%L_sum + u(i)*ident, inv)
  TrInv = inv(1,1)+inv(2,2)+inv(3,3)
  mult = matmul(bjmt_com%Lv,inv)
  TrMult = mult(1,1)+mult(2,2)+mult(3,3)

  bjmt_int_v(i)=(u(i)**.5)/(det**.5)*(TrLi*TrInv-3.0*TrMult)
ENDDO

END FUNCTION bjmt_int_v
  
END SUBROUTINE bjmt

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  Subroutine bane(ring, mode, rates)
!
!  This is a private subroutine. To access this subroutine, call
!  ibs_rates.
!
!  This is an implementation of equations 10-15 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a high energy approximation of the Bjorken-Mtingwa IBS
!  formulation.
!-

SUBROUTINE bane(ring, mode, rates)

IMPLICIT NONE
 
! Elpha is a variable common between the bane subroutine and
! the integrand it contains.

TYPE(modes_struct), INTENT(IN) :: mode
TYPE(ring_struct), INTENT(IN), target :: ring
TYPE(ibs_struct), INTENT(OUT) :: rates
TYPE(ele_struct), pointer :: ele

REAL(rp) Elpha
REAL(rp) sigma_p, emit_x, emit_y, sigma_z, E_tot
REAL(rp) gamma, KE, beta, beta_x, beta_y
REAL(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
REAL(rp) Dx, Dy, Dxp, Dyp
REAL(rp) alpha_x, alpha_y, coulomb_log
REAL(rp) a, b, g_bane
REAL(rp) NB, big_A
REAL(rp) sigma_H, Hx, Hy
REAL(rp) sum_inv_Tz, sum_inv_Tx, sum_inv_Ty
REAL(rp) inv_Tz, inv_Tx, inv_Ty
REAL(rp) length_multiplier
INTEGER i

!

NB = ring%param%n_part
sigma_p = mode%sigE_E
emit_x = mode%a%emittance
emit_y = mode%b%emittance
sigma_z = mode%sig_z

E_tot = ring%ele_(0)%value(beam_energy$)
CALL convert_total_energy_to(E_tot, ring%param%particle, gamma, KE, beta)

sum_inv_Tz = 0.0
sum_inv_Tx = 0.0
sum_inv_Ty = 0.0
big_A=(r_e**2)*c_light*NB/16.0/(gamma**3)/(emit_x**(3./4.))/(emit_y**(3./4.))/sigma_z/(sigma_p**3)
DO i=1,ring%n_ele_use
  ele => ring%ele_(i)
                                                                                                     
  beta_x = ele%x%beta
  beta_y = ele%y%beta
  alpha_x = ele%x%alpha
  alpha_y = ele%y%alpha
  Dxp = ele%x%etap
  Dyp = ele%y%etap
  Dx = ele%x%eta
  Dy = ele%y%eta
  sigma_x_beta = SQRT(beta_x * emit_x)
  sigma_y_beta = SQRT(beta_y * emit_y)
  sigma_x = SQRT(sigma_x_beta**2 + (Dx**2)*(sigma_p**2))
  sigma_y = SQRT(sigma_y_beta**2 + (Dy**2)*(sigma_p**2))

  coulomb_log = LOG( (gamma**2)*sigma_y*emit_x/r_e/beta_x )

  Hx = ( (Dx**2) + (beta_x*Dxp + alpha_x*Dx)**2 ) / beta_x
  Hy = ( (Dy**2) + (beta_y*Dyp + alpha_y*Dy)**2 ) / beta_y
  sigma_H = 1.0/SQRT(1.0/(sigma_p**2) + Hx/emit_x + Hy/emit_y)

  a = sigma_H/gamma*SQRT(beta_x/emit_x)
  b = sigma_H/gamma*SQRT(beta_y/emit_y)

  Elpha = a/b
  g_bane = 2.*SQRT(Elpha)/pi*qromo(integrand, 0._rp, 9999._rp, midexp)

  inv_Tz = big_A*coulomb_log*sigma_H*g_bane*((beta_x*beta_y)**(-1./4.))
  inv_Tx = (sigma_p**2)*Hx/emit_x*inv_Tz
  inv_Ty = (sigma_p**2)*Hy/emit_y*inv_Tz

  length_multiplier = ring%ele_(i)%value(l$)/2.0 + ring%ele_(i+1)%value(l$)/2.0
  sum_inv_Tz = sum_inv_Tz + inv_Tz * length_multiplier
  sum_inv_Tx = sum_inv_Tx + inv_Tx * length_multiplier
  sum_inv_Ty = sum_inv_Ty + inv_Ty * length_multiplier
ENDDO
  
rates%inv_Tz = sum_inv_Tz / ring%param%total_length
rates%inv_Tx = sum_inv_Tx / ring%param%total_length
rates%inv_Ty = sum_inv_Ty / ring%param%total_length

!------------------------------------------------------------------------
CONTAINS

FUNCTION integrand(u)

real(rp), DIMENSION(:), INTENT(IN) :: u
real(rp), DIMENSION(size(u)) :: integrand

integer i

DO i=1, size(u)
  integrand(i) = 1 / (SQRT(1+u(i)**2)*SQRT(Elpha**2+u(i)**2))
ENDDO

END FUNCTION integrand

END SUBROUTINE bane

!------------------------------------------------------------------------
!------------------------------------------------------------------------
!------------------------------------------------------------------------
!+
!  subroutine cimp(ring, mode, rates)
!
!  This is a private subroutine. To access this subroutine, call
!  ibs_rates.
!
!  This is an implementation of equations 34,38-40 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a modified version of the Piwinski IBS formulation.
!  The integral (34) is handled with a piecewise interpolation generated
!  in mathematica.  The interpolation is accurate beyond 1% through it's
!  effective range (.0001 - 3000).
!
!  This is the quickest of the three IBS formuations in this module. 
!-

SUBROUTINE cimp(ring, mode, rates)

IMPLICIT NONE

TYPE(modes_struct), INTENT(IN) :: mode
TYPE(ring_struct), INTENT(IN), target :: ring
TYPE(ibs_struct), INTENT(OUT) :: rates
TYPE(ele_struct), pointer :: ele

REAL(rp) sigma_p, emit_x, emit_y, sigma_z, E_tot
REAL(rp) gamma, KE, beta, beta_x, beta_y
REAL(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
REAL(rp) Dx, Dy, Dxp, Dyp
REAL(rp) alpha_x, alpha_y, coulomb_log
REAL(rp) a, b
REAL(rp) NB, big_A
REAL(rp) sigma_H, Hx, Hy
REAL(rp) sum_inv_Tz, sum_inv_Tx, sum_inv_Ty
REAL(rp) inv_Tz, inv_Tx, inv_Ty
REAL(rp) length_multiplier
REAL(rp) g_ab,g_ba
REAL(rp) distance_s
INTEGER  i

NB = ring%param%n_part
sigma_p = mode%sigE_E
emit_x = mode%a%emittance
emit_y = mode%b%emittance
sigma_z = mode%sig_z

E_tot = ring%ele_(0)%value(beam_energy$)
CALL convert_total_energy_to(E_tot, ring%param%particle, gamma, KE, beta)

sum_inv_Tz = 0.0
sum_inv_Tx = 0.0
sum_inv_Ty = 0.0
distance_s = 0.0
  
big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(beta**3)/(gamma**4)/emit_x/emit_y/sigma_z/sigma_p
  
DO i=1,ring%n_ele_use   
  ele => ring%ele_(i)

  alpha_x = ele%x%alpha
  alpha_y = ele%y%alpha
  beta_x = ele%x%beta
  beta_y = ele%y%beta
  sigma_x_beta = SQRT(beta_x * emit_x)
  sigma_y_beta = SQRT(beta_y * emit_y)
  Dx = ele%x%eta
  Dy = ele%y%eta
  Dxp = ele%x%etap
  Dyp = ele%y%etap
  sigma_x = SQRT(sigma_x_beta**2 + (Dx**2)*(sigma_p**2))
  sigma_y = SQRT(sigma_y_beta**2 + (Dy**2)*(sigma_p**2))

  Hx = ( Dx**2 + (beta_x*Dxp + alpha_x*Dx)**2 ) / beta_x
  Hy = ( Dy**2 + (beta_y*Dyp + alpha_y*Dy)**2 ) / beta_y

  sigma_H = 1.0/SQRT( 1.0/(sigma_p**2)+ Hx/emit_x + Hy/emit_y )

  coulomb_log = LOG( (gamma**2)*sigma_y*emit_x/r_e/beta_x )

  a = sigma_H/gamma*SQRT(beta_x/emit_x)
  b = sigma_H/gamma*SQRT(beta_y/emit_y)

  g_ba = g(b/a)
  g_ab = g(a/b)

  inv_Tz = 2.*(pi**(3./2.))*big_A*(sigma_H**2)/(sigma_p**2) * &
      coulomb_log * ( g_ba/a + g_ab/b ) 
  inv_Tx = 2.*(pi**(3./2.))*big_A*coulomb_log*&
      (-a*g_ba + Hx*(sigma_H**2)/emit_x* &
      ( g_ba/a + g_ab/b ) ) 
  inv_Ty = 2.*(pi**(3./2.))*big_A*coulomb_log*&
      (-b*g_ab + Hy*(sigma_H**2)/emit_y* &
      ( g_ba/a + g_ab/b ) )

  length_multiplier = ring%ele_(i)%value(l$)/2.0 + ring%ele_(i+1)%value(l$)/2.0
  distance_s = distance_s + length_multiplier

  sum_inv_Tz = sum_inv_Tz + inv_Tz * length_multiplier
  sum_inv_Tx = sum_inv_Tx + inv_Tx * length_multiplier
  sum_inv_Ty = sum_inv_Ty + inv_Ty * length_multiplier
ENDDO
  
rates%inv_Tz = sum_inv_Tz / ring%param%total_length
rates%inv_Tx = sum_inv_Tx / ring%param%total_length
rates%inv_Ty = sum_inv_Ty / ring%param%total_length

!------------------------------------------------------------------------
CONTAINS

!+
!  function g(u)
!
!  This is an 13-degree piecewise polynomial interpolation of the
!  integral for the CIMP ibs formulation (equation 34 in "Intrabeam 
!  Scattering Formulas for High Energy Beams").
!  The segments were generated in mathematica and each is accurate
!  to beyond 1%.  The effective range of this interpolation is
!  .0001 to 3000, and it prints an error message when that range
!  is exceeded.
!-

FUNCTION g(u)
                                                                                       
REAL(rp), INTENT(IN) :: u
REAL(rp) :: g
REAL(rp), DIMENSION(0:10) ::  o,p,pa
REAL(rp), DIMENSION(0:13) ::  pb,qb,qa,q,ra,rb,r
 
o = (/ 0.233741, -0.00096173, 2.39076E-6, -3.82308E-9, 4.11144E-12, -3.039E-15, &
       1.5480865E-18, -5.3409137E-22, 1.1915518E-25, -1.550695E-29, 8.9365157E-34 /)
  
p = (/ 1.21744, -0.0453023, 0.00109043, -0.000017262, 1.85617E-7, -1.3792E-9, &
       7.08546E-12, -2.47044E-14, 5.57813E-17, -7.35484E-20, 4.2976E-23 /)
   
pa = (/ 2.4384485, -0.3342655, 0.0377896, -0.0032197625, .000203807, -9.5103569E-6, &
       3.22441912E-7, -7.7195818E-9, 1.2365680E-10, -1.1888946E-12, 5.1862247E-15 /)

pb = (/ 0.2509305915339517, 3.655027986910339, -3.538178517122364, 2.0378515154003622, &
       -0.8131215410970869, 0.23600479506508648, -0.0508941392205539, 0.008214142113823073, &
       -0.000988882504551898, 0.00008751565040446301, -5.526659141312207E-6, &
       2.3564088895171855E-7, -6.0770418816204626E-9, 7.15762060347684E-11 /)
 
qb = (/ -4.5128859976303835, 44.038623289844345, -201.50265131380843, 685.7422921643122, &
       -1733.9460775721304, 3272.7093350374225, -4637.724347309373, 4940.122242846527, &
       -3931.6350463124686, 2301.047404938234, -960.7211908043258, 270.6618330819612, &
       -46.088423620894815, 3.581345541476106 /)
 
qa= (/ -7.9542527228302164, 207.51519155414007, -4406.431743148087, 75592.0699189963, &
       -992747.0161009694, 9.983494244578896E6, -7.72264540055081E7, 4.593887542691004E8, &
       -2.0856924957059484E9, 7.103526742990105E9, -1.758262228582309E10, 2.9881056304916122E10, &
       -3.11970246734184E10, 1.5092353688380434E10 /)
  
q = (/ -10.427660034546474, 614.5741414073021, -37752.0301748355, 1.8467024380709291E6, &
       -6.787056202900247E7, 1.87800146779508E9, -3.9344782189063065E10, 6.245185690585099E11, &
       -7.460242970344108E12, 6.596837840619067E13, -4.186282233147678E14, 1.8023052362000775E15, &
       -4.713009373619726E15, 5.649407049035322E15 /)

ra = (/ -13.151132987122772, 2076.0767607226853, -439312.2888781302, 7.549171067948443E7, &
       -9.918579267112797E9, 9.97612075477821E11, -7.717556597617397E13, 4.591087739048907E15, &
       -2.0844941927298595E17, 7.099636713666537E18, -1.757338267407412E20, 2.986592440773308E21, &
       -3.1181759869417024E22, 1.5085206738847978E23 /)
  
rb = (/ -16.158250737444753, 7653.567961683397, -5.677248220833501E6, 3.2736621037820387E9, &
       -1.3885723919300298E12, 4.356123741560311E14, -1.018979723977809E17, 1.782083486703034E19, &
       -2.31837909173977E21, 2.209799968074478E23, -1.4978388984105023E25, 6.831861840241021E26, &
       -1.878873461985899E28, 2.3529610201488045E29 /)
 
r = (/ -21.460085854969584, 80412.06002488811, -6.294777676054858E8, 3.84371292991909E12, &
       -1.7311478313121348E16, 5.7791722579515556E19, -1.4411700564165495E23, +2.6910329100895636E26, &
       -3.742607580627404E29, 3.8178627111416235E32, -2.7722017774216816E35, 1.3556698189607729E38, &
       -4.000249119089512E30, 5.378513442816346E32 /) ! E10 added to these last two in formula

IF    (u .gt. 3000.0) THEN
  WRITE(*,*) "CRITICAL WARNING: interpolation range exceeded"
  g = 0.
ELSEIF(u .gt. 300.0) THEN
  g = o(0)+o(1)*u+o(2)*(u**2)+o(3)*(u**3)+o(4)*(u**4)+o(5)*(u**5)+ &
      o(6)*(u**6)+o(7)*(u**7)+o(8)*(u**8)+o(9)*(u**9)+o(10)*(u**10)
ELSEIF(u .gt. 30.0)  THEN
  g = p(0)+p(1)*u+p(2)*(u**2)+p(3)*(u**3)+p(4)*(u**4)+p(5)*(u**5)+ &
      p(6)*(u**6)+p(7)*(u**7)+p(8)*(u**8)+p(9)*(u**9)+p(10)*(u**10)
ELSEIF(u .gt. 10.0)  THEN
  g = pa(0)+pa(1)*u+pa(2)*(u**2)+pa(3)*(u**3)+pa(4)*(u**4)+pa(5)*(u**5)+ &
      pa(6)*(u**6)+pa(7)*(u**7)+pa(8)*(u**8)+pa(9)*(u**9)+pa(10)*(u**10)
ELSEIF(u .gt. 1.6)   THEN
  g = pb(0)+pb(1)*u+pb(2)*(u**2)+pb(3)*(u**3)+pb(4)*(u**4)+pb(5)*(u**5)+ &
      pb(6)*(u**6)+pb(7)*(u**7)+pb(8)*(u**8)+pb(9)*(u**9)+pb(10)*(u**10)+ &
      pb(11)*(u**11)+pb(12)*(u**12)+pb(13)*(u**13)
ELSEIF(u .gt. .20)   THEN
  g = qb(0)+qb(1)*u+qb(2)*(u**2)+qb(3)*(u**3)+qb(4)*(u**4)+qb(5)*(u**5)+ &
      qb(6)*(u**6)+qb(7)*(u**7)+qb(8)*(u**8)+qb(9)*(u**9)+qb(10)*(u**10)+ &
      qb(11)*(u**11)+qb(12)*(u**12)+qb(13)*(u**13)
ELSEIF(u .gt. .10)   THEN
  g = qa(0)+qa(1)*u+qa(2)*(u**2)+qa(3)*(u**3)+qa(4)*(u**4)+qa(5)*(u**5)+ &
      qa(6)*(u**6)+qa(7)*(u**7)+qa(8)*(u**8)+qa(9)*(u**9)+qa(10)*(u**10)+ &
      qa(11)*(u**11)+qa(12)*(u**12)+qa(13)*(u**13)
ELSEIF(u .gt. .02)  THEN
  g = q(0)+q(1)*u+q(2)*(u**2)+q(3)*(u**3)+q(4)*(u**4)+q(5)*(u**5)+ &
      q(6)*(u**6)+q(7)*(u**7)+q(8)*(u**8)+q(9)*(u**9)+q(10)*(u**10)+ &
      q(11)*(u**11)+q(12)*(u**12)+q(13)*(u**13)
ELSEIF(u .gt. .01)  THEN
  g = ra(0)+ra(1)*u+ra(2)*(u**2)+ra(3)*(u**3)+ra(4)*(u**4)+ra(5)*(u**5)+ &
      ra(6)*(u**6)+ra(7)*(u**7)+ra(8)*(u**8)+ra(9)*(u**9)+ra(10)*(u**10)+ &
      ra(11)*(u**11)+ra(12)*(u**12)+ra(13)*(u**13)
ELSEIF(u .gt. .001)  THEN
  g = rb(0)+rb(1)*u+rb(2)*(u**2)+rb(3)*(u**3)+rb(4)*(u**4)+rb(5)*(u**5)+ &
      rb(6)*(u**6)+rb(7)*(u**7)+rb(8)*(u**8)+rb(9)*(u**9)+rb(10)*(u**10)+ &
      rb(11)*(u**11)+rb(12)*(u**12)+rb(13)*(u**13)
ELSEIF(u .gt. .0001) THEN
  g = r(0)+r(1)*u+r(2)*(u**2)+r(3)*(u**3)+r(4)*(u**4)+r(5)*(u**5)+ &
      r(6)*(u**6)+r(7)*(u**7)+r(8)*(u**8)+r(9)*(u**9)+r(10)*(u**10)+ &
      r(11)*(u**11)+r(12)*(1.0E10)*(u**12)+r(13)*(1.0E10)*(u**13)
ELSE
  WRITE(*,*) "CRITICAL WARNING: interpolation range exceeded"
  g = 0.
ENDIF

END FUNCTION g

END SUBROUTINE cimp

END MODULE ibs_mod






