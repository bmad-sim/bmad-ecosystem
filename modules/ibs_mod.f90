#include "CESR_platform.inc"

MODULE ibs_mod

USE bmad_struct
USE bmad_interface
USE nr, only: qromo, midexp, midpnt


TYPE ibs_struct
  REAL(rp) inv_Tx
  REAL(rp) inv_Ty
  REAL(rp) inv_Tz
END TYPE

TYPE ibs_lifetime_struct
  REAL(rp) Tlx
  REAL(rp) Tly
  REAL(rp) Tlp
END TYPE

TYPE ibs_maxratio_struct
  REAL(rp) rx
  REAL(rp) ry
  REAL(rp) r_p
END TYPE

TYPE bjmt_struct
  REAL(rp) Lp(3,3)
  REAL(rp) Lh(3,3)
  REAL(rp) Lv(3,3)
  REAL(rp) L_sum(3,3)
END TYPE
!bjmt_struct is a bjmt common struct for passing matrices
!to the integrand functions of the bjmt subroutine
TYPE(bjmt_struct) :: bjmt_com

!Elpha is a variable common between the bane subroutine and
!the integrand function.
REAL(rp) Elpha

PUBLIC ibs_equilibrium
PUBLIC ibsequilibrium2
PUBLIC ibs_rates
PRIVATE cimp
PRIVATE bane
PRIVATE bjmt
PRIVATE mtto

CONTAINS
!+
!  Subroutine ibsequilibrium2(lat,inmode,ibsmode,formula,ratio,initial_blow_up)
!
!  Computes the equilibrium beam size using the equilibrium equations given in 
!  'Intrabeam scattering formulas for high energy beams' by Kubo et al.
!  In contrast to the ibs_equilibrium subroutine, ibsequilibrium2 takes into
!  account nominal vertical emittance due to dispersion.  The parameter
!  ratio is a statement about how much of the vertical beam size is expected
!  to be due to dispersion, and how much from coupling.  Ratios between
!  .50 and .85 are commonly used.
!  This method requires that the initial beam parameters be larger
!  than the expected equilibrium.  An initial_blow_up of 3 to 5 is 
!  a good place to start.
!
!  Modules needed:
!    use ibs_mod
!
!  Input:
!    lat             -- lat_struct: lattice for tracking
!      %param$n_part  -- Real: number of particles in bunch
!    inmode           -- normal_modes_struct: natural beam parameters 
!    formula          -- character*4: IBS formulation to use (see ibs_rates)
!    ratio            -- Real: Ratio of vert_emit_coupling / vert_emit_total
!    initial_blow_up  -- Real: Factor multiplied to all 3 bunch dimensions
!                              prior to starting iteration.
!
!  Output:
!    ibsmode          -- normal_modes_struct: beam parameters after IBS effects
!
!  See ibs_rates subroutine for available IBS rate formulas.
!-
SUBROUTINE ibsequilibrium2(lat,inmode,ibsmode,formula,ratio,initial_blow_up)

  IMPLICIT NONE
                                                                                                         
  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(normal_modes_struct), INTENT(IN) :: inmode
  TYPE(normal_modes_struct), INTENT(OUT) :: ibsmode
  CHARACTER*4, INTENT(IN) :: formula
  REAL(rp), INTENT(IN) :: ratio
  REAL(rp), INTENT(IN) :: initial_blow_up
  TYPE(ele_struct), pointer :: ele
  TYPE(ibs_struct) rates

  REAL(rp) time_for_one_turn
  REAL(rp) tau_x, tau_y, tau_z
  REAL(rp) Tx, Ty, Tz
  REAL(rp) xfactor,yfactor,zfactor
  REAL(rp) emit_x, emit_y, emit_z
  REAL(rp) emit_x0, emit_y0, emit_z0
  REAL(rp) advance, threshold
  REAL(rp) sigE_E, sigma_z, L_ratio
  REAL(rp) sigE_E0, sigma_z0
  REAL(rp) fpw
  REAL(rp) dex, d2ex, dey, d2ey, deE, d2eE
  REAL(rp) running_emit_x(1:10), running_emit_y(1:10), running_sigE_E(1:10)
  !REAL(rp) running_emit_x(1:3), running_emit_y(1:3), running_sigE_E(1:3)
  INTEGER half
  REAL(rp) runavg_emit_x_A, runavg_emit_y_A, runavg_sigE_E_A
  REAL(rp) runavg_emit_x_B, runavg_emit_y_B, runavg_sigE_E_B
  LOGICAL converged
  INTEGER counter, i

  CHARACTER(20) :: r_name = 'ibs_equilibrium2'

  !compute the SR betatron damping times
  time_for_one_turn = lat%param%total_length / c_light
  tau_x = time_for_one_turn / inmode%a%alpha_damp
  tau_y = time_for_one_turn / inmode%b%alpha_damp
  tau_z = time_for_one_turn / inmode%z%alpha_damp 

  !fpw is a simple way of modeling potential well bunch lengthing.
  !This factor is multiplied by the zero current L_ratio when 
  !determining bunch length from energy spread.
  fpw = 1.0 !+ .21/(12.0E9) * lat%param%n_part
  sigma_z0 = inmode%sig_z
  sigE_E0 = inmode%sigE_E 
  L_ratio = sigma_z0 / sigE_E0
  emit_x0 = inmode%a%emittance
  emit_y0 = inmode%b%emittance 

  ibsmode%a%emittance = emit_x0 * initial_blow_up
  ibsmode%b%emittance = emit_y0 * initial_blow_up
  ibsmode%sig_z = sigma_z0      * sqrt(initial_blow_up)
  ibsmode%sigE_E = sigE_E0      * sqrt(initial_blow_up)

  !compute equilibrium
  converged = .false.
  counter = 0
  !Advance is what percent of the way from the current emittance 
  !towards the equilibrium emittance the beam should be advanced on
  !each iteration.
  advance = .05
  !Changes in emittance between iterations less than threshold
  !indicate convergence.
  threshold = advance * .001

  DO i=1, SIZE(running_emit_x)
    running_emit_x(i) = 0.0
    running_emit_y(i) = 0.0
    running_sigE_E(i) = 0.0
  ENDDO
  half = SIZE(running_emit_x) / 2

  DO WHILE(.not.converged)
    CALL ibs_rates(lat,ibsmode,rates,formula)
    counter = counter + 1
    !It is possible that this method can give negative emittances
    !at some point in the iterative process, in which case the case
    !structure below will terminate the program.  If this happens, try
    !using different values for initial_blow_up.
    IF( rates%inv_Tx .gt. 0.0 ) THEN
      Tx = 1.0/rates%inv_Tx
      xfactor = 1.0/(1.0-(tau_x/Tx))
    ELSEIF( rates%inv_Tx .eq. 0.0 ) THEN
      xfactor = 1.0
    ELSE
      CALL out_io(s_abort$, r_name, &
         'FATAL ERROR: Negative emittance encountered: ', &
         'Try adjusting initial_blow_up.')
      STOP
    ENDIF
    IF( rates%inv_Ty .gt. 0.0 ) THEN
      Ty = 1.0/rates%inv_Ty
      yfactor = 1.0/(1.0-(tau_y/Ty))
    ELSEIF( rates%inv_Ty .eq. 0.0 ) THEN
      yfactor = 1.0
    ELSE
      CALL out_io(s_abort$, r_name, &
         'FATAL ERROR: Negative emittance encountered: ', &
         'Try adjusting initial_blow_up.')
      STOP
    ENDIF
    IF( rates%inv_Tz .gt. 0.0 ) THEN
      Tz = 1.0/rates%inv_Tz
      zfactor = 1/(1-(tau_z/Tz))
    ELSEIF( rates%inv_Tz .eq. 0.0 ) THEN
      zfactor = 1.0
    ELSE
      CALL out_io(s_abort$, r_name, &
         'FATAL ERROR: Negative emittance encountered: ', &
         'Try adjusting initial_blow_up.')
      STOP
    ENDIF

    emit_y = ibsmode%b%emittance + advance*( &
             ((1.0-ratio)*yfactor+ratio*xfactor)*emit_y0 - ibsmode%b%emittance )
    emit_x = ibsmode%a%emittance + advance*( xfactor*emit_x0 - ibsmode%a%emittance )
    sigE_E = ibsmode%sigE_E      + advance*( zfactor*sigE_E0 - ibsmode%sigE_E )

    running_emit_x = CSHIFT(running_emit_x, -1)
    running_emit_y = CSHIFT(running_emit_y, -1)
    running_sigE_E = CSHIFT(running_sigE_E, -1)
    running_emit_x(1) = emit_x
    running_emit_y(1) = emit_y
    running_sigE_E(1) = sigE_E
    runavg_emit_x_A = SUM(running_emit_x(1:half))/(half)
    runavg_emit_y_A = SUM(running_emit_y(1:half))/(half)
    runavg_sigE_E_A = SUM(running_sigE_E(1:half))/(half)
    runavg_emit_x_B = SUM(running_emit_x(half+1:))/(half)
    runavg_emit_y_B = SUM(running_emit_y(half+1:))/(half)
    runavg_sigE_E_B = SUM(running_sigE_E(half+1:))/(half)

    IF(counter .gt. SIZE(running_emit_x)) THEN
      IF( abs(runavg_emit_x_A/runavg_emit_x_B - 1.) .lt. threshold ) THEN
        IF( abs(runavg_emit_y_A/runavg_emit_y_B - 1.) .lt. threshold ) THEN
          IF( abs(runavg_sigE_E_A/runavg_sigE_E_B - 1.) .lt. threshold ) THEN
            converged = .true.
          ENDIF
        ENDIF
      ENDIF        
    ENDIF        

    !dex = running_emit_x(1) - running_emit_x(2)
    !dey = running_emit_y(1) - running_emit_y(2)
    !deE = running_sigE_E(1) - running_sigE_E(2)
    !d2ex = running_emit_x(1) - 2.0*running_emit_x(2) + running_emit_x(3)
    !d2ey = running_emit_y(1) - 2.0*running_emit_y(2) + running_emit_y(3)
    !d2eE = running_sigE_E(1) - 2.0*running_sigE_E(2) + running_sigE_E(3)
    !IF(counter .gt. SIZE(running_emit_x)) THEN
    !  IF( abs(dex) .lt. threshold ) THEN
    !  IF( abs(dey) .lt. threshold ) THEN
    !  IF( abs(deE) .lt. threshold ) THEN
    !  IF( abs(d2ex) .lt. threshold ) THEN
    !  IF( abs(d2ey) .lt. threshold ) THEN
    !  IF( abs(d2eE) .lt. threshold ) THEN
    !    converged = .true.
    !  ENDIF
    !  ENDIF
    !  ENDIF        
    !  ENDIF
    !  ENDIF
    !  ENDIF        
    !ENDIF        

    IF(.not.converged) THEN
      ibsmode%a%emittance = emit_x
      ibsmode%b%emittance = emit_y
      ibsmode%sigE_E = sigE_E 
      ibsmode%sig_z = L_ratio * sigE_E * fpw
    ENDIF

!    WRITE(55,'(4e20.11)') emit_x, emit_y, &
!               sigE_E, L_ratio * sigE_E * fpw

  ENDDO

END SUBROUTINE ibsequilibrium2

!+
!  Subroutine ibs_equilibrium(lat,inmode,ibsmode,formula,coupling)
!
!  Computes equilibrium beam sizes taking into account
!  radiation damping, IBS growth rates, and coupling.
!
!  Modules needed:
!    use ibs_mod
!
!  Input:
!    lat             -- lat_struct: lattice for tracking
!      %param$n_part  -- Real: number of particles in bunch
!    inmode           -- normal_modes_struct: natural beam parameters 
!    formula          -- character*4: IBS formulation to use (see ibs_rates)
!    coupling         -- real: horizontal to vertical emittanc coupling
!
!  Output:
!    ibsmode          -- normal_modes_struct: beam parameters after IBS effects
!
!  See ibs_rates subroutine for available IBS rate formulas.
!-
SUBROUTINE ibs_equilibrium(lat,inmode,ibsmode,formula,coupling)

  IMPLICIT NONE
                                                                                                         
  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(normal_modes_struct), INTENT(IN) :: inmode
  TYPE(normal_modes_struct), INTENT(OUT) :: ibsmode
  CHARACTER*4, INTENT(IN) :: formula
  REAL(rp), INTENT(IN), OPTIONAL :: coupling
  TYPE(ele_struct), pointer :: ele
  TYPE(ibs_struct) rates
  TYPE(normal_modes_struct) :: naturalmode

  REAL(rp) time_for_one_turn
  REAL(rp) tau_x, tau_y, tau_z
  REAL(rp) Tx, Ty, Tz
  REAL(rp) emit_x, emit_y, emit_z
  REAL(rp) emit_x0, emit_y0, emit_z0
  REAL(rp) dxdt, dydt, dzdt, dT, reps
  REAL(rp) threshold
  REAL(rp) sigE_E0, sigma_z0, L_ratio
  REAL(rp) ka, ka_one, ka_small
  LOGICAL converged
  INTEGER counter

  !natural mode here means emittances before IBS effects
  naturalmode = inmode
  IF(present(coupling)) THEN
    ka = coupling
    naturalmode%b%emittance = coupling * inmode%a%emittance
  ELSE
    ka = inmode%b%emittance / inmode%a%emittance
  ENDIF

  !compute the SR betatron damping times
  time_for_one_turn = lat%param%total_length / c_light
  tau_x = time_for_one_turn / naturalmode%a%alpha_damp
  tau_y = time_for_one_turn / naturalmode%b%alpha_damp
  tau_z = time_for_one_turn / naturalmode%z%alpha_damp 

  sigma_z0 = naturalmode%sig_z
  sigE_E0 = naturalmode%sigE_E
  emit_x0 = naturalmode%a%emittance
  emit_y0 = naturalmode%b%emittance
  emit_z0 = sigma_z0 * sigE_E0
  L_ratio = sigma_z0 / sigE_E0
  threshold = .00001 !fractional changes in emittance smaller than this
                     !indicate convergence
  converged = .false.
  dT = tau_x / 40.0 !Time to advance per iteration
  ka_one = 1. / (1.+ka) !Used in determining de_dt
  ka_small = ka / (1.+ka) !Used in determining de_dt

  counter = 0
  ibsmode = naturalmode
  DO WHILE(.not.converged)
    CALL ibs_rates(lat,ibsmode,rates,formula)
    counter = counter + 1
    Tx = 1.0/rates%inv_Tx
    Ty = 1.0/rates%inv_Ty
    Tz = 1.0/rates%inv_Tz
    emit_x = ibsmode%a%emittance
    emit_y = ibsmode%b%emittance
    emit_z = ibsmode%sigE_E*ibsmode%sig_z

    !Compute change in emittance per time for x,y,z dimensions, taking
    !into account radiation damping, IBS blow-up, and x-y coupling
    dxdt = -(emit_x-ka_one*emit_x0)*2./tau_x + emit_x*2./Tx
    dydt = -(emit_y-ka_small*emit_x)*2./tau_y + emit_y*2./Ty
    dzdt = -(emit_z-emit_z0)*2./tau_z + 2./Tz*emit_z

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

!+
!  Subroutine ibs_lifetime(lat, mode, maxratio, lifetime, formula)
!
!  This module computes the beam lifetime due to
!  the diffusion process according to equation 12
!  from page 129 of The Handbook for Accelerator
!  Physics and Engineering 2nd edition.
!
!  Input:
!    lat             -- lat_struct: lattice for tracking
!    mode             -- normal_modes_struct: beam parameters 
!    maxratio(ibs_maxratio_struct)  Ax,y,p/sigma_x,y,p where Ax,y,p
!                 is the maximum sigma.  Note that this quantity
!                 is just the ratio, not the ratio squared.  For
!                 example, maxratio%Rx = 1.1 says that the maximum
!                 acceptable beamsize is 10% larger than the beamsize
!                 before IBS effects.
!    formula(CHAR*4) See the ibs_rates subroutine for available formulas.
!                 formula = 'cimp' is a good choice.
!  Output:
!    lifetime(ibs_lifetime_struct) --structure returning IBS lifetimes
!-
SUBROUTINE ibs_lifetime(lat,mode,maxratio,lifetime,formula)
  IMPLICIT NONE

  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(normal_modes_struct), INTENT(IN) :: mode
  TYPE(ibs_maxratio_struct), INTENT(IN) :: maxratio
  TYPE(ibs_lifetime_struct), INTENT(OUT) :: lifetime
  CHARACTER*4, INTENT(IN) :: formula

  TYPE(ibs_struct) rates

  REAL(rp) Rx, Ry, R_p

  Rx = maxratio%rx**2
  Ry = maxratio%ry**2
  R_p = maxratio%r_p**2

  CALL ibs_rates(lat, mode, rates, formula)

  lifetime%Tlx = exp(Rx)/2/Rx/rates%inv_Tx
  lifetime%Tly = exp(Ry)/2/Ry/rates%inv_Ty
  lifetime%Tlp = exp(R_p)/2/R_p/rates%inv_Tz
END SUBROUTINE ibs_lifetime

!+
!  Subroutine ibs_rates(lat, mode, rates, formula)
!
!  Calculates IBS risetimes for given lat and mode.
!  This is basically a front-end for the various formulas 
!  available in this module of calculating IBS rates.
!
!  Available IBS formulas:
!    cimp - Modified Piwinski
!    bane - Bane approximation of Bjorken-Mtingwa formulation
!    bjmt - Bjorken-Mtingwa formulation general to bunched beams (time consuming)
!    mtto - Mtingwa-Tollerstrup formulation
!
!  Input:
!    lat             -- lat_struct: lattice for tracking
!      %param$n_part  -- Real: number of particles in bunch
!    mode             -- normal_modes_struct: beam parameters 
!    formula          -- character*4: IBS formulation to use
!
!  Output:
!    rates          -- ibs_struct: ibs rates in x,y,z 
!-
SUBROUTINE ibs_rates(lat, mode, rates, formula)
  IMPLICIT NONE

  TYPE(normal_modes_struct), INTENT(IN) :: mode
  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(ibs_struct), INTENT(OUT) :: rates
  CHARACTER*4, INTENT(IN) ::  formula

  IF(mode%a%emittance .le. 0.0) THEN
    WRITE(*,*) "Negative or zero emittance detected: stopping"
    STOP
  ENDIF
  IF(mode%b%emittance .le. 0.0) THEN
    WRITE(*,*) "Negative or zero emittance detected: stopping"
    STOP
  ENDIF
  IF(mode%sig_z*mode%sigE_E .le. 0.0) THEN
    WRITE(*,*) "Negative or zero emittance detected: stopping"
    STOP
  ENDIF

  IF(formula == 'cimp') THEN
    CALL cimp(lat, mode, rates)
  ELSEIF(formula == 'bjmt') THEN
    CALL bjmt(lat, mode, rates)
  ELSEIF(formula == 'bane') THEN
    CALL bane(lat, mode, rates)
  ELSEIF(formula == 'mtto') THEN
    CALL mtto(lat, mode, rates)
  ELSE
    WRITE(*,*) "Invalid IBS formula selected ... returning zero"
    rates%inv_Tx = 0.0
    rates%inv_Ty = 0.0
    rates%inv_Tz = 0.0
  ENDIF
END SUBROUTINE ibs_rates

!+
!  subroutine bjmt(lat, mode, rates)
!
!  This is a private subroutine.  To access this subroutine, call
!  ibs_rates.
!
!  This is an implementation of equations 1-9 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is the most general form (for bunched beams) of the 
!  Bjorken-Mtingwa IBS formulation.
!
!  This formulation takes a very long time to evaluate.
!-
SUBROUTINE bjmt(lat, mode, rates)
  IMPLICIT NONE

  TYPE(normal_modes_struct), INTENT(IN) :: mode
  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(ibs_struct), INTENT(OUT) :: rates
  TYPE(ele_struct), pointer :: ele

  REAL(rp) sigma_p, emit_x, emit_y, sigma_z, E_tot
  REAL(rp) gamma, KE, beta, beta_a, beta_b
  REAL(rp) sigma_y
  REAL(rp) Dx, Dy, Dxp, Dyp
  REAL(rp) alpha_a, alpha_b, coulomb_log
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

  NB = lat%param%n_part
  sigma_p = mode%sigE_E
  emit_x = mode%a%emittance
  emit_y = mode%b%emittance
  sigma_z = mode%sig_z

  E_tot = lat%ele(0)%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)

  sum_inv_Tz = 0.0
  sum_inv_Tx = 0.0
  sum_inv_Ty = 0.0
  big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(beta**3)/(gamma**4)/emit_x/emit_y/sigma_z/sigma_p
  DO i=1,lat%n_ele_track
    ele => lat%ele(i)

    alpha_a = ele%a%alpha
    alpha_b = ele%b%alpha
    beta_a = ele%a%beta
    beta_b = ele%b%beta
    Dx = ele%a%eta
    Dy = ele%b%eta
    Dxp = ele%a%etap
    Dyp = ele%b%etap
    sigma_y = SQRT(beta_b*emit_y + (Dy*sigma_p)**2)

    coulomb_log = LOG( (gamma**2)*sigma_y*emit_x/r_e/beta_a )

    Hx = ( Dx**2 + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
    Hy = ( Dy**2 + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b

    phi_h = Dxp + alpha_a*Dx/beta_a
    phi_v = Dyp + alpha_b*Dy/beta_b

    bjmt_com%Lp = 0.0
    bjmt_com%Lp(2,2) = 1.0
    bjmt_com%Lp = (gamma**2)/(sigma_p**2)*bjmt_com%Lp

    bjmt_com%Lh = 0.0
    bjmt_com%Lh(1,1) = 1.0
    bjmt_com%Lh(1,2) = -1.0*gamma*phi_h
    bjmt_com%Lh(2,1) = -1.0*gamma*phi_h
    bjmt_com%Lh(2,2) = (gamma**2)*Hx/beta_a
    bjmt_com%Lh = beta_a/emit_x*bjmt_com%Lh

    bjmt_com%Lv = 0.0
    bjmt_com%Lv(2,2) = (gamma**2)*Hy/beta_b
    bjmt_com%Lv(2,3) = -1.0*gamma*phi_v
    bjmt_com%Lv(3,2) = -1.0*gamma*phi_v
    bjmt_com%Lv(3,3) = 1.0
    bjmt_com%Lv = beta_b/emit_y*bjmt_com%Lv

    bjmt_com%L_sum = bjmt_com%Lp + bjmt_com%Lh + bjmt_com%Lv

    !The bjmt integrals are slow to converge.  To help out qromo, the
    !integrals are evaluated in segments, each segment 10 times longer
    !than its predecessor.  When the integral of the segment contributes
    !less than 1% to the total thus far, we call it converged.  This typically
    !results in the integrals being evaluated out to 10^19 and higher.

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

    length_multiplier = lat%ele(i)%value(l$)/2.0 + lat%ele(i+1)%value(l$)/2.0
    sum_inv_Tz = sum_inv_Tz + inv_Tz * length_multiplier
    sum_inv_Tx = sum_inv_Tx + inv_Tx * length_multiplier
    sum_inv_Ty = sum_inv_Ty + inv_Ty * length_multiplier
  ENDDO

  rates%inv_Tz = sum_inv_Tz / lat%param%total_length
  rates%inv_Tx = sum_inv_Tx / lat%param%total_length
  rates%inv_Ty = sum_inv_Ty / lat%param%total_length

END SUBROUTINE bjmt

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
!-
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
!-
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
  

!+
!  subroutine bane(lat, mode, rates)
!
!  This is a private subroutine. To access this subroutine, call
!  ibs_rates.
!
!  This is an implementation of equations 10-15 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a high energy approximation of the Bjorken-Mtingwa IBS
!  formulation.
!-
SUBROUTINE bane(lat, mode, rates)

  IMPLICIT NONE
 

  TYPE(normal_modes_struct), INTENT(IN) :: mode
  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(ibs_struct), INTENT(OUT) :: rates
  TYPE(ele_struct), pointer :: ele

  REAL(rp) sigma_p, emit_x, emit_y, sigma_z, E_tot
  REAL(rp) gamma, KE, beta, beta_a, beta_b
  REAL(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
  REAL(rp) Dx, Dy, Dxp, Dyp
  REAL(rp) alpha_a, alpha_b, coulomb_log
  REAL(rp) a, b, g_bane
  REAL(rp) NB, big_A
  REAL(rp) sigma_H, Hx, Hy
  REAL(rp) sum_inv_Tz, sum_inv_Tx, sum_inv_Ty
  REAL(rp) inv_Tz, inv_Tx, inv_Ty
  REAL(rp) length_multiplier
  INTEGER i

  NB = lat%param%n_part
  sigma_p = mode%sigE_E
  emit_x = mode%a%emittance
  emit_y = mode%b%emittance
  sigma_z = mode%sig_z

  E_tot = lat%ele(0)%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)

  sum_inv_Tz = 0.0
  sum_inv_Tx = 0.0
  sum_inv_Ty = 0.0
  big_A=(r_e**2)*c_light*NB/16.0/(gamma**3)/(emit_x**(3./4.))/(emit_y**(3./4.))/sigma_z/(sigma_p**3)
  DO i=1,lat%n_ele_track
    ele => lat%ele(i)
                                                                                                     
    beta_a = ele%a%beta
    beta_b = ele%b%beta
    alpha_a = ele%a%alpha
    alpha_b = ele%b%alpha
    Dxp = ele%a%etap
    Dyp = ele%b%etap
    Dx = ele%a%eta
    Dy = ele%b%eta
    sigma_x_beta = SQRT(beta_a * emit_x)
    sigma_y_beta = SQRT(beta_b * emit_y)
    sigma_x = SQRT(sigma_x_beta**2 + (Dx**2)*(sigma_p**2))
    sigma_y = SQRT(sigma_y_beta**2 + (Dy**2)*(sigma_p**2))
                                                                                                       
    coulomb_log = LOG( (gamma**2)*sigma_y*emit_x/r_e/beta_a )
                                                                                                     
    Hx = ( (Dx**2) + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
    Hy = ( (Dy**2) + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b
    !-for atf-! Hy = 5.13E-7
    sigma_H = 1.0/SQRT(1.0/(sigma_p**2) + Hx/emit_x + Hy/emit_y)

    a = sigma_H/gamma*SQRT(beta_a/emit_x)
    b = sigma_H/gamma*SQRT(beta_b/emit_y)
                                                                                                     
    Elpha = a/b
    g_bane = 2.*SQRT(Elpha)/pi*qromo(integrand, 0._rp, 9999._rp, midexp)
                                                                                                 
    inv_Tz = big_A*coulomb_log*sigma_H*g_bane*((beta_a*beta_b)**(-1./4.))
    inv_Tx = (sigma_p**2)*Hx/emit_x*inv_Tz
    inv_Ty = (sigma_p**2)*Hy/emit_y*inv_Tz
                                                                                                 
    length_multiplier = lat%ele(i)%value(l$)/2.0 + lat%ele(i+1)%value(l$)/2.0
    sum_inv_Tz = sum_inv_Tz + inv_Tz * length_multiplier
    sum_inv_Tx = sum_inv_Tx + inv_Tx * length_multiplier
    sum_inv_Ty = sum_inv_Ty + inv_Ty * length_multiplier
  ENDDO
  
  rates%inv_Tz = sum_inv_Tz / lat%param%total_length
  rates%inv_Tx = sum_inv_Tx / lat%param%total_length
  rates%inv_Ty = sum_inv_Ty / lat%param%total_length

END SUBROUTINE bane

FUNCTION integrand(u)
real(rp), DIMENSION(:), INTENT(IN) :: u
real(rp), DIMENSION(size(u)) :: integrand

integer i

DO i=1, size(u)
  integrand(i) = 1 / (SQRT(1+u(i)**2)*SQRT(Elpha**2+u(i)**2))
ENDDO
END FUNCTION integrand


!+
!  subroutine cimp(lat, mode, rates)
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
SUBROUTINE cimp(lat, mode, rates)

  IMPLICIT NONE

  TYPE(normal_modes_struct), INTENT(IN) :: mode
  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(ibs_struct), INTENT(OUT) :: rates
  TYPE(ele_struct), pointer :: ele

  REAL(rp) sigma_p, emit_x, emit_y, sigma_z, E_tot
  REAL(rp) gamma, KE, beta, beta_a, beta_b
  REAL(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
  REAL(rp) Dx, Dy, Dxp, Dyp
  REAL(rp) alpha_a, alpha_b, coulomb_log
  REAL(rp) a, b
  REAL(rp) NB, big_A
  REAL(rp) sigma_H, Hx, Hy
  REAL(rp) sum_inv_Tz, sum_inv_Tx, sum_inv_Ty
  REAL(rp) inv_Tz, inv_Tx, inv_Ty
  REAL(rp) length_multiplier
  REAL(rp) g_ab,g_ba
  REAL(rp) distance_s
  INTEGER  i

  NB = lat%param%n_part
  sigma_p = mode%sigE_E
  emit_x = mode%a%emittance
  emit_y = mode%b%emittance
  sigma_z = mode%sig_z

  E_tot = lat%ele(0)%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)

  sum_inv_Tz = 0.0
  sum_inv_Tx = 0.0
  sum_inv_Ty = 0.0
  distance_s = 0.0

  big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(beta**3)/(gamma**4)/emit_x/emit_y/sigma_z/sigma_p

  !OPEN(88, FILE='contribs.dat',STATUS='REPLACE')
  DO i=1,lat%n_ele_track   
    ele => lat%ele(i)

    alpha_a = ele%a%alpha
    alpha_b = ele%b%alpha
    beta_a = ele%a%beta
    beta_b = ele%b%beta
    sigma_x_beta = SQRT(beta_a * emit_x)
    sigma_y_beta = SQRT(beta_b * emit_y)
    Dx = ele%a%eta
    Dy = ele%b%eta
    Dxp = ele%a%etap
    Dyp = ele%b%etap
    sigma_x = SQRT(sigma_x_beta**2 + (Dx**2)*(sigma_p**2))
    sigma_y = SQRT(sigma_y_beta**2 + (Dy**2)*(sigma_p**2))

    Hx = ( Dx**2 + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
    Hy = ( Dy**2 + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b

    sigma_H = 1.0/SQRT( 1.0/(sigma_p**2)+ Hx/emit_x + Hy/emit_y )

    coulomb_log = LOG( (gamma**2)*sigma_y*emit_x/r_e/beta_a )

    a = sigma_H/gamma*SQRT(beta_a/emit_x)
    b = sigma_H/gamma*SQRT(beta_b/emit_y)

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

    length_multiplier = lat%ele(i)%value(l$)/2.0 + lat%ele(i+1)%value(l$)/2.0
    distance_s = distance_s + length_multiplier

!    WRITE(88,FMT="(ES8.2,' ',ES11.4,' ',ES11.4,' ',ES11.4)") &
!        distance_s,inv_Tx,beta_b
    sum_inv_Tz = sum_inv_Tz + inv_Tz * length_multiplier
    sum_inv_Tx = sum_inv_Tx + inv_Tx * length_multiplier
    sum_inv_Ty = sum_inv_Ty + inv_Ty * length_multiplier
  ENDDO
!  CLOSE(88)

  rates%inv_Tz = sum_inv_Tz / lat%param%total_length
  rates%inv_Tx = sum_inv_Tx / lat%param%total_length
  rates%inv_Ty = sum_inv_Ty / lat%param%total_length
END SUBROUTINE cimp

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
REAL(rp), DIMENSION(0:13) ::  o,p,pa,pb,qb,qa,q,ra,rb,r

o = (/ 0.27537408880308967,-0.0014411559647655005,4.760237316963749E-6,-1.056200254057877E-8, &
       1.652957676252096E-11,-1.8807607644579788E-14,1.582497859744359E-17,-9.911417850544026E-21, &
       4.606212267581107E-24,-1.5661683056801226E-27,3.7836532348260456E-31,-6.148053455724941E-35, &
       6.021954146486214E-29,-2.6856345839264216E-33 /) !last two terms have E-10 multiplied in function

p = (/ 1.3728514796955633,-0.06329159139596323,0.0019843996282577487,-0.00004282113130209425, &
       6.588299268041797E-7,-7.411644739022654E-9,6.186473337975396E-11,-3.851874112809632E-13, &
       1.7820926145729545E-15,-6.038143148108221E-18,1.454670683108795E-20,-2.358356775328005E-23, &
       2.305700018891372E-26,-1.026694954987225E-29 /)

pa = (/ 2.4560269972177395,-0.34483385959575796,0.04060452194930618,-0.0036548931452558006, &
        0.00024651396218017436,-0.000012205596049761074,4.224156486889531E-7,-8.746709793365587E-9, &
	2.0010958003696878E-11,5.293109492480939E-12,-1.8692549451553078E-13,3.313011636321362E-15, &
        -3.191257818162888E-17,1.333629848206322E-19 /)

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
      o(6)*(u**6)+o(7)*(u**7)+o(8)*(u**8)+o(9)*(u**9)+o(10)*(u**10)+ &
      o(11)*(u**11)+o(12)*(1.0E-10)*(u**12)+o(13)*(1.0E-10)*(u**13)
ELSEIF(u .gt. 30.0)  THEN
  g = p(0)+p(1)*u+p(2)*(u**2)+p(3)*(u**3)+p(4)*(u**4)+p(5)*(u**5)+ &
      p(6)*(u**6)+p(7)*(u**7)+p(8)*(u**8)+p(9)*(u**9)+p(10)*(u**10)+ &
      p(11)*(u**11)+p(12)*(u**12)+p(13)*(u**13)
ELSEIF(u .gt. 10.0)  THEN
  g = pa(0)+pa(1)*u+pa(2)*(u**2)+pa(3)*(u**3)+pa(4)*(u**4)+pa(5)*(u**5)+ &
      pa(6)*(u**6)+pa(7)*(u**7)+pa(8)*(u**8)+pa(9)*(u**9)+pa(10)*(u**10)+ &
      pa(11)*(u**11)+pa(12)*(u**12)+pa(13)*(u**13)
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


!+
!  subroutine mtto(lat, mode, rates)
!
!  This is a private subroutine. To access this subroutine, call
!  ibs_rates.
!
!  This subroutine is an implementation of equations 54-56 from
!  "Intrabeam Scattering Formulae for Asymptotic Beams with
!  Unequal Horizontal and Vertical Emittances" by Mtingwa and
!  Tollerstrup.  This formula is useful in that it can be
!  nicely transparent when you're figuring out how the lattice
!  properties contribute to the IBS rates.
!
!-
SUBROUTINE mtto(lat, mode, rates)

  IMPLICIT NONE

  TYPE(normal_modes_struct), INTENT(IN) :: mode
  TYPE(lat_struct), INTENT(IN), target :: lat
  TYPE(ibs_struct), INTENT(OUT) :: rates
  TYPE(ele_struct), pointer :: ele

  REAL(rp) sigma_p, emit_x, emit_y, sigma_z, E_tot
  REAL(rp) gamma, KE, beta, beta_a, beta_b
  REAL(rp) sigma_x_beta, sigma_y_beta
  REAL(rp) sigma_x_eta, sigma_y_eta
  REAL(rp) Dx
  REAL(rp) a, b, c, T, big_A
  REAL(rp) NB
  REAL(rp) inv_Tz, inv_Tx, inv_Ty
  REAL(rp) length_multiplier
  REAL(rp) g_b,g_1_b
  INTEGER  i

  NB = lat%param%n_part
  sigma_p = mode%sigE_E
  emit_x = mode%a%emittance
  emit_y = mode%b%emittance
  sigma_z = mode%sig_z

  E_tot = lat%ele(0)%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, beta)

  big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(beta**3)/(gamma**4)/emit_x/emit_y/sigma_z/sigma_p

  beta_a = 0.0
  beta_b = 0.0
  Dx = 0.0
 
  DO i=1,lat%n_ele_track   
    ele => lat%ele(i)
    length_multiplier = lat%ele(i)%value(l$)/2.0 + lat%ele(i+1)%value(l$)/2.0
    !Compute average lattice beta_a, beta_b, and Dx
    beta_a = beta_a + (ele%a%beta*length_multiplier)
    beta_b = beta_b + (ele%b%beta*length_multiplier)
    Dx = Dx + (ele%a%eta*length_multiplier)
  ENDDO

  beta_a = beta_a / lat%param%total_length
  beta_b = beta_b / lat%param%total_length
  Dx = Dx / lat%param%total_length

  T = beta_a*emit_x / ( beta_a*emit_x + (Dx**2)*(sigma_p**2) )
  a = 1./Sqrt(T)/sigma_p*Sqrt((gamma**2)*emit_y/beta_b)
  b = Sqrt(emit_y*beta_b/emit_x/beta_a)
  c = ((gamma**2)*emit_y/beta_b/r_e)**(1./2.) * &
      (Sqrt(6*pi)*sigma_z/NB)**(1./6.) * &
      (36.*(pi*pi)*(gamma**2)*emit_x*emit_y*beta_a*beta_b/T)**(1./12.)

  !This is the same g used by the cimp formula.
  g_b = g(b)
  g_1_b = g(1./b)

  rates%inv_Ty = -8.*(pi**(3./2.))*big_A/a*g_b*Log(c)
  rates%inv_Tx =  8.*(pi**(3./2.))*big_A*a*(1.-T)*Log(c) * &
                  ( g_b + (1./b)*g_1_b)
  rates%inv_Tz =  8.*(pi**(3./2.))*big_A*a*T*Log(c) * &
                  ( g_b + (1./b)*g_1_b)

END SUBROUTINE mtto

END MODULE ibs_mod






