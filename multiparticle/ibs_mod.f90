MODULE ibs_mod

USE bmad
USE nr
USE fgsl
USE, INTRINSIC :: iso_c_binding

TYPE ibs_sim_param_struct
  REAL(rp) tau_a        ! horizontal damping rate (needed for coulomb log tail cut)
  INTEGER clog_to_use   ! 1=Wolski, 2=Raubenheimer Tail Cut, 3=Bane Tail Cut
  LOGICAL set_dispersion
  REAL(rp) eta_set
  REAL(rp) etap_set
  REAL(rp) inductance
  REAL(rp) resistance
!  REAL(rp) clog         ! Smuggle Coulomb Log out of this depths of the calculation.
END TYPE

TYPE ibs_struct
  REAL(rp) inv_Ta
  REAL(rp) inv_Tb
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

REAL(fgsl_double), PARAMETER :: eps7 = 1.0d-7
INTEGER(fgsl_size_t), PARAMETER :: limit = 1000_fgsl_size_t

PUBLIC ibs_equib_rlx
PUBLIC ibs_equib_der
PUBLIC ibs_rates1turn
PUBLIC ibs_blowup1turn
PUBLIC cimp1
PUBLIC bane1
PUBLIC bjmt1

CONTAINS

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
!  Subroutine ibsequilibrium2(lat,inmode,ibsmode,formula,ratio,initial_blow_up,granularity)
!  Iterates to equilibrium beam conditions using relaxation method
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
SUBROUTINE ibs_equib_rlx(lat,ibs_sim_params,inmode,ibsmode,formula,ratio,initial_blow_up,granularity)

  IMPLICIT NONE
                                                                                                         
  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) ibs_sim_params
  TYPE(normal_modes_struct) :: inmode
  TYPE(normal_modes_struct) :: ibsmode
  REAL(rp) :: granularity
  CHARACTER*4 :: formula
  REAL(rp) :: ratio
  REAL(rp) :: initial_blow_up
  TYPE(ibs_struct) rates

  REAL(rp) time_for_one_turn
  REAL(rp) cur
  REAL(rp) alpha_c
  REAL(rp) pwd_ratio
  REAL(rp) sigma_z_nat
  REAL(rp) sigma_z_vlassov
  REAL(rp) tau_a, tau_b, tau_z
  REAL(rp) Ta, Tb, Tz
  REAL(rp) afactor,bfactor,zfactor
  REAL(rp) emit_a, emit_b
  REAL(rp) emit_a0, emit_b0
  REAL(rp) advance, threshold
  REAL(rp) sigE_E, sigE_E0, sigma_z0, L_ratio
  REAL(rp) running_emit_x(1:10), running_emit_y(1:10), running_sigE_E(1:10)
  INTEGER half
  REAL(rp) runavg_emit_x_A, runavg_emit_y_A, runavg_sigE_E_A
  REAL(rp) runavg_emit_x_B, runavg_emit_y_B, runavg_sigE_E_B
  REAL(rp) E_tot, gamma, KE, rbeta
  LOGICAL converged
  INTEGER counter, i

  REAL(rp) Vrf
  REAL(rp) U0
  REAL(rp) rf_freq, rf_omega

  CHARACTER(20) :: r_name = 'ibsequilibrium2'

  ! Used to plug in vertical beam size data
  !cur = lat%param%n_part*e_charge/(2.56E-6) * 1.0E3
  !(2.1 GeV Run 16) vbs = ( 0.057*cur**3 + 0.025*cur**2 - 0.552*cur + 21.833 ) * 1.0E-6
  !(1.8 GeV Run k) vbs = ( 0.130*cur**3 - 0.370*cur**2 + 0.997*cur + 22.655 ) * 1.0E-6
  !vbs = ( 0.363*cur**3 - 2.571*cur**2 + 5.668*cur + 30.854 ) * 1.0E-6
  !v_emit_data = vbs**2/40.0
  !ibsmode%b%emittance = v_emit_data

  Vrf = 0.0d0
  DO i=1,lat%n_ele_track
    IF(lat%ele(i)%key .eq. rfcavity$) THEN
      Vrf = Vrf + lat%ele(i)%value(voltage$)
      rf_freq = lat%ele(i)%value(rf_frequency$)
    ENDIF
  ENDDO
  rf_omega = twopi*rf_freq
  U0 = inmode%e_loss

  !compute the SR betatron damping times
  time_for_one_turn = lat%param%total_length / c_light
  tau_a = time_for_one_turn / inmode%a%alpha_damp
  tau_b = time_for_one_turn / inmode%b%alpha_damp
  tau_z = time_for_one_turn / inmode%z%alpha_damp 
  ibs_sim_params%tau_a = tau_a  !needed for tail cut calculation

  sigma_z0 = inmode%sig_z
  sigE_E0 = inmode%sigE_E 
  L_ratio = sigma_z0 / sigE_E0
  emit_a0 = inmode%a%emittance
  emit_b0 = inmode%b%emittance 

  cur = lat%param%n_part*e_charge/time_for_one_turn
  alpha_c = inmode%synch_int(1)/lat%param%total_length
  E_tot = lat%ele(0)%value(E_TOT$)

  ibsmode%a%emittance = emit_a0 * initial_blow_up
  ibsmode%b%emittance = emit_b0 * initial_blow_up
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
    DO i=1, lat%n_ele_track
      lat%ele(i)%a%emit = ibsmode%a%emittance
      lat%ele(i)%b%emit = ibsmode%b%emittance
      lat%ele(i)%z%sigma = ibsmode%sig_z
      lat%ele(i)%z%sigma_p = ibsmode%sigE_E
      lat%ele(i)%z%emit = ibsmode%sig_z * ibsmode%sigE_E
    ENDDO
    CALL ibs_rates1turn(lat,ibs_sim_params,rates,formula,granularity)
    counter = counter + 1
    !It is possible that this method can give negative emittances
    !at some point in the iterative process, in which case the case
    !structure below will terminate the program.  If this happens, try
    !using different values for initial_blow_up.
    IF( rates%inv_Ta .gt. 0.0 ) THEN
      Ta = 1.0/rates%inv_Ta
      afactor = 1.0/(1.0-(tau_a/Ta))
    ELSEIF( rates%inv_Ta .eq. 0.0 ) THEN
      afactor = 1.0
    ELSE
      CALL out_io(s_abort$, r_name, &
         'FATAL ERROR: Negative emittance encountered: ', &
         'Try adjusting initial_blow_up.')
      STOP
    ENDIF
    IF( rates%inv_Tb .gt. 0.0 ) THEN
      Tb = 1.0/rates%inv_Tb
      bfactor = 1.0/(1.0-(tau_b/Tb))
    ELSEIF( rates%inv_Tb .eq. 0.0 ) THEN
      bfactor = 1.0
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

    emit_a = ibsmode%a%emittance + advance*( afactor*emit_a0 - ibsmode%a%emittance ) 
    emit_b = ibsmode%b%emittance + advance*( &
            ((1.0-ratio)*bfactor+ratio*afactor)*emit_b0 - ibsmode%b%emittance )
    sigE_E = ibsmode%sigE_E      + advance*( zfactor*sigE_E0 - ibsmode%sigE_E )

    running_emit_x = CSHIFT(running_emit_x, -1)
    running_emit_y = CSHIFT(running_emit_y, -1)
    running_sigE_E = CSHIFT(running_sigE_E, -1)
    running_emit_x(1) = emit_a
    running_emit_y(1) = emit_b
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

    IF(.not.converged) THEN
      ibsmode%a%emittance = emit_a
      ibsmode%b%emittance = emit_b
      ibsmode%sigE_E = sigE_E 
      ! ibsmode%sig_z = L_ratio * sigE_E
      ! sigma_z_nat = L_ratio * sigE_E
      ! CALL potential_well_distortion(cur,alpha_c,E_tot,lat%z%tune/twopi,lat%param%total_length,inductance,sigma_z_nat,pwd_ratio)
      ! ibsmode%sig_z = sigma_z_nat * pwd_ratio
      CALL bl_via_vlassov(cur,alpha_c,E_tot,sigE_E,Vrf,rf_omega,U0,lat%param%total_length,ibs_sim_params%resistance, &
                          ibs_sim_params%inductance,sigma_z_vlassov)
      ibsmode%sig_z = sigma_z_vlassov
    ENDIF

  ENDDO

END SUBROUTINE ibs_equib_rlx

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
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
SUBROUTINE ibs_equib_der(lat,ibs_sim_params,inmode,ibsmode,formula,ratio,granularity)
  ! Iterates to equilibrium beam conditions using derivatives

  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  TYPE(normal_modes_struct), INTENT(IN) :: inmode
  TYPE(normal_modes_struct), INTENT(OUT) :: ibsmode
  REAL(rp), INTENT(IN) :: granularity
  CHARACTER*4, INTENT(IN) :: formula
  TYPE(ibs_struct) rates
  TYPE(normal_modes_struct) :: naturalmode

  REAL(rp) time_for_one_turn
  REAL(rp) tau_a, tau_b, tau_z
  REAL(rp) Ta, Tb, Tz
  REAL(rp) emit_a, emit_b, sigE_E
  REAL(rp) emit_a0, emit_b0, sigE_E0
  REAL(rp) dadt, dbdt, dsEdt, dT
  REAL(rp) threshold
  REAL(rp) sigma_z0, L_ratio
  REAL(rp) sigma_z_vlassov
  REAL(rp) ratio
  REAL(rp) cur
  REAL(rp) alpha_c
  REAL(rp) E_tot
  REAL(rp) sigma_z_nat
  REAL(rp) pwd_ratio

  REAL(rp) Vrf
  REAL(rp) U0
  REAL(rp) rf_freq, rf_omega

  LOGICAL converged
  INTEGER counter, i

  !natural mode here means emittances before IBS effects
  naturalmode = inmode

  !compute the SR betatron damping times
  time_for_one_turn = lat%param%total_length / c_light
  tau_a = time_for_one_turn / naturalmode%a%alpha_damp
  tau_b = time_for_one_turn / naturalmode%b%alpha_damp
  tau_z = time_for_one_turn / naturalmode%z%alpha_damp

  cur = lat%param%n_part*e_charge/time_for_one_turn
  alpha_c = inmode%synch_int(1)/lat%param%total_length
  E_tot = lat%ele(0)%value(E_TOT$)

  Vrf = 0.0d0
  DO i=1,lat%n_ele_track
    IF(lat%ele(i)%key .eq. rfcavity$) THEN
      Vrf = Vrf + lat%ele(i)%value(voltage$)
      rf_freq = lat%ele(i)%value(rf_frequency$)
    ENDIF
  ENDDO
  rf_omega = twopi*rf_freq
  U0 = inmode%e_loss

  sigma_z0 = naturalmode%sig_z
  sigE_E0 = naturalmode%sigE_E
  emit_a0 = naturalmode%a%emittance
  emit_b0 = naturalmode%b%emittance
  L_ratio = sigma_z0 / sigE_E0
  threshold = .00001 !fractional changes in emittance smaller than this
                     !indicate convergence
  converged = .false.
  !dT = tau_a / 240.0 !Time to advance per iteration
  dT = tau_a / 10.0 !Time to advance per iteration

  counter = 0
  ibsmode = naturalmode
  DO WHILE(.not.converged)
    DO i=1, lat%n_ele_track
      lat%ele(i)%a%emit = ibsmode%a%emittance
      lat%ele(i)%b%emit = ibsmode%b%emittance
      lat%ele(i)%z%sigma = ibsmode%sig_z
      lat%ele(i)%z%sigma_p = ibsmode%sigE_E
      lat%ele(i)%z%emit = ibsmode%sig_z * ibsmode%sigE_E
    ENDDO
    CALL ibs_rates1turn(lat,ibs_sim_params,rates,formula,granularity)
    counter = counter + 1
    Ta = 1.0/rates%inv_Ta
    Tb = 1.0/rates%inv_Tb
    Tz = 1.0/rates%inv_Tz
    emit_a = ibsmode%a%emittance
    emit_b = ibsmode%b%emittance
    sigE_E = ibsmode%sigE_E

    !Compute change in emittance per time for x,y,z dimensions, taking
    !into account radiation damping and IBS blow-up
    dadt = -(emit_a-emit_a0)*2./tau_a + emit_a*2./Ta
    dbdt = -(emit_b-emit_b0)*2./tau_b + emit_b*2./Tb 
    dsEdt= -(sigE_E-sigE_E0)/tau_z + sigE_E/Tz

    IF( (dadt*dT)/emit_a .lt. threshold ) THEN
      IF( (dbdt*dT)/emit_b .lt. threshold ) THEN
        IF( (dsEdt*dT)/ibsmode%sigE_E .lt. threshold ) THEN
          converged = .true.
        ENDIF
      ENDIF
    ENDIF
    IF(.not.converged) THEN
      ibsmode%a%emittance = emit_a + dadt*dT
      ibsmode%b%emittance = emit_b + dbdt*dT
      ibsmode%sigE_E = ibsmode%sigE_E + dsEdt*dT
      !sigma_z_nat = L_ratio * sigE_E
      !CALL potential_well_distortion(cur,alpha_c,E_tot,lat%z%tune/twopi,lat%param%total_length,inductance,sigma_z_nat,pwd_ratio)
      !ibsmode%sig_z = sigma_z_nat * pwd_ratio
      CALL bl_via_vlassov(cur,alpha_c,E_tot,ibsmode%sigE_E,Vrf,rf_omega,U0,lat%param%total_length, &
                          ibs_sim_params%resistance,ibs_sim_params%inductance,sigma_z_vlassov)
      ibsmode%sig_z = sigma_z_vlassov
    ENDIF
  ENDDO

!  WRITE(*,'(A,F11.3,F12.5)') "Final Coulomb Log: ", ibs_sim_params%clog, cur

END SUBROUTINE ibs_equib_der

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
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
SUBROUTINE ibs_lifetime(lat,tau_a,maxratio,lifetime,formula,granularity,clog_to_use)
  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  REAL(rp) tau_a
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  TYPE(ibs_maxratio_struct), INTENT(IN) :: maxratio
  TYPE(ibs_lifetime_struct), INTENT(OUT) :: lifetime
  REAL(rp), INTENT(IN) :: granularity
  CHARACTER*4, INTENT(IN) :: formula
  INTEGER clog_to_use

  TYPE(ibs_struct) rates

  REAL(rp) Rx, Ry, R_p

  Rx = maxratio%rx**2
  Ry = maxratio%ry**2
  R_p = maxratio%r_p**2

  ibs_sim_params%tau_a = tau_a
  ibs_sim_params%clog_to_use = clog_to_use
  CALL ibs_rates1turn(lat, ibs_sim_params, rates, formula, granularity)

  lifetime%Tlx = exp(Rx)/2/Rx/rates%inv_Ta
  lifetime%Tly = exp(Ry)/2/Ry/rates%inv_Tb
  lifetime%Tlp = exp(R_p)/2/R_p/rates%inv_Tz
END SUBROUTINE ibs_lifetime


!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
!  Subroutine ibs_delta_calc (lat, tau_a, ix, formula, clog_to_use, 
!                              delta_sigma_energy, delta_emit_a, delta_emit_b)
!
!  Calculates change in energy spread and emittances due to IBS.
!
!  Available IBS formulas:
!    cimp - Modified Piwinski
!    bane - Bane approximation of Bjorken-Mtingwa formulation
!    bjmt - Bjorken-Mtingwa formulation general to bunched beams (time consuming)
!    mtto - Mtingwa-Tollerstrup formulation
!
!  Input:
!    lat              -- lat_struct: lattice for tracking
!      %param%n_part  -- real(rp): number of particles in bunch
!    ix               -- integer: index of element to use: lat%ele(ix)
!    tau_a            -- real(rp) :: Horizontal damping rate. See See ibs_sim_param_struct.
!                                    Only used if clog_to_use == (Raubenheimer tail cut)
!    formula          -- character*4: IBS formulation to use
!    clog_to_use      -- integer: Coulomb log to use. See ibs_sim_param_struct
!
!
!  Output:
!    delta_sigma_energy -- real(rp), optional: change in energy spread in eV
!    delta_emit_a       -- real(rp), optional: change in a-mode emittance (geometric)
!    delta_emit_b       -- real(rp), optional: change in b-mode emittance (geometric)
!-
SUBROUTINE ibs_delta_calc (lat, ix, tau_a, formula, clog_to_use, delta_sigma_energy, delta_emit_a, delta_emit_b)
  IMPLICIT NONE

  TYPE(lat_struct), INTENT(IN) :: lat
  REAL(rp) tau_a
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  INTEGER, INTENT(IN) :: ix
  CHARACTER*4, INTENT(IN) ::  formula
  real(rp), optional :: delta_sigma_energy, delta_emit_a, delta_emit_b
  TYPE(ibs_struct) rates1ele
  INTEGER clog_to_use

  ibs_sim_params%tau_a = tau_a
  ibs_sim_params%clog_to_use = clog_to_use

  IF(lat%ele(ix)%value(l$) .gt. 0.0) THEN
    IF(formula == 'cimp') THEN
      CALL cimp1(lat, ibs_sim_params, rates1ele, i=ix)
    ELSEIF(formula == 'bjmt') THEN
      CALL bjmt1(lat, ibs_sim_params, rates1ele, i=ix)
    ELSEIF(formula == 'bane') THEN
      CALL bane1(lat, ibs_sim_params, rates1ele, i=ix)
    ELSEIF(formula == 'mpzt') THEN
      CALL mpzt1(lat, ibs_sim_params, rates1ele, i=ix)
    ELSEIF(formula == 'mpxx') THEN
      CALL mpxx1(lat, ibs_sim_params, rates1ele, i=ix)
    ELSE
      WRITE(*,*) "Invalid IBS formula selected ... returning zero"
      rates1ele%inv_Ta = 0.0
      rates1ele%inv_Tb = 0.0
      rates1ele%inv_Tz = 0.0
      RETURN
    ENDIF
  ENDIF

  if (present(delta_sigma_energy)) delta_sigma_energy = lat%ele(ix)%value(l$)/c_light*rates1ele%inv_Tz*lat%ele(ix)%z%sigma_p*lat%ele(ix)%value(E_TOT$)
  if (present(delta_emit_a))       delta_emit_a     = 2*lat%ele(ix)%value(l$)/c_light*rates1ele%inv_Ta*lat%ele(ix)%a%emit
  if (present(delta_emit_b))       delta_emit_b     = 2*lat%ele(ix)%value(l$)/c_light*rates1ele%inv_Tb*lat%ele(ix)%b%emit

END SUBROUTINE ibs_delta_calc

!+
!  Subroutine ibs_rates1turn(lat, rates1turn, formula, granularity)
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
SUBROUTINE ibs_rates1turn(lat, ibs_sim_params, rates1turn, formula, granularity)
  IMPLICIT NONE

  REAL(rp) sum_inv_Tz, sum_inv_Ta, sum_inv_Tb
  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  TYPE(ibs_struct) :: rates1ele
  TYPE(ibs_struct), INTENT(OUT) :: rates1turn
  REAL(rp), INTENT(IN) :: granularity
  CHARACTER*4, INTENT(IN) ::  formula
  INTEGER i
  INTEGER n_steps
  REAL(rp), ALLOCATABLE :: steps(:)
  REAL(rp) step_size
  REAL(rp) length_multiplier
!  REAL(rp) clog_avg !FOO


  sum_inv_Tz = 0.0
  sum_inv_Ta = 0.0
  sum_inv_Tb = 0.0

!  clog_avg = 0.0d0 !FOO

  IF( granularity .lt. 0.0 ) THEN
    DO i=1,lat%n_ele_track
      IF(lat%ele(i)%value(l$) .eq. 0.0) THEN
        CYCLE
      ENDIF
      IF(formula == 'cimp') THEN
        CALL cimp1(lat, ibs_sim_params, rates1ele, i=i)
      ELSEIF(formula == 'bjmt') THEN
        CALL bjmt1(lat, ibs_sim_params, rates1ele, i=i)
      ELSEIF(formula == 'bane') THEN
        CALL bane1(lat, ibs_sim_params, rates1ele, i=i)
      ELSEIF(formula == 'mpzt') THEN
        CALL mpzt1(lat, ibs_sim_params, rates1ele, i=i)
      ELSEIF(formula == 'mpxx') THEN
        CALL mpxx1(lat, ibs_sim_params, rates1ele, i=i)
      ELSEIF(formula == 'kubo') THEN
        CALL kubo1(lat, ibs_sim_params, rates1ele, ix=i)
      ELSE
        WRITE(*,*) "Invalid IBS formula selected ... returning zero"
        rates1turn%inv_Ta = 0.0
        rates1turn%inv_Tb = 0.0
        rates1turn%inv_Tz = 0.0
        RETURN
      ENDIF

      !length_multiplier = lat%ele(i)%value(l$)/2.0 + lat%ele(i+1)%value(l$)/2.0
      length_multiplier = lat%ele(i)%value(l$)
      sum_inv_Tz = sum_inv_Tz + rates1ele%inv_Tz * length_multiplier
      sum_inv_Ta = sum_inv_Ta + rates1ele%inv_Ta * length_multiplier
      sum_inv_Tb = sum_inv_Tb + rates1ele%inv_Tb * length_multiplier
!      clog_avg = clog_avg + ibs_sim_params%clog * length_multiplier !FOO
    ENDDO
  ELSE
    n_steps = CEILING(lat%param%total_length / granularity)
    step_size = lat%param%total_length / n_steps

    ALLOCATE(steps(1:n_steps))
    DO i=1,n_steps-1
      steps(i) = i*step_size
    ENDDO
    steps(n_steps) = lat%param%total_length

    DO i=1,n_steps
      IF(formula == 'cimp') THEN
        CALL cimp1(lat, ibs_sim_params, rates1ele, s=steps(i))
      ELSEIF(formula == 'bjmt') THEN
        CALL bjmt1(lat, ibs_sim_params, rates1ele, s=steps(i))
      ELSEIF(formula == 'bane') THEN
        CALL bane1(lat, ibs_sim_params, rates1ele, s=steps(i))
      ELSEIF(formula == 'mpzt') THEN
        CALL mpzt1(lat, ibs_sim_params, rates1ele, s=steps(i))
      ELSEIF(formula == 'mpxx') THEN
        CALL mpxx1(lat, ibs_sim_params, rates1ele, s=steps(i))
      ELSEIF(formula == 'kubo') THEN
        CALL kubo1(lat, ibs_sim_params, rates1ele, s=steps(i))
      ELSE
        WRITE(*,*) "Invalid IBS formula selected ... returning zero"
        rates1turn%inv_Ta = 0.0
        rates1turn%inv_Tb = 0.0
        rates1turn%inv_Tz = 0.0
        RETURN
      ENDIF

      sum_inv_Tz = sum_inv_Tz + rates1ele%inv_Tz * step_size
      sum_inv_Ta = sum_inv_Ta + rates1ele%inv_Ta * step_size
      sum_inv_Tb = sum_inv_Tb + rates1ele%inv_Tb * step_size
    ENDDO
    DEALLOCATE(steps)
  ENDIF

!  ibs_sim_params%clog = clog_avg / ibs_sim_params%lat%param%total_length !FOO

  rates1turn%inv_Tz = sum_inv_Tz / lat%param%total_length
  rates1turn%inv_Ta = sum_inv_Ta / lat%param%total_length
  rates1turn%inv_Tb = sum_inv_Tb / lat%param%total_length
END SUBROUTINE ibs_rates1turn

!+
!  Subroutine ibs_blowup1turn(lat,formula)
!
!  Available IBS formulas:
!    cimp - Modified Piwinski
!    bane - Bane approximation of Bjorken-Mtingwa formulation
!    bjmt - Bjorken-Mtingwa formulation general to bunched beams (time consuming)
!
!  Input:
!    lat             -- lat_struct: lattice for tracking
!      %param$n_part  -- Real: number of particles in bunch
!    formula          -- character*4: IBS formulation to use
!
!  Output:
!    rates          -- ibs_struct: ibs rates in x,y,z 
!-
SUBROUTINE ibs_blowup1turn(lat, ibs_sim_params, formula)
  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  CHARACTER*4, INTENT(IN) ::  formula
  TYPE(ibs_struct) :: rates1ele
  INTEGER i
  REAL(rp) delta_t, pp, gg, gamma1, gamma2

  DO i=1,lat%n_ele_track
    delta_t = lat%ele(i)%value(l$)/c_light
    IF(delta_t .gt. 0.) THEN
      rates1ele%inv_Ta = 0.0_rp
      rates1ele%inv_Tb = 0.0_rp
      rates1ele%inv_Tz = 0.0_rp

      IF(lat%ele(i)%value(l$) .ne. 0.0) THEN
        IF(formula == 'cimp') THEN
          CALL cimp1(lat, ibs_sim_params, rates1ele, i=i)
        ELSEIF(formula == 'bjmt') THEN
          CALL bjmt1(lat, ibs_sim_params, rates1ele, i=i)
        ELSEIF(formula == 'bane') THEN
          CALL bane1(lat, ibs_sim_params, rates1ele, i=i)
        ELSEIF(formula == 'mpzt') THEN
          CALL mpzt1(lat, ibs_sim_params, rates1ele, i=i)
        ELSEIF(formula == 'mpxx') THEN
          CALL mpxx1(lat, ibs_sim_params, rates1ele, i=i)
        ELSEIF(formula == 'kubo') THEN
          CALL kubo1(lat, ibs_sim_params, rates1ele, ix=i)
        ELSE
          WRITE(*,*) "Invalid IBS formula selected ... terminating"
          STOP
        ENDIF
      ENDIF

      CALL convert_total_energy_to(lat%ele(i)%value(E_TOT$), lat%param%particle, gamma1)
      CALL convert_total_energy_to(lat%ele(i+1)%value(E_TOT$), lat%param%particle, gamma2)
      gg = gamma1/gamma2
      pp = lat%ele(i)%value(p0c$)/lat%ele(i+1)%value(p0c$)
      lat%ele(i+1)%a%emit = lat%ele(i)%a%emit * (1 + delta_t*rates1ele%inv_Ta) * gg
      lat%ele(i+1)%b%emit = lat%ele(i)%b%emit * (1 + delta_t*rates1ele%inv_Tb) * gg
      lat%ele(i+1)%z%sigma_p = lat%ele(i)%z%sigma_p * (1 + delta_t*rates1ele%inv_Tz*2.0) * pp
    ELSE
      lat%ele(i+1)%a%emit = lat%ele(i)%a%emit 
      lat%ele(i+1)%b%emit = lat%ele(i)%b%emit 
      lat%ele(i+1)%z%sigma_p = lat%ele(i)%z%sigma_p 
    ENDIF
  ENDDO
END SUBROUTINE ibs_blowup1turn

!+
!  subroutine bjmt1(lat, ibs_sim_params, rates, i, s)
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
SUBROUTINE bjmt1(lat, ibs_sim_params, rates, i, s)

  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  INTEGER, INTENT(IN), OPTIONAL :: i
  REAL(rp), INTENT(IN), OPTIONAL :: s
  TYPE(ibs_struct), INTENT(OUT) :: rates
  TYPE(ele_struct), pointer :: ele
  TYPE(ele_struct), TARGET :: stubele

  REAL(rp) sigma_p, emit_a, emit_b, sigma_z, E_tot
  REAL(rp) gamma, KE, rbeta, beta_a, beta_b
  REAL(rp) sigma_y
  REAL(rp) Dx, Dy, Dxp, Dyp
  REAL(rp) alpha_a, alpha_b, coulomb_log
  REAL(rp) NB, big_A
  REAL(rp) Hx, Hy
  REAL(rp) inv_Tz, inv_Ta, inv_Tb

  REAL(rp) phi_h, phi_v

  REAL(c_double) :: Lp(3,3), Lh(3,3), Lv(3,3), L(3,3)
  REAL(c_double), TARGET :: mats(2,3,3)

  REAL(fgsl_double), TARGET :: Elpha
  TYPE(fgsl_function) :: integrand_ready
  REAL(fgsl_double) :: integration_result
  REAL(fgsl_double) :: abserr
  TYPE(c_ptr) :: ptr
  TYPE(fgsl_integration_workspace) :: integ_wk
  INTEGER(fgsl_int) :: fgsl_status

  NB = lat%param%n_part

  IF(PRESENT(i) .and. .not.PRESENT(s)) THEN
    ele => lat%ele(i)
  ELSEIF(PRESENT(s) .and. .not.PRESENT(i)) THEN
    CALL twiss_and_track_at_s(lat,s,stubele)
    ele => stubele
  ELSE
    WRITE(*,*) "FATAL ERROR IN ibs_mod: Either i or s (and not both) must be specified"
    STOP
  ENDIF

  E_tot = ele%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  Dx = ele%a%eta
  Dy = ele%b%eta
  Dxp = ele%a%etap
  Dyp = ele%b%etap
  sigma_y = SQRT(beta_b*emit_b + (Dy*sigma_p)**2)

  !coulomb_log = LOG( (gamma**2)*sigma_y*emit_a/r_e/beta_a )
  CALL multi_coulomb_log(lat, ibs_sim_params, ele, coulomb_log)

  Hx = ( Dx**2 + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
  Hy = ( Dy**2 + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b

  phi_h = Dxp + alpha_a*Dx/beta_a
  phi_v = Dyp + alpha_b*Dy/beta_b

  Lp = 0.0
  Lp(2,2) = 1.0
  Lp = (gamma**2)/(sigma_p**2)*Lp

  Lh = 0.0
  Lh(1,1) = 1.0
  Lh(1,2) = -1.0*gamma*phi_h
  Lh(2,1) = -1.0*gamma*phi_h
  Lh(2,2) = (gamma**2)*Hx/beta_a
  Lh = beta_a/emit_a*Lh

  Lv = 0.0
  Lv(2,2) = (gamma**2)*Hy/beta_b
  Lv(2,3) = -1.0*gamma*phi_v
  Lv(3,2) = -1.0*gamma*phi_v
  Lv(3,3) = 1.0
  Lv = beta_b/emit_b*Lv

  L = Lp + Lh + Lv

  mats(1,:,:) = L

  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(mats)
  integrand_ready = fgsl_function_init(bjmt_integrand, ptr)

  mats(2,:,:) = Lp
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 100.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  inv_Tz=4.0*pi*big_A*coulomb_log*integration_result

  mats(2,:,:) = Lh
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 100.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  inv_Ta=4.0*pi*big_A*coulomb_log*integration_result

  mats(2,:,:) = Lv
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 100.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  inv_Tb=4.0*pi*big_A*coulomb_log*integration_result

  CALL fgsl_integration_workspace_free(integ_wk)
  CALL fgsl_function_free(integrand_ready)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

END SUBROUTINE bjmt1

FUNCTION bjmt_integrand(xx,params) BIND(c)
  REAL(c_double), VALUE :: xx
  REAL(c_double) :: x
  TYPE(c_ptr), VALUE :: params
  REAL(c_double) :: bjmt_integrand

  REAL(c_double), POINTER :: mats(:,:,:)

  REAL(c_double) :: im(3,3)
  REAL(c_double) :: L(3,3)
  REAL(c_double) :: Li(3,3)
  REAL(c_double) :: Li_im(3,3)

  REAL(c_double) :: Tr_Li
  REAL(c_double) :: Tr_im
  REAL(c_double) :: Tr_Li_im
  REAL(c_double) :: Det_L_l

  CALL c_f_pointer(params,mats,[2,3,3])

  L = mats(1,:,:)
  Li = mats(2,:,:)

  x = EXP(xx) !change of variables

  Tr_Li = Li(1,1) + Li(2,2) + Li(3,3)

  im(1,1) = (L(2,2)+x)*(L(3,3)+x)-L(2,3)*L(3,2)
  im(1,2) = -L(1,2)*(L(3,3)+x)
  im(1,3) = L(1,2)*L(2,3)
  im(2,1) = -L(2,1)*(L(3,3)+x)
  im(2,2) = (L(1,1)+x)*(L(3,3)+x)
  im(2,3) = -(L(1,1)+x)*L(2,3)
  im(3,1) = L(2,3)*L(3,2)
  im(3,2) = -(L(1,1)+x)*L(3,2)
  im(3,3) = (L(1,1)+x)*(L(2,2)+x)-L(1,2)*L(2,1)

  im(:,:) = im(:,:) / ( (L(1,1)+x)*(L(2,2)+x)*(L(3,3)+x) - (L(1,1)+x)*L(3,2)*L(2,3) - L(1,2)*L(2,1)*(L(3,3)+x) )

  Tr_im = im(1,1)+im(2,2)+im(3,3)

  Li_im = matmul(Li,im)
  Tr_Li_im = Li_im(1,1)+Li_im(2,2)+Li_im(3,3)

  Det_L_l = (L(1,1)+x)*( (L(2,2)+x)*(L(3,3)+x) - L(3,2)*L(2,3) ) - L(1,2)*L(2,1)*(L(3,3)+x)

  bjmt_integrand = SQRT(x/Det_L_l) * ( Tr_Li*Tr_im - 3.0_rp*Tr_Li_im )
  bjmt_integrand = bjmt_integrand * x  !COV
END FUNCTION bjmt_integrand

!+
! This is a sigma matrix based IBS calculation.  It returns the growth rates of the normal mode emittances.
!-
SUBROUTINE kubo1(lat, ibs_sim_params, rates, ix, s)
  ! Some parts of this subroutine are patterned from the SAD accelerator code.

  USE mode3_mod
  USE eigen_mod

  IMPLICIT NONE

  ! Ps.sigma_mat rearranges the sigma matrix so that all xx terms are in top left block,
  ! all px terms are in top right block, and pp terms are in bottom left block.  Lower right
  ! block is px'.
  INTEGER, PARAMETER :: Ps(6,6) = RESHAPE( [1,0,0,0,0,0, 0,0,1,0,0,0, 0,0,0,0,1,0, &
                                            0,1,0,0,0,0, 0,0,0,1,0,0, 0,0,0,0,0,1],[6,6] )

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) ibs_sim_params
  TYPE(ibs_struct) :: rates
  INTEGER, OPTIONAL :: ix
  REAL(rp), OPTIONAL :: s

  REAL(rp) Tboost(6,6)
  REAL(rp) Tboost_inv(6,6)
  REAL(rp) Tspat(6,6)
  REAL(rp) Tspat_inv(6,6)
  REAL(rp) spatial_sigma_update(6,6)
  REAL(rp) sigma_update(6,6)
  REAL(rp) sigma_mat(6,6)
  REAL(rp) spatial_sigma_mat(6,6)
  REAL(rp) boosted_sigma_mat(6,6)
  REAL(rp) ar_sigma_mat(6,6)
  REAL(rp) sig_xx(3,3)
  REAL(rp) sig_xx_inv(3,3)
  REAL(rp) sig_xp(3,3)
  REAL(rp) sig_pp(3,3)
  REAL(rp) sig_pp_update(3,3)
  REAL(rp) sig_pl(3,3)
  REAL(rp) u(3)
  REAL(rp) R(3,3)
  INTEGER etypes(3)
  REAL(rp) g1, g2, g3
  REAL(rp) NB
  REAL(rp) clog
  REAL(rp) cI
  REAL(rp) Dw(3,3)
  REAL(rp) normal(3)
  REAL(rp) vol, vol1, pvol
  REAL(rp) ptrans, bn, bm, bmin, bmax
  REAL(rp) bmin1, bmin2
  REAL(rp) t6(6,6)
  REAL(rp) W(6,6)

  INTEGER error
  LOGICAL lerror
  INTEGER i, j, k
  REAL(rp) dt

  REAL(rp), PARAMETER :: pi_2 = pi/2.0d0

  TYPE(fgsl_integration_workspace) :: integ_wk
  TYPE(c_ptr) :: ptr
  TYPE(fgsl_function) :: integrand_ready
  INTEGER(fgsl_int) :: fgsl_status
  REAL(fgsl_double), TARGET :: args(3)
  REAL(fgsl_double) :: abserr
  REAL(fgsl_double) :: integration_result

  TYPE(normal_modes_struct) mode
  TYPE(ele_struct), POINTER :: ele
  REAL(rp) E_tot, gamma, KE, rbeta

  LOGICAL lerr, ok

  IF( ibs_sim_params%set_dispersion ) THEN
    W = 0.0_rp
    DO i=1,6
      W(i,i) = 1.0_rp
    ENDDO
    W(3,6) = -ibs_sim_params%eta_set
    W(4,6) = -ibs_sim_params%etap_set
    W(5,3) =  ibs_sim_params%etap_set
    W(5,4) = -ibs_sim_params%eta_set
  ENDIF

  NB = lat%param%n_part

  IF(PRESENT(ix) .and. .not.PRESENT(s)) THEN
    ele => lat%ele(ix)
  ELSEIF(PRESENT(s) .and. .not.PRESENT(ix)) THEN
    WRITE(*,*) "FAIL: Calculation at location s not currently supported by kubo1"
    STOP
    !CALL twiss_and_track_at_s(lat,s,stubele)
    !ele => stubele
  ELSE
    WRITE(*,*) "FATAL ERROR IN ibs_mod: Either ix or s (and not both) must be specified"
    STOP
  ENDIF

  E_tot = ele%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, rbeta)

  ! make sigma matrix in x,px,y,py,z,pz form
  mode%a%emittance = ele%a%emit
  mode%b%emittance = ele%b%emit
  mode%z%emittance = ele%z%emit
  CALL transfer_matrix_calc (lat, .true., t6, ix1=ix)
  IF( ibs_sim_params%set_dispersion ) THEN
    t6 = MATMUL(t6,W) 
  ENDIF
  CALL make_smat_from_abc(t6, mode, sigma_mat, lerr)
  IF( lerr ) THEN
    WRITE(*,'(A)') "BAD: make_smat_from_abc failed"
    rates%inv_Ta = 0.0d0
    rates%inv_Tb = 0.0d0
    rates%inv_Tz = 0.0d0
    RETURN
  ENDIF

  ! make transfer matrix from canonical to spatial coordinates
  Tspat = 0.0d0
  DO i=1,6
    Tspat(i,i) = 1.0d0
  ENDDO
!  Tspat(1,2) = ele%value(l$)/4
!  Tspat(3,4) = ele%value(l$)/4
!  Tspat(5,6) = ele%value(l$)/4/gamma/gamma
  spatial_sigma_mat = MATMUL(Tspat,MATMUL(sigma_mat,TRANSPOSE(Tspat)))

  ! boost sigma matrix to COM frame of bunch
  Tboost = 0.0d0
  do i=1,6
    Tboost(i,i) = 1.0d0
  enddo
  Tboost(5,5) = 1.0d0 + gamma*gamma/(1.0d0+gamma)
  Tboost(6,6) = 1.0d0 / gamma
  boosted_sigma_mat = MATMUL(Tboost,MATMUL(spatial_sigma_mat,TRANSPOSE(Tboost)))

  ! permute sigma matrix to x,y,z,px,py,pz form
  ar_sigma_mat = MATMUL(MATMUL(TRANSPOSE(Ps),boosted_sigma_mat),Ps)
  sig_xx = ar_sigma_mat(1:3,1:3)
  sig_xp = ar_sigma_mat(1:3,4:6)
  sig_pp = ar_sigma_mat(4:6,4:6)
  CALL eigensys(sig_xx, u, R, etypes, 3, error)
  vol1 = SQRT(u(1)*u(2)*u(3))
  vol = SQRT(4.0d0*pi)**3 * vol1
  bm = SQRT(MIN( u(1), u(2), u(3) ))  !minimum beam dimension
  CALL mat_inverse(sig_xx,sig_xx_inv,ok)
  IF( .not. ok ) THEN
    WRITE(*,*) "BAD: Could not invert sig_xx"
    rates%inv_Ta = 0.0d0
    rates%inv_Tb = 0.0d0
    rates%inv_Tz = 0.0d0
    RETURN
  ENDIF

  !Get local momentum matrix
  sig_pl = sig_pp - MATMUL(TRANSPOSE(sig_xp),MATMUL(sig_xx_inv,sig_xp))

  !Get eigen vectors of local momentum matrix
  CALL eigensys(sig_pl, u, R, etypes, 3, error)
  IF( error .ne. 0 ) THEN
    WRITE(*,'(A,I6," ",A)') "BAD: Eigenvectors of transfer matrix not found for element ", ix, ele%name
    rates%inv_Ta = 0.0d0
    rates%inv_Tb = 0.0d0
    rates%inv_Tz = 0.0d0
    RETURN
  ENDIF
  R=TRANSPOSE(R) !R is defined as inverse of eigenvector matrix, and tr(R) = inv(r)
  ptrans = SQRT(u(1)+u(2)+u(3))
  pvol = SQRT(u(1)*u(2)*u(3))
  
  !- Integration using fgsl
  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(kubo_integrand, ptr)

  args = (/u(1),u(2),u(3)/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, pi_2, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  g1 = integration_result

  args = (/u(2),u(1),u(3)/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, pi_2, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  g2 = integration_result

  args = (/u(3),u(1),u(2)/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, pi_2, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  g3 = integration_result

  CALL fgsl_integration_workspace_free(integ_wk)
  CALL fgsl_function_free(integrand_ready)
  !- End integration with fgsl

  bn = (vol/NB)**(1.0d0/3.0d0)
  bmax = MIN( bm, bn )  !minimum dimension or debye radius
  bmin1 = r_e/(ptrans*gamma)**2
  bmin2 = SQRT(ABS(vol/NB/pi/(ptrans*c_light)/ibs_sim_params%tau_a))
  bmin = MAX( bmin1, bmin2 )
  !bmin = bmin1 !FOO 
  clog = LOG(bmax/bmin)
!  ibs_sim_params%clog = clog  !smuggle out Coulomb Log so we can computes its average for curiosity' sake. !FOO

  !cI = (r_e**2)*NB*clog/4.0d0/pi/(gamma**4)/mode%a%emittance/mode%b%emittance/mode%z%emittance * ele%value(l$)
  cI = (r_e**2)*NB*clog/4.0d0/pi/(gamma**4)/vol1/pvol * ele%value(l$)

  Dw = 0.0d0
  Dw(1,1) = cI*(g2-g1+g3-g1)
  Dw(2,2) = cI*(g1-g2+g3-g2)
  Dw(3,3) = cI*(g1-g3+g2-g3)

  !- Build update matrix
  sig_pp_update = MATMUL(MATMUL(R,Dw),TRANSPOSE(R))
  sigma_update = 0.0d0
  DO i=1,3
    DO j=1,3
      sigma_update(i*2,j*2) = sig_pp_update(i,j)
    ENDDO
  ENDDO
  
  ! boost updates to lab frame
  CALL mat_inverse(Tboost,Tboost_inv,ok)
  IF( .not. ok ) THEN
    WRITE(*,*) "BAD: Could not invert Tboost"
    rates%inv_Ta = 0.0d0
    rates%inv_Tb = 0.0d0
    rates%inv_Tz = 0.0d0
    RETURN
  ENDIF
  spatial_sigma_update = MATMUL(Tboost_inv,MATMUL(sigma_update,TRANSPOSE(Tboost_inv)))

  ! back to canonical coordinates
  CALL mat_inverse(Tspat, Tspat_inv, ok)
  sigma_update = MATMUL(Tspat_inv,MATMUL(spatial_sigma_update,TRANSPOSE(Tspat_inv)))

  !- Update sigma matrix
  sigma_mat = sigma_mat + sigma_update

!  CALL get_abc_from_updated_smat(lat, ix, sigma_mat, normal, lerror)
  CALL normal_sigma_mat(sigma_mat,normal)

  !- Here we convert the updated sigma matrix into a time rate of change of the emittances.  This is done to
  !- make the kubo1 subroutine "drop in compatible" alongside the other IBS routines.
  dt = ele%value(l$) / c_light
  rates%inv_Ta = ((normal(1) - ele%a%emit) / dt)/ele%a%emit/2.0d0
  rates%inv_Tb = ((normal(2) - ele%b%emit) / dt)/ele%b%emit/2.0d0
  rates%inv_Tz = ((normal(3) - ele%z%emit) / dt)/ele%z%emit/2.0d0
END SUBROUTINE kubo1

FUNCTION kubo_integrand(s, params) BIND(c)
    REAL(c_double), VALUE :: s
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: kubo_integrand

    REAL(c_double), POINTER :: args(:)
    REAL(c_double) u1,u2,u3

    CALL c_f_pointer(params,args,[3])
    u1 = args(1)
    u2 = args(2)
    u3 = args(3)

    kubo_integrand = (2.0d0*u1*SIN(s)**2*COS(s)) / &
                     SQRT((SIN(s)**2 + u1/u2*COS(s)**2)*(SIN(s)**2+u1/u3*COS(s)**2))
END FUNCTION kubo_integrand

!+
!  subroutine bane1(lat, rates, i, s)
!
!  This is a private subroutine. To access this subroutine, call
!  ibs_rates.
!
!  This is an implementation of equations 10-15 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a high energy approximation of the Bjorken-Mtingwa IBS
!  formulation.
!-
SUBROUTINE bane1(lat, ibs_sim_params, rates, i, s)

  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  INTEGER, INTENT(IN), OPTIONAL :: i
  REAL(rp), INTENT(IN), OPTIONAL :: s
  TYPE(ibs_struct), INTENT(OUT) :: rates
  TYPE(ele_struct), pointer :: ele
  TYPE(ele_struct), TARGET :: stubele

  REAL(rp) sigma_p, emit_a, emit_b, sigma_z, E_tot
  REAL(rp) gamma, KE, rbeta, beta_a, beta_b
  REAL(rp) sigma_b, sigma_b_beta
  REAL(rp) Da, Db, Dap, Dbp
  REAL(rp) alpha_a, alpha_b, coulomb_log
  REAL(rp) a, b, g_bane
  REAL(rp) NB, big_A
  REAL(rp) sigma_H, Ha, Hb
  REAL(rp) inv_Tz, inv_Ta, inv_Tb

  REAL(fgsl_double), TARGET :: Elpha
  TYPE(fgsl_function) :: integrand_ready
  REAL(fgsl_double) :: integration_result
  REAL(fgsl_double) :: abserr
  TYPE(c_ptr) :: ptr
  TYPE(fgsl_integration_workspace) :: integ_wk
  INTEGER(fgsl_int) :: fgsl_status

  NB = lat%param%n_part

  IF(PRESENT(i) .and. .not.PRESENT(s)) THEN
    ele => lat%ele(i)
  ELSEIF(PRESENT(s) .and. .not.PRESENT(i)) THEN
    CALL twiss_and_track_at_s(lat,s,stubele)
    ele => stubele
  ELSE
    WRITE(*,*) "FATAL ERROR IN ibs_mod: Either i or s (and not both) must be specified"
    STOP
  ENDIF

  E_tot = ele%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*NB/16.0/(gamma**3)/(emit_a**(3./4.))/(emit_b**(3./4.))/sigma_z/(sigma_p**3)

  beta_a = ele%a%beta
  beta_b = ele%b%beta
  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  Dap = ele%a%etap
  Dbp = ele%b%etap
  Da = ele%a%eta
  Db = ele%b%eta
  sigma_b_beta = SQRT(beta_b * emit_b)
  sigma_b = SQRT(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  !coulomb_log = LOG( (gamma**2)*sigma_b*emit_a/r_e/beta_a )

  CALL multi_coulomb_log(lat, ibs_sim_params, ele, coulomb_log)

  Ha = ( (Da**2) + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( (Db**2) + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b
  sigma_H = 1.0/SQRT(1.0/(sigma_p**2) + Ha/emit_a + Hb/emit_b)

  a = sigma_H/gamma*SQRT(beta_a/emit_a)
  b = sigma_H/gamma*SQRT(beta_b/emit_b)

  Elpha = a/b
  ptr = c_loc(Elpha)
  integ_wk = fgsl_integration_workspace_alloc(limit)
  integrand_ready = fgsl_function_init(integrand, ptr)
  fgsl_status = fgsl_integration_qagiu(integrand_ready, 0.0d0, eps7, eps7, &
                                       limit, integ_wk, integration_result, abserr)
  g_bane = 2.0d0 * SQRT(Elpha) / pi * integration_result

  CALL fgsl_function_free(integrand_ready)
  CALL fgsl_integration_workspace_free(integ_wk)

  inv_Tz = big_A*coulomb_log*sigma_H*g_bane*((beta_a*beta_b)**(-1./4.))
  inv_Ta = (sigma_p**2)*Ha/emit_a*inv_Tz
  inv_Tb = (sigma_p**2)*Hb/emit_b*inv_Tz

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

END SUBROUTINE bane1

FUNCTION integrand(x, params) BIND(c)
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: integrand

    REAL(c_double), POINTER :: alpha
    CALL c_f_pointer(params, alpha)

    integrand = 1.0_c_double / ( SQRT(1+x*x)*SQRT(alpha*alpha+x*x) )
END FUNCTION integrand

!+
!  SUBROUTINE mpxx1(lat, rates, i, s)
!
!  This is a private subroutine. To access this subroutine, call
!  ibs_rates.
!
!-
SUBROUTINE mpxx1(lat, ibs_sim_params, rates, i, s)

  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  INTEGER, INTENT(IN), OPTIONAL :: i
  REAL(rp), INTENT(IN), OPTIONAL :: s
  TYPE(ibs_struct), INTENT(OUT), OPTIONAL :: rates
  TYPE(ele_struct), pointer :: ele
  TYPE(ele_struct), TARGET :: stubele

  REAL(rp) sigma_p, emit_a, emit_b, sigma_z, E_tot
  REAL(rp) gamma, KE, rbeta, beta_a, beta_b
  REAL(rp) sigma_a, sigma_b, sigma_a_beta, sigma_b_beta
  REAL(rp) Da, Db, Dap, Dbp
  REAL(rp) alpha_a, alpha_b, coulomb_log
  REAL(rp) a,b,q
  REAL(rp) NB, big_A
  REAL(rp) sigma_H, Ha, Hb
  REAL(rp) inv_Tz, inv_Ta, inv_Tb
  REAL(rp) fab, f1b, f1a

  REAL(fgsl_double), TARGET :: args(1:2)
  TYPE(fgsl_function) :: integrand_ready
  REAL(fgsl_double) :: integration_result
  REAL(fgsl_double) :: abserr
  TYPE(c_ptr) :: ptr
  TYPE(fgsl_integration_workspace) :: integ_wk
  INTEGER(fgsl_int) :: fgsl_status

  NB = lat%param%n_part

  IF(PRESENT(i) .and. .not.PRESENT(s)) THEN
    ele => lat%ele(i)
  ELSEIF(PRESENT(s) .and. .not.PRESENT(i)) THEN
    CALL twiss_and_track_at_s(lat,s,stubele)
    ele => stubele
  ELSE
    WRITE(*,*) "FATAL ERROR IN ibs_mod: Either i or s (and not both) must be specified"
    STOP
  ENDIF

  E_tot = ele%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_a_beta = SQRT(beta_a * emit_a)
  sigma_b_beta = SQRT(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = SQRT(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = SQRT(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( Db**2 + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b

  sigma_H = 1.0/SQRT( 1.0/(sigma_p**2)+ Ha/emit_a + Hb/emit_b )

  a = sigma_H/gamma*SQRT(beta_a/emit_a)
  b = sigma_H/gamma*SQRT(beta_b/emit_b)

  !------------------------Begin calls to GSL integrator
  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(mpxx_integrand, ptr)

  args = (/a,b/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  fab = integration_result

  args = (/ 1.0_rp/a, b/a /) 
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1b = integration_result

  args = (/ 1.0_rp/b, a/b /)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1a = integration_result

  CALL fgsl_integration_workspace_free(integ_wk)
  CALL fgsl_function_free(integrand_ready)
  !------------------------End calls to GSL integrator

  CALL multi_coulomb_log(lat, ibs_sim_params, ele, coulomb_log)

  inv_Tz = coulomb_log * big_A * sigma_H**2 / sigma_p**2 * fab
  inv_Ta = coulomb_log * big_A * (f1b + Ha*sigma_H**2/emit_a*fab)
  inv_Tb = coulomb_log * big_A * (f1a + Hb*sigma_H**2/emit_b*fab)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb
END SUBROUTINE mpxx1

FUNCTION mpxx_integrand(x, params) BIND(c)
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: mpxx_integrand

    REAL(c_double), POINTER :: args(:)
    REAL(c_double) av, bv

    REAL(c_double) u

    CALL c_f_pointer(params,args,[2])
    av = args(1)
    bv = args(2)

    mpxx_integrand = 8.0d0*pi*(1-3.0d0*x*x) / &
                     SQRT(av*av+(1-av*av)*x*x) / &
                     SQRT(bv*bv+(1-bv*bv)*x*x)

END FUNCTION mpxx_integrand

!+
!  SUBROUTINE mpzt1(lat, rates, i, s)
!
!  This is a private subroutine. To access this subroutine, call
!  ibs_rates.
!
!-
SUBROUTINE mpzt1(lat, ibs_sim_params, rates, i, s)

  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  INTEGER, INTENT(IN), OPTIONAL :: i
  REAL(rp), INTENT(IN), OPTIONAL :: s
  TYPE(ibs_struct), INTENT(OUT), OPTIONAL :: rates
  TYPE(ele_struct), pointer :: ele
  TYPE(ele_struct), TARGET :: stubele

  REAL(rp) sigma_p, emit_a, emit_b, sigma_z, E_tot
  REAL(rp) gamma, KE, rbeta, beta_a, beta_b
  REAL(rp) sigma_a, sigma_b, sigma_a_beta, sigma_b_beta
  REAL(rp) Da, Db, Dap, Dbp
  REAL(rp) alpha_a, alpha_b, coulomb_log
  REAL(rp) a,b,q
  REAL(rp) NB, big_A
  REAL(rp) sigma_H, Ha, Hb
  REAL(rp) inv_Tz, inv_Ta, inv_Tb
  REAL(rp) fabq, f1bq, f1aq

  REAL(fgsl_double), TARGET :: args(1:3)
  TYPE(fgsl_function) :: integrand_ready
  REAL(fgsl_double) :: integration_result
  REAL(fgsl_double) :: abserr
  TYPE(c_ptr) :: ptr
  TYPE(fgsl_integration_workspace) :: integ_wk
  INTEGER(fgsl_int) :: fgsl_status

  NB = lat%param%n_part

  IF(PRESENT(i) .and. .not.PRESENT(s)) THEN
    ele => lat%ele(i)
  ELSEIF(PRESENT(s) .and. .not.PRESENT(i)) THEN
    CALL twiss_and_track_at_s(lat,s,stubele)
    ele => stubele
  ELSE
    WRITE(*,*) "FATAL ERROR IN ibs_mod: Either i or s (and not both) must be specified"
    STOP
  ENDIF

  E_tot = ele%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_a_beta = SQRT(beta_a * emit_a)
  sigma_b_beta = SQRT(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = SQRT(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = SQRT(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( Db**2 + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b

  sigma_H = 1.0/SQRT( 1.0/(sigma_p**2)+ Ha/emit_a + Hb/emit_b )

  a = sigma_H/gamma*SQRT(beta_a/emit_a)
  b = sigma_H/gamma*SQRT(beta_b/emit_b)
  q = sigma_H*rbeta*SQRT(2.0_rp*sigma_b/r_e)
  !---- q = (gamma**2)*sigma_b*emit_a/r_e/beta_a  !effective coulomb log 

  !------------------------Begin calls to GSL integrator
  integ_wk = fgsl_integration_workspace_alloc(limit)
  ptr = c_loc(args)
  integrand_ready = fgsl_function_init(zot_integrand, ptr)

  args = (/a,b,q/)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  fabq = integration_result

  args = (/ 1.0_rp/a, b/a, q/a /) 
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1bq = integration_result

  args = (/ 1.0_rp/b, a/b, q/b /)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 1.0d0, eps7, eps7, &
                                     limit, 3, integ_wk, integration_result, abserr)
  f1aq = integration_result

  CALL fgsl_integration_workspace_free(integ_wk)
  CALL fgsl_function_free(integrand_ready)
  !------------------------End calls to GSL integrator

  inv_Tz = big_A * sigma_H**2 / sigma_p**2 * fabq
  inv_Ta = big_A * (f1bq + Ha*sigma_H**2/emit_a*fabq)
  inv_Tb = big_A * (f1aq + Hb*sigma_H**2/emit_b*fabq)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb
END SUBROUTINE mpzt1

FUNCTION zot_integrand(x, params) BIND(c)
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: zot_integrand

    REAL(c_double), POINTER :: args(:)
    REAL(c_double) av, bv, qv, P, Q

    REAL(c_double) u

    CALL c_f_pointer(params,args,[3])
    av = args(1)
    bv = args(2)
    qv = args(3)

    !COV to remove singularity at endpoints
    u=30.*((x**5.)/5. - (x**4.)/2. + (x**3.)/3.)

    P = SQRT(av*av + (1.-av*av)*u*u)
    Q = SQRT(bv*bv + (1.-bv*bv)*u*u)
    zot_integrand = 8.0_rp*pi * (1.0_rp-3.0_rp*u*u)/P/Q * &
                (2.0_rp*LOG(qv/2.0_rp*(1/P+1/Q)) - 0.577215665)

    zot_integrand = zot_integrand * 30.*(x**4. - 2.*(x**3.) + x**2.)  !COV to remove singularity

    !without COV
    !P = SQRT(av*av + (1-av*av)*x*x)
    !Q = SQRT(bv*bv + (1-bv*bv)*x*x)
    !zot_integrand = 8.0_rp*pi * (1.0_rp-3.0_rp*x*x)/P/Q * &
    !            (2.0_rp*LOG(qv/2.0_rp*(1/P+1/Q)) - 0.577215665)

END FUNCTION zot_integrand

!+
!  SUBROUTINE cimp1(lat, rates, i, s)
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
SUBROUTINE cimp1(lat, ibs_sim_params, rates, i, s)

  IMPLICIT NONE

  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  INTEGER, INTENT(IN), OPTIONAL :: i
  REAL(rp), INTENT(IN), OPTIONAL :: s
  TYPE(ibs_struct), INTENT(OUT), OPTIONAL :: rates
  TYPE(ele_struct), pointer :: ele
  TYPE(ele_struct), TARGET :: stubele

  REAL(rp) sigma_p, emit_a, emit_b, sigma_z, E_tot
  REAL(rp) gamma, KE, rbeta, beta_a, beta_b
  REAL(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
  REAL(rp) Dx, Dy, Dxp, Dyp
  REAL(rp) alpha_a, alpha_b, coulomb_log
  REAL(rp) a, b
  REAL(rp) NB, big_A
  REAL(rp) sigma_H, Hx, Hy
  REAL(rp) inv_Tz, inv_Ta, inv_Tb
  REAL(rp) g_ab,g_ba
  REAL(rp) bminstar, bmax

  !- Code specific to alternitive log representation
  REAL(rp) q, lnqa, lnqb
  !-

  NB = lat%param%n_part

  IF(PRESENT(i) .and. .not.PRESENT(s)) THEN
    ele => lat%ele(i)
  ELSEIF(PRESENT(s) .and. .not.PRESENT(i)) THEN
    CALL twiss_and_track_at_s(lat,s,stubele)
    ele => stubele
  ELSE
    WRITE(*,*) "FATAL ERROR IN ibs_mod: Either i or s (and not both) must be specified"
    STOP
  ENDIF

  E_tot = ele%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*NB/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_x_beta = SQRT(beta_a * emit_a)
  sigma_y_beta = SQRT(beta_b * emit_b)
  Dx =ele%a%eta
  Dy = ele%b%eta
  Dxp = ele%a%etap
  Dyp = ele%b%etap
  sigma_x = SQRT(sigma_x_beta**2 + (Dx**2)*(sigma_p**2))
  sigma_y = SQRT(sigma_y_beta**2 + (Dy**2)*(sigma_p**2))

  Hx = ( Dx**2 + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
  Hy = ( Dy**2 + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b

  sigma_H = 1.0/SQRT( 1.0/(sigma_p**2)+ Hx/emit_a + Hy/emit_b )

  CALL multi_coulomb_log(lat, ibs_sim_params, ele, coulomb_log)

  a = sigma_H/gamma*SQRT(beta_a/emit_a)
  b = sigma_H/gamma*SQRT(beta_b/emit_b)

  g_ba = g(b/a)
  g_ab = g(a/b)

  inv_Tz = 2.*(pi**(3./2.))*big_A*(sigma_H**2)/(sigma_p**2) * &
    coulomb_log * ( g_ba/a + g_ab/b ) 
  inv_Ta = 2.*(pi**(3./2.))*big_A*coulomb_log*&
    (-a*g_ba + Hx*(sigma_H**2)/emit_a* &
    ( g_ba/a + g_ab/b ) ) 
  inv_Tb = 2.*(pi**(3./2.))*big_A*coulomb_log*&
    (-b*g_ab + Hy*(sigma_H**2)/emit_b* &
    ( g_ba/a + g_ab/b ) )

  !!- Code specific to alternitive log representation
  !inv_Tz = 2.*(pi**(3./2.))*big_A*(sigma_H**2)/(sigma_p**2) * &
  !  ( lnqa*g_ba/a + lnqb*g_ab/b ) 
  !inv_Ta = 2.*(pi**(3./2.))*big_A* &
  !  (-a*lnqa*g_ba + Hx*(sigma_H**2)/emit_a* &
  !  ( lnqa*g_ba/a + lnqb*g_ab/b ) ) 
  !inv_Tb = 2.*(pi**(3./2.))*big_A* &
  !  (-b*lnqb*g_ab + Hy*(sigma_H**2)/emit_b* &
  !  ( lnqa*g_ba/a + lnqb*g_ab/b ) )
  !!-

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb
END SUBROUTINE cimp1

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
!-
SUBROUTINE multi_coulomb_log(lat,ibs_sim_params,ele,coulomb_log)
  TYPE(lat_struct) :: lat
  TYPE(ibs_sim_param_struct) :: ibs_sim_params
  REAL(rp) coulomb_log

  TYPE(ele_struct) ele
  REAL(rp) gamma, E_tot, g2
  REAL(rp) sigma_a, sigma_b, sigma_z, sigma_p, sp2
  REAL(rp) sigma_a_beta, sigma_b_beta
  REAL(rp) emit_a, emit_b
  REAL(rp) beta_a, beta_b
  REAL(rp) alpha_a, alpha_b
  REAL(rp) Da, Db, Dap, Dbp
  REAL(rp) NB
  REAL(rp) bminstar
  REAL(rp) bmax
  REAL(rp) gamma_a, gamma_b
  REAL(rp) Ha, Hb
  REAL(rp) Bbar
  REAL(rp) u, v, w  !used for Raubenheimer's calculation
  REAL(rp) qmin, qmax
  REAL(rp) vol, debye

  !fgsl variables
  TYPE(c_ptr) :: ptr
  REAL(fgsl_double), TARGET :: args(1:3)
  TYPE(fgsl_integration_workspace) :: integ_wk
  TYPE(fgsl_function) :: integrand_ready
  INTEGER(fgsl_int) :: fgsl_status
  REAL(fgsl_double) :: integration_result
  REAL(fgsl_double) :: abserr

  NB = lat%param%n_part

  E_tot = ele%value(E_TOT$)
  CALL convert_total_energy_to(E_tot, lat%param%particle, gamma=gamma)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  beta_a = ele%a%beta
  beta_b = ele%b%beta
  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  gamma_a = ele%a%gamma
  gamma_b = ele%b%gamma
  sigma_a_beta = SQRT(beta_a * emit_a)
  sigma_b_beta = SQRT(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = SQRT(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = SQRT(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  IF( ibs_sim_params%clog_to_use == 1 ) THEN
    !Classic Coulomb Log.
    coulomb_log = LOG( (gamma**2)*sigma_b*emit_a/r_e/beta_a )
  ELSEIF( ibs_sim_params%clog_to_use == 2 ) THEN
    !Tail cut Raubenheimer gstar (ignores vertical dispersion !)
    Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
    g2 = gamma*gamma
    sp2 = sigma_p*sigma_p
    u = g2*( Ha/emit_a + 1.0/sp2 + beta_a/g2/emit_a + beta_b/g2/emit_b )
    v = g2*( Ha*beta_b/emit_a/emit_b + beta_a*beta_b/g2/emit_a/emit_b + beta_a/sp2/emit_a + beta_b/sp2/emit_b + Da*Da/emit_a/emit_a ) 
    w = g2*( Da*Da*beta_b/emit_a/emit_a/emit_b + beta_a*beta_b/sp2/emit_a/emit_b )

    Bbar = gamma*SQRT( gamma_a*emit_a + gamma_b*emit_b + sp2/g2)
    qmin = 2.0 * r_e / sigma_b / Bbar

    !FGSL integration
    ptr = c_loc(args)
    args = (/u,v,w/)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    integrand_ready = fgsl_function_init(rclog_integrand, ptr)
    fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 40.0d0, eps7, eps7, &
                                           limit, 3, integ_wk, integration_result, abserr)

    CALL fgsl_function_free(integrand_ready)
    CALL fgsl_integration_workspace_free(integ_wk)

    qmax = SQRT( ibs_sim_params%tau_a*NB*c_light*r_e*r_e/4.0/pi/g2/emit_a/emit_b/sigma_z/sigma_p * integration_result )

    coulomb_log = LOG(qmax/qmin)
  ELSEIF( ibs_sim_params%clog_to_use == 3 ) THEN
    !Tail cut Bane form: SLAC-PUB-9227
    vol = (4.0d0*pi)**(3.0d0/2.0d0)*sigma_a*sigma_b*sigma_z*gamma
    bminstar = SQRT(vol/NB/pi/(ibs_sim_params%tau_a/gamma)) / SQRT(c_light*gamma*SQRT(emit_a/beta_a))
    debye = (vol/NB)**(1.0d0/3.0d0)
    bmax = MIN(MIN(sigma_a,sigma_b),debye)
    coulomb_log = LOG(bmax/bminstar)
  ELSEIF( ibs_sim_params%clog_to_use == 4) THEN
    !Tail cut Bane, but using relative velocity in all 3 dims, rather than just horizontal
    vol = (4.0d0*pi)**(3.0d0/2.0d0)*sigma_a*sigma_b*sigma_z*gamma
    Bbar = SQRT( gamma_a*emit_a + gamma_b*emit_b + sigma_p*sigma_p/gamma/gamma)  !beta in the rest frame
    bminstar = SQRT(vol/NB/pi/(ibs_sim_params%tau_a/gamma)) / SQRT(c_light*gamma*Bbar)
    debye = (vol/NB)**(1.0d0/3.0d0)
    bmax = MIN(MIN(sigma_a,sigma_b),debye)
    coulomb_log = LOG(bmax/bminstar)
  ENDIF

END SUBROUTINE multi_coulomb_log

FUNCTION rclog_integrand(x, params) BIND(c)
    REAL(c_double), VALUE :: x
    TYPE(c_ptr), VALUE :: params
    REAL(c_double) :: rclog_integrand

    REAL(c_double), POINTER :: args(:)
    REAL(c_double) u, v, w

    CALL c_f_pointer(params, args, [3])
    u = args(1)
    v = args(2)
    w = args(3)

    rclog_integrand = 2.0_rp/SQRT(EXP(6.0_rp*x) + u*EXP(4.0_rp*x) + v*EXP(2.0_rp*x) + w)*EXP(x)
END FUNCTION rclog_integrand

!+
!-
SUBROUTINE bl_via_vlassov(current,alpha,Energy,sigma_p,Vrf,omega,U0,circ,R,L,sigma_z)
  USE longitudinal_profile_mod

  IMPLICIT none

  REAL(rp) current
  REAL(rp) alpha
  REAL(rp) Energy
  REAL(rp) sigma_p
  REAL(rp) Vrf
  REAL(rp) omega
  REAL(rp) U0
  REAL(rp) circ
  REAL(rp) R
  REAL(rp) L
  REAL(rp) sigma_z

  REAL(rp) delta_e
  REAL(rp) A
  REAL(rp) Q
  REAL(rp) T0
  REAL(rp) phi
  REAL(rp), PARAMETER :: bound = 500.0d-12

  REAL(rp) args(1:8)

  delta_e = sigma_p * Energy
  T0 = circ/c_light
  Q = current * T0
  phi = -ACOS(U0/Vrf)
  A = Energy/delta_e/delta_e/alpha/T0

  args(1) = A
  args(2) = Vrf
  args(3) = Q
  args(4) = omega
  args(5) = phi
  args(6) = R
  args(7) = L
  args(8) = U0

  CALL get_bl_from_fwhm(bound,args,sigma_z)
END SUBROUTINE bl_via_vlassov

!+
!-
!SUBROUTINE potential_well_distortion(I,alpha,Energy,tune,circ,inductance,sig_nat,ratio)
!  ! Computes the ratio between energy spread and bunch length, taking potential well distortion into account.
!  real(rp) i             ! current
!  real(rp) alpha         ! momentum compaction
!  real(rp) energy        ! beam energy in ev
!  real(rp) tune          ! fractional part of synchrotron tune
!  real(rp) circ          ! circumpherence of storage ring
!  real(rp) inductance    ! longitudinal inductance of maching (typically 10s of nh)
!  real(rp) sig_nat       ! bunch length in absence of pwd
!  real(rp) ratio         ! use to calculate bunch length:  sigma_pwd(i) = sigma_nat * ratio(i)
!
!  REAL(rp) :: onethird = 0.333333333333333333333_rp
!  REAL(rp) :: twothird = 0.666666666666666666666_rp
!  REAL(rp) :: roottwopi = 2.506628274631_rp
!  COMPLEX(rp) :: pioverthree = (0.5_rp,0.8660254037844_rp)
!  REAL(rp) Zn
!  REAL(rp) R
!  REAL(rp) A
!  REAL(rp) term
!  COMPLEX(rp) arg, argcubed
!
!  R = circ / twopi
!  Zn = twopi/(circ/c_light) * inductance
!  A = alpha/Energy/tune/tune*((R/sig_nat)**3)*Zn/roottwopi
!
!  term = 9.0_rp*A*I
!  argcubed = SQRT(CMPLX(term*term-12.0)) - CMPLX(term)
!  arg = argcubed**onethird
!
!  ratio = -1.0*( (-2.,0)*pioverthree*(3.,0)**onethird + pioverthree**2*(2.,0)**onethird*arg*arg ) / (6.,0)**twothird / arg
!
!END SUBROUTINE potential_well_distortion

END MODULE ibs_mod






