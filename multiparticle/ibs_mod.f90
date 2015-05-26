module ibs_mod

use bmad
use fgsl
use, intrinsic :: iso_c_binding

type ibs_sim_param_struct
  real(rp) tau_a           ! horizontal damping rate (needed for coulomb log tail cut)
  integer clog_to_use      ! see multi_coulomb_log subroutine for valid settings.  Set to 1 to disable tail-cut.  Set to 1 for linacs.
  logical set_dispersion   ! True: add vertical dispersion to transfer matrix.  Valid for kubo method.
  real(rp) eta_set         ! If set_dispersion, then this value is used to add y-z coupling to the transfer matrix.
  real(rp) etap_set        ! If set_dispersion, then this value is used to add y-z coupling to the transfer matrix.
  logical do_pwd           ! If true, then use potential well distortion to calculate bunch lengths.  If false,
                           ! bunch length is proportional to energy spread.
  real(rp) inductance      ! Inductive part of impedance for pwd calc.  
  real(rp) resistance      ! Resistive part of impedance for pwd calc.
  character(4) formula     ! Which IBS formulation to use.  See subroutine ibs1 for a list.
  logical use_t6_cache     ! use ele%r(:,:,1) for one turn mat, rather than calculating fresh every turn.  Relevant only for kubo. 
end type

type ibs_struct  !these are betatron growth rates.  To get emittance growth rate use:
                 ! demit_a/dt = 2.0*emit_*inv_Ta
  real(rp) inv_Ta
  real(rp) inv_Tb
  real(rp) inv_Tz
end type

type ibs_lifetime_struct  ! Beam lifetime based on IBS.  Useful for linacs.  These quantities are populated with time
                          ! required for beam size to increase by come amount.
  real(rp) Tlx
  real(rp) Tly
  real(rp) Tlp
end type

type ibs_maxratio_struct  ! Parameters for IBS lifetime calculation
  real(rp) rx
  real(rp) ry
  real(rp) r_p
end type

real(fgsl_double), parameter :: eps7 = 1.0d-7
integer(fgsl_size_t), parameter :: limit = 1000_fgsl_size_t

contains

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  Subroutine ibs_equib_rlx(lat,ibs_sim_params,inmode,ibsmode,ratio,initial_blow_up,granularity)
!  Iterates to equilibrium beam conditions using relaxation method
!
!  This method requires that the initial beam size be larger than the equilibrium beam size.
!  An initial_blow_up of 3 to 5 is a good place to start.
!
!  See ibs_rates subroutine for available IBS rate formulas.
!
!  Modules needed:
!    use ibs_mod
!
!  Input:
!    lat             -- lat_struct: lattice for tracking
!      %param%n_part  -- Real: number of particles in bunch
!    ibs_sim_params   -- ibs_sim_param_struct: parameters for IBS calculation
!    inmode           -- normal_modes_struct: natural beam parameters 
!    ratio            -- Real: Ratio of vert_emit_coupling / vert_emit_total
!    initial_blow_up  -- Real: Factor multiplied to all thre bunch dimensions prior to starting iteration.
!    granularity      -- Real: Step size for slicing lattice.  i.e. set to 1 to calculate IBS rates every 1 meter.
!
!  Output:
!    ibsmode          -- normal_modes_struct: beam parameters after IBS effects
!-

subroutine ibs_equib_rlx(lat,ibs_sim_params,inmode,ibsmode,ratio,initial_blow_up,granularity)

  implicit none
                                                                                                         
  type(lat_struct) :: lat
  type(ibs_sim_param_struct) ibs_sim_params
  type(normal_modes_struct) :: inmode
  type(normal_modes_struct) :: ibsmode
  real(rp), intent(in) :: granularity
  real(rp) :: ratio
  real(rp) :: initial_blow_up
  type(ibs_struct) rates

  real(rp) time_for_one_turn
  real(rp) cur
  real(rp) alpha_c
  real(rp) pwd_ratio
  real(rp) sigma_z_nat
  real(rp) sigma_z_vlassov
  real(rp) tau_a, tau_b, tau_z
  real(rp) Ta, Tb, Tz
  real(rp) afactor,bfactor,zfactor
  real(rp) emit_a, emit_b
  real(rp) emit_a0, emit_b0
  real(rp) advance, threshold
  real(rp) sigE_E, sigE_E0, sigma_z0, L_ratio
  real(rp) running_emit_x(1:40), running_emit_y(1:40), running_sigE_E(1:40)
  integer half
  real(rp) runavg_emit_x_A, runavg_emit_y_A, runavg_sigE_E_A
  real(rp) runavg_emit_x_B, runavg_emit_y_B, runavg_sigE_E_B
  real(rp) energy, gamma, KE, rbeta
  logical converged
  integer counter, i

  real(rp) Vrf
  real(rp) U0
  real(rp) rf_freq, rf_omega
  real(rp) sig_z

  character(20) :: r_name = 'ibs_equib_rlx'

  ! Used to plug in vertical beam size data
  !cur = lat%param%n_part*e_charge/(2.56E-6) * 1.0E3
  !(2.1 GeV Run 16) vbs = ( 0.057*cur**3 + 0.025*cur**2 - 0.552*cur + 21.833 ) * 1.0E-6
  !(1.8 GeV Run k) vbs = ( 0.130*cur**3 - 0.370*cur**2 + 0.997*cur + 22.655 ) * 1.0E-6
  !vbs = ( 0.363*cur**3 - 2.571*cur**2 + 5.668*cur + 30.854 ) * 1.0E-6
  !v_emit_data = vbs**2/40.0
  !ibsmode%b%emittance = v_emit_data

  Vrf = 0.0d0
  do i=1,lat%n_ele_track
    if(lat%ele(i)%key .eq. rfcavity$) then
      Vrf = Vrf + lat%ele(i)%value(voltage$)
      rf_freq = lat%ele(i)%value(rf_frequency$)
    endif
  enddo
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
  energy = lat%ele(0)%value(E_TOT$)

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
  advance = .02
  !Changes in emittance between iterations less than threshold
  !indicate convergence.
  threshold = advance * 0.001

  do i=1, SIZE(running_emit_x)
    running_emit_x(i) = 0.0
    running_emit_y(i) = 0.0
    running_sigE_E(i) = 0.0
  enddo
  half = SIZE(running_emit_x) / 2

  do WHILE(.not.converged)
    if( ibs_sim_params%do_pwd ) then
      ! call bl_via_vlassov(cur,alpha_c,energy,ibsmode%sigE_E,Vrf,rf_omega,U0,lat%param%total_length, &
      !                     ibs_sim_params%resistance,ibs_sim_params%inductance,sigma_z_vlassov)
      call bl_via_mat(lat, ibs_sim_params, ibsmode, sig_z)
      ibsmode%sig_z = sig_z
    else
      ibsmode%sig_z = L_ratio * ibsmode%sigE_E
    endif

    do i=1, lat%n_ele_track
      lat%ele(i)%a%emit = ibsmode%a%emittance
      lat%ele(i)%b%emit = ibsmode%b%emittance
      lat%ele(i)%z%sigma = ibsmode%sig_z
      lat%ele(i)%z%sigma_p = ibsmode%sigE_E
      lat%ele(i)%z%emit = ibsmode%sig_z * ibsmode%sigE_E
    enddo

    call ibs_rates1turn(lat,ibs_sim_params,rates,granularity)

    counter = counter + 1
    !It is possible that this method can give negative emittances
    !at some point in the iterative process, in which case the case
    !structure below will terminate the program.  If this happens, try
    !using different values for initial_blow_up.
    if( rates%inv_Ta .eq. 0.0 ) then
      afactor = 1.0
    else
      Ta = 1.0/rates%inv_Ta
      afactor = 1.0/(1.0-(tau_a/Ta))
      if( afactor .lt. 0.0 ) then
        call out_io(s_abort$, r_name, &
             'FATAL ERROR: Negative emittance encountered: ', &
             'Try adjusting initial_blow_up or switch to derivatives method "der".')
        STOP
      endif
    endif
    if( rates%inv_Tb .eq. 0.0 ) then
      bfactor = 1.0
    else
      Tb = 1.0/rates%inv_Tb
      bfactor = 1.0/(1.0-(tau_b/Tb))
      if( bfactor .lt. 0.0 ) then
        call out_io(s_abort$, r_name, &
             'FATAL ERROR: Negative emittance encountered: ', &
             'Try adjusting initial_blow_up or switch to derivatives method "der".')
        STOP
      endif
    endif
    if( rates%inv_Tz .eq. 0.0 ) then
      zfactor = 1.0
    else
      Tz = 1.0/rates%inv_Tz
      zfactor = 1.0/(1.0-(tau_z/Tz))
      if( zfactor .lt. 0.0 ) then
        call out_io(s_abort$, r_name, &
             'FATAL ERROR: Negative emittance encountered: ', &
             'Try adjusting initial_blow_up or switch to derivatives method "der".')
        STOP
      endif
    endif

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

    if(counter .gt. SIZE(running_emit_x)) then
      if( abs(runavg_emit_x_A/runavg_emit_x_B - 1.) .lt. threshold ) then
        if( abs(runavg_emit_y_A/runavg_emit_y_B - 1.) .lt. threshold ) then
          if( abs(runavg_sigE_E_A/runavg_sigE_E_B - 1.) .lt. threshold ) then
            converged = .true.
          endif
        endif
      endif        
    endif        

    if(.not.converged) then
      ibsmode%a%emittance = emit_a
      ibsmode%b%emittance = emit_b
      ibsmode%sigE_E = sigE_E 
      ibsmode%z%emittance = ibsmode%sig_z * ibsmode%sigE_E 
    endif
  enddo
end subroutine ibs_equib_rlx

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  Subroutine ibs_equib_der(lat,ibs_sim_params,inmode,ibsmode,granularity)
!
!  Computes equilibrium beam sizes by calculating emittance growth rates from IBS growth rates.
!  Steps beam size through time till equilibrium is reached.
!
!  Modules needed:
!    use ibs_mod
!
!  Input:
!    lat             -- lat_struct: lattice for tracking
!      %param%n_part  -- Real: number of particles in bunch
!    ibs_sim_params   -- ibs_sim_param_struct: parameters for IBS calculation
!    inmode           -- normal_modes_struct: natural beam parameters 
!    granularity      -- Real: Step size for slicing lattice.  i.e. set to 1 to calculate IBS rates every 1 meter.
!                              Set to -1 to calculate element-by-element.
!
!  Output:
!    ibsmode          -- normal_modes_struct: beam parameters after IBS effects
!-

subroutine ibs_equib_der(lat,ibs_sim_params,inmode,ibsmode,granularity)
  ! Iterates to equilibrium beam conditions using derivatives

  implicit none

  type(lat_struct) :: lat
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(normal_modes_struct), intent(in) :: inmode
  type(normal_modes_struct), intent(out) :: ibsmode
  real(rp), intent(in) :: granularity
  type(ibs_struct) rates
  type(normal_modes_struct) :: naturalmode

  real(rp) time_for_one_turn
  real(rp) tau_a, tau_b, tau_z
  real(rp) Ta, Tb, Tz
  real(rp) emit_a, emit_b, sigE_E
  real(rp) emit_a0, emit_b0, sigE_E0
  real(rp) dadt, dbdt, dsEdt, dT
  real(rp) threshold
  real(rp) sigma_z0, L_ratio
  real(rp) sigma_z_vlassov
  real(rp) cur
  real(rp) alpha_c
  real(rp) energy
  real(rp) sigma_z_nat
  real(rp) pwd_ratio
  real(rp) sig_z

  real(rp) Vrf
  real(rp) U0
  real(rp) rf_freq, rf_omega

  logical converged
  integer i

  !natural mode here means emittances before IBS effects
  naturalmode = inmode

  !compute the SR betatron damping times
  time_for_one_turn = lat%param%total_length / c_light
  tau_a = time_for_one_turn / naturalmode%a%alpha_damp
  tau_b = time_for_one_turn / naturalmode%b%alpha_damp
  tau_z = time_for_one_turn / naturalmode%z%alpha_damp

  cur = lat%param%n_part*e_charge/time_for_one_turn
  alpha_c = inmode%synch_int(1)/lat%param%total_length
  energy = lat%ele(0)%value(E_TOT$)

  Vrf = 0.0d0
  do i=1,lat%n_ele_track
    if(lat%ele(i)%key .eq. rfcavity$) then
      Vrf = Vrf + lat%ele(i)%value(voltage$)
      rf_freq = lat%ele(i)%value(rf_frequency$)
    endif
  enddo
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
  dT = tau_a / 10.0 !Time to advance per iteration

  ibsmode = naturalmode
  do WHILE(.not.converged)
    if( ibs_sim_params%do_pwd ) then
      ! call bl_via_vlassov(cur,alpha_c,energy,ibsmode%sigE_E,Vrf,rf_omega,U0,lat%param%total_length, &
      !                     ibs_sim_params%resistance,ibs_sim_params%inductance,sigma_z_vlassov)
      call bl_via_mat(lat, ibs_sim_params, ibsmode, sig_z)
      ibsmode%sig_z = sig_z
    else
      ibsmode%sig_z = L_ratio * ibsmode%sigE_E
    endif

    do i=1, lat%n_ele_track
      lat%ele(i)%a%emit = ibsmode%a%emittance
      lat%ele(i)%b%emit = ibsmode%b%emittance
      lat%ele(i)%z%sigma = ibsmode%sig_z
      lat%ele(i)%z%sigma_p = ibsmode%sigE_E
      lat%ele(i)%z%emit = ibsmode%sig_z * ibsmode%sigE_E
    enddo

    call ibs_rates1turn(lat,ibs_sim_params,rates,granularity)

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

    if( (dadt*dT)/emit_a .lt. threshold ) then
      if( (dbdt*dT)/emit_b .lt. threshold ) then
        if( (dsEdt*dT)/ibsmode%sigE_E .lt. threshold ) then
          converged = .true.
        endif
      endif
    endif
    if(.not.converged) then
      ibsmode%a%emittance = emit_a + dadt*dT
      ibsmode%b%emittance = emit_b + dbdt*dT
      ibsmode%sigE_E = ibsmode%sigE_E + dsEdt*dT
      ibsmode%z%emittance = ibsmode%sig_z * ibsmode%sigE_E 
    endif
  enddo
end subroutine ibs_equib_der

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  Subroutine ibs_lifetime(lat,ibs_sim_params,maxratio,lifetime,granularity)
!
!  This module computes the beam lifetime due to
!  the diffusion process according to equation 12
!  from page 129 of The Handbook for Accelerator
!  Physics and Engineering 2nd edition.
!
!  Input:
!    lat                            -- lat_struct: lattice for tracking.
!    ibs_sim_params                 -- ibs_sim_param_struct: parameters for calculation of IBS rates.
!    maxratio(ibs_maxratio_struct)  -- Ax,y,p/sigma_x,y,p where Ax,y,p
!                                      is the maximum sigma.  Note that this quantity
!                                      is just the ratio, not the ratio squared.  For
!                                      example, maxratio%Rx = 1.1 says that the maximum
!                                      acceptable beamsize is 10% larger than the beamsize
!                                      before IBS effects.
!    granularity                    -- Step size when slicing lattice.  -1 for element-by-element.
!
!  Output:
!    lifetime(ibs_lifetime_struct) --structure returning IBS lifetimes
!-

subroutine ibs_lifetime(lat,ibs_sim_params,maxratio,lifetime,granularity)
  implicit none

  type(lat_struct) :: lat
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_maxratio_struct), intent(in) :: maxratio
  type(ibs_lifetime_struct), intent(out) :: lifetime
  real(rp), intent(in) :: granularity

  type(ibs_struct) rates

  real(rp) Rx, Ry, R_p

  Rx = maxratio%rx**2
  Ry = maxratio%ry**2
  R_p = maxratio%r_p**2

  call ibs_rates1turn(lat, ibs_sim_params, rates, granularity)

  lifetime%Tlx = exp(Rx)/2/Rx/rates%inv_Ta
  lifetime%Tly = exp(Ry)/2/Ry/rates%inv_Tb
  lifetime%Tlp = exp(R_p)/2/R_p/rates%inv_Tz
end subroutine ibs_lifetime

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  Subroutine ibs_delta_calc (lat, ix, ibs_sim_params, sigma_mat, delta_sigma_energy, delta_emit_a, delta_emit_b)
! 
!  Calculates change in energy spread and emittances due to IBS for a single element.
!
!  Input:
!    lat              -- lat_struct: lattice for tracking
!      %param%n_part  -- real(rp): number of particles in bunch
!    ix               -- integer: index of element to use: lat%ele(ix)
!    ibs_sim_params   -- ibs_sim_params_struct: parameters for calculation of IBS rates.
!    sigma_mat(6,6)   -- real(rp), optional: Beam's sigma matrix. Required for 'kubo' method.
!
!  Output:
!    delta_sigma_energy -- real(rp), optional: change in energy spread in eV
!    delta_emit_a       -- real(rp), optional: change in a-mode emittance (geometric)
!    delta_emit_b       -- real(rp), optional: change in b-mode emittance (geometric)
!-

subroutine ibs_delta_calc (lat, ix, ibs_sim_params, sigma_mat, delta_sigma_energy, delta_emit_a, delta_emit_b)
use mode3_mod, only: get_emit_from_sigma_mat

implicit none

type(lat_struct), intent(in) :: lat
type (ele_struct), pointer :: ele
type(ibs_sim_param_struct) :: ibs_sim_params
type(ibs_struct) :: rates1ele
real(rp), optional :: sigma_mat(6,6)
real(rp) :: dt, emit(3), emit_updated(3), sigma_mat_updated(6,6)
real(rp), optional :: delta_sigma_energy, delta_emit_a, delta_emit_b
integer :: ix
character(20) :: r_name = 'ibs_delta_calc'
logical err_flag
!
ele => lat%ele(ix)

if ( ibs_sim_params%formula == 'kubo') then
  if (.not. present(sigma_mat)) call out_io(s_fatal$, r_name, 'sigma_mat must be present for kubo calc')
  ! Extracted from kubo1_twiss_wrapper. Use sigma matrix only. 
  call kubo1(sigma_mat, ibs_sim_params, sigma_mat_updated, ele%value(E_TOT$), lat%param%n_part)
  ! Get emittances
  call get_emit_from_sigma_mat(sigma_mat, emit, err_flag=err_flag)
  call get_emit_from_sigma_mat(sigma_mat_updated, emit_updated, err_flag=err_flag)
  if (present(delta_emit_a))  delta_emit_a = emit_updated(1) - emit(1)
  if (present(delta_emit_b))  delta_emit_b = emit_updated(2) - emit(2)
  ! Slice energy spread is sqrt(s_66 - s_56^2/s_55)
  if (present(delta_sigma_energy)) delta_sigma_energy = (sqrt(sigma_mat_updated(6,6) - sigma_mat_updated(5,6)**2 / sigma_mat_updated(5,5)) &
        -  sqrt(sigma_mat(6,6) - sigma_mat(5,6)**2 / sigma_mat(5,5)))*ele%value(E_TOT$)
else
  call ibs1(lat, ibs_sim_params, rates1ele, i=ix)
  if (present(delta_sigma_energy)) delta_sigma_energy = ele%value(l$)/c_light*rates1ele%inv_Tz*ele%z%sigma_p*ele%value(E_TOT$)
  if (present(delta_emit_a))       delta_emit_a     = 2*ele%value(l$)/c_light*rates1ele%inv_Ta*ele%a%emit
  if (present(delta_emit_b))       delta_emit_b     = 2*ele%value(l$)/c_light*rates1ele%inv_Tb*ele%b%emit
endif

end subroutine ibs_delta_calc

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  Subroutine ibs_rates1turn(lat, ibs_sim_params, rates1turn, granularity)
!
!  Calculates IBS risetimes for given lat
!  This is basically a front-end for the various formulas 
!  available in this module of calculating IBS rates.
!
!  Input:
!    lat              -- lat_struct: lattice for tracking.
!      %param$n_part  -- Real: number of particles in bunch.
!    ibs_sim_params   -- ibs_sim_param_struct: parameters for IBS calculation.
!    granularity      -- real(rp): slice length.  -1 for element-by-element.
!
!  Output:
!    rates1turn       -- ibs_struct: ibs rates for onr turn on the lattice.
!-

subroutine ibs_rates1turn(lat, ibs_sim_params, rates1turn, granularity)
  implicit none

  real(rp) sum_inv_Tz, sum_inv_Ta, sum_inv_Tb
  type(lat_struct) :: lat
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct) :: rates1ele
  type(ibs_struct), intent(out) :: rates1turn
  real(rp), intent(in) :: granularity
  integer i
  integer n_steps
  real(rp), ALLOCATABLE :: steps(:)
  real(rp) step_size
  real(rp) length_multiplier


  sum_inv_Tz = 0.0
  sum_inv_Ta = 0.0
  sum_inv_Tb = 0.0


  if( granularity .lt. 0.0 ) then
    do i=1,lat%n_ele_track
      if(lat%ele(i)%value(l$) .eq. 0.0) then
        CYCLE
      endif
      
      call ibs1(lat, ibs_sim_params, rates1ele, i=i)

      !length_multiplier = lat%ele(i)%value(l$)/2.0 + lat%ele(i+1)%value(l$)/2.0
      length_multiplier = lat%ele(i)%value(l$)
      sum_inv_Tz = sum_inv_Tz + rates1ele%inv_Tz * length_multiplier
      sum_inv_Ta = sum_inv_Ta + rates1ele%inv_Ta * length_multiplier
      sum_inv_Tb = sum_inv_Tb + rates1ele%inv_Tb * length_multiplier
    enddo
  else
    n_steps = CEILING(lat%param%total_length / granularity)
    step_size = lat%param%total_length / n_steps

    ALLOCATE(steps(1:n_steps))
    do i=1,n_steps-1
      steps(i) = i*step_size
    enddo
    steps(n_steps) = lat%param%total_length


    do i=1,n_steps
      call ibs1(lat, ibs_sim_params, rates1ele, s=steps(i))

      sum_inv_Tz = sum_inv_Tz + rates1ele%inv_Tz * step_size
      sum_inv_Ta = sum_inv_Ta + rates1ele%inv_Ta * step_size
      sum_inv_Tb = sum_inv_Tb + rates1ele%inv_Tb * step_size
    enddo
    DEALLOCATE(steps)
  endif

  rates1turn%inv_Tz = sum_inv_Tz / lat%param%total_length
  rates1turn%inv_Ta = sum_inv_Ta / lat%param%total_length
  rates1turn%inv_Tb = sum_inv_Tb / lat%param%total_length
end subroutine ibs_rates1turn

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  Subroutine ibs_blowup1turn(lat, ibs_sim_params)
!
!  Updates beam emittances with effect of IBS for
!  one turn on the lattice.
!
!  Input:
!    lat                  -- lat_struct: lattice
!       %ele(:)%a%emit    -- real(rp): initial a-mode emittance
!       %ele(:)%b%emit    -- real(rp): initial b-mode emittance
!       %ele(:)%z%sigma_p -- real(rp): initial energy spread
!    ibs_sim_params       -- ibs_sim_param_struct: Parameters for calculation of IBS rates
!
!  Output:
!    lat                  -- lat_struct: lattice
!       %ele(:)%a%emit    -- real(rp): a-mode emittance after 1 turn
!       %ele(:)%b%emit    -- real(rp): b-mode emittance after 1 turn
!       %ele(:)%z%sigma_p -- real(rp): energy spread after 1 turn
!-

subroutine ibs_blowup1turn(lat, ibs_sim_params)
  implicit none

  type(lat_struct) :: lat
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct) :: rates1ele
  integer i
  real(rp) delta_t, pp, gg, gamma1, gamma2

  do i=1,lat%n_ele_track
    delta_t = lat%ele(i)%value(l$)/c_light
    if(delta_t .gt. 0.) then
      rates1ele%inv_Ta = 0.0_rp
      rates1ele%inv_Tb = 0.0_rp
      rates1ele%inv_Tz = 0.0_rp

      call ibs1(lat, ibs_sim_params, rates1ele, i=i)

      call convert_total_energy_to(lat%ele(i)%value(E_TOT$), -1, gamma1)
      call convert_total_energy_to(lat%ele(i+1)%value(E_TOT$), -1, gamma2)
      gg = gamma1/gamma2
      pp = lat%ele(i)%value(p0c$)/lat%ele(i+1)%value(p0c$)
      lat%ele(i+1)%a%emit = lat%ele(i)%a%emit * (1 + 2.0*delta_t*rates1ele%inv_Ta) * gg
      lat%ele(i+1)%b%emit = lat%ele(i)%b%emit * (1 + 2.0*delta_t*rates1ele%inv_Tb) * gg
      lat%ele(i+1)%z%sigma_p = lat%ele(i)%z%sigma_p * (1 + delta_t*rates1ele%inv_Tz) * pp
    else
      lat%ele(i+1)%a%emit = lat%ele(i)%a%emit 
      lat%ele(i+1)%b%emit = lat%ele(i)%b%emit 
      lat%ele(i+1)%z%sigma_p = lat%ele(i)%z%sigma_p 
    endif
  enddo
end subroutine ibs_blowup1turn

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! subroutine ibs1(lat, ibs_sim_params, rates, i, s)
!
! Calculates IBS growth rates at some location in a lattice.
! The IBS rates are betatron growth rates.  That is, they are the rate of
! change in sigma_x, sigma_y, and sigma_p.  The emittance growth
! rate is twice the betatron growth rate.
! 1/T_emit = 2/T_betatron.
! eg  emit(t) = emit_0 * exp(-2*t/T_betatron) because emit = sigma^2/beta
!
!  Available IBS formulas (ibs_sim_params%formula):
!    cimp - Completely Integrated Modified Piwinski
!    bjmt - Bjorken-Mtingwa formulation general to bunched beams (time consuming)
!    bane - Bane approximation of Bjorken-Mtingwa formulation
!    mpzt - Modified Piwinski with Zotter's Integral
!    mpxx - Modified Piwinski with a constant Coulomb log.
!    kubo - Kubo and Oide's sigma matrix-based
!
! Either i or s, but not both, must be specified.
!
! Input:
!   lat                       - lat_struct
!      %ele(i)%a%emit         - each must be populated with a-mode emittance.
!      %ele(i)%b%emit         - each must be populated with a-mode emittance.
!      %ele(i)%z%sigma_p      - each must be populated with a-mode emittance.
!   ibs_sim_params            - ibs_sim_param_struct: parameters for IBS calculation
!   i                         - integer: element index of location to calculate IBS rates.
!   s                         - real(rp): location in meters to calculate IBS rates.
!
! Output:
!   rates$inv_Ta              - real(rp): 1/Ta, where Ta is rise time of a betatron mode
!   rates$inv_Tb              - real(rp): 1/Ta, where Ta is rise time of b betatron mode
!   rates$inv_Tz              - real(rp): 1/Ta, where Ta is rise time of energy spread
!-

subroutine ibs1(lat, ibs_sim_params, rates, i, s)

  implicit none

  type(lat_struct) :: lat
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct) :: rates
  integer, intent(in), OPTIONAL :: i
  real(rp), intent(in), OPTIONAL :: s

  type(ele_struct), pointer :: ele
  type(ele_struct), target :: stubele
  real(rp) n_part

  n_part = lat%param%n_part

  if(ibs_sim_params%formula == 'kubo') then
    if( PRESENT(i) ) then
      call kubo1_twiss_wrapper(lat, ibs_sim_params, rates, ix=i)
    elseif( PRESENT(s) ) then
      call kubo1_twiss_wrapper(lat, ibs_sim_params, rates, s=s)
    endif
  else
    if(PRESENT(i) .and. .not.PRESENT(s)) then
      if(lat%ele(i)%value(l$) .eq. 0.0) then
        rates%inv_Tz = 0.0
        rates%inv_Ta = 0.0
        rates%inv_Tb = 0.0
        return
      else
        ele => lat%ele(i)
      endif
    elseif(PRESENT(s) .and. .not.PRESENT(i)) then
      call twiss_and_track_at_s(lat,s,stubele)
      ele => stubele
    else
      write(*,*) "FATAL ERROR in ibs_mod: Either i or s (and not both) must be specified"
      STOP
    endif

    if( ibs_sim_params%set_dispersion ) then
      ele%b%eta = ibs_sim_params%eta_set
      ele%b%etap = ibs_sim_params%etap_set
    endif

    if(ibs_sim_params%formula == 'cimp') then
      call cimp1(ele, ibs_sim_params, rates, n_part)
    elseif(ibs_sim_params%formula == 'bjmt') then
      call bjmt1(ele, ibs_sim_params, rates, n_part)
    elseif(ibs_sim_params%formula == 'bane') then
      call bane1(ele, ibs_sim_params, rates, n_part)
    elseif(ibs_sim_params%formula == 'mpzt') then
      call mpzt1(ele, ibs_sim_params, rates, n_part)
    elseif(ibs_sim_params%formula == 'mpxx') then
      call mpxx1(ele, ibs_sim_params, rates, n_part)
    else
      write(*,*) "Invalid IBS formula selected ... terminating"
      STOP
    endif
  endif

end subroutine ibs1

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine bjmt1(ele, ibs_sim_params, rates, n_part)
!
!  This is an implementation of equations 1-9 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is the most general form (for bunched beams) of the 
!  Bjorken-Mtingwa IBS formulation.
!
!  This formulation is one of the slowest methods for calculating IBS rates.
!
!  rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    ibs_sim_params    - ibs_sim_params_struct: parameters for IBS calculation
!    n_part            - real(rp): number of particles in the bunch.
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine bjmt1(ele, ibs_sim_params, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct), intent(out) :: rates
  type(ele_struct), target :: stubele
  real(rp) n_part

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_y
  real(rp) Dx, Dy, Dxp, Dyp
  real(rp) alpha_a, alpha_b, coulomb_log
  real(rp) big_A
  real(rp) Hx, Hy
  real(rp) inv_Tz, inv_Ta, inv_Tb

  real(rp) phi_h, phi_v

  real(c_double) :: Lp(3,3), Lh(3,3), Lv(3,3), L(3,3)
  real(c_double), target :: mats(2,3,3)

  real(fgsl_double), target :: Elpha
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, -1, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  Dx = ele%a%eta
  Dy = ele%b%eta
  Dxp = ele%a%etap
  Dyp = ele%b%etap
  sigma_y = sqrt(beta_b*emit_b + (Dy*sigma_p)**2)

  !coulomb_log = log( (gamma**2)*sigma_y*emit_a/r_e/beta_a )
  call multi_coulomb_log(ibs_sim_params, ele, coulomb_log, n_part)

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

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

end subroutine bjmt1

function bjmt_integrand(xx,params) bind(c)
  real(c_double), value :: xx
  real(c_double) :: x
  type(c_ptr), value :: params
  real(c_double) :: bjmt_integrand

  real(c_double), pointer :: mats(:,:,:)

  real(c_double) :: im(3,3)
  real(c_double) :: L(3,3)
  real(c_double) :: Li(3,3)
  real(c_double) :: Li_im(3,3)

  real(c_double) :: Tr_Li
  real(c_double) :: Tr_im
  real(c_double) :: Tr_Li_im
  real(c_double) :: Det_L_l

  call c_f_pointer(params,mats,[2,3,3])

  L = mats(1,:,:)
  Li = mats(2,:,:)

  x = exp(xx) !change of variables

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

  bjmt_integrand = sqrt(x/Det_L_l) * ( Tr_Li*Tr_im - 3.0_rp*Tr_Li_im )
  bjmt_integrand = bjmt_integrand * x  !COV
end function bjmt_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! Subroutine kubo1_twiss_wrapper(lat, ibs_sim_params, rates, ix, s)
!
! This is a wrapper for the sigma-matrix based kubo1 IBS subroutine.
! This subrouine generates a sigma-matrix from the one-turn transfer matrix
! using the normal mode emittances in lat%ele(ix).  It then calls
! kubo1, which updates the sigma-matrix with the effect of IBS over the 
! length of the element, and computes a time rate of change of the emittances.
!
! Input:
!   lat                   - lat_struct
!      %ele(ix)           - ele_struct: element at which calculation is performed
!      %ele(ix)%r(6,6,1)  - real(rp) (optional): if associated, then assumed to contain one-turn
!                                                transfer matrix.  Faster than recomputing one-turn
!                                                map on every turn.
!   ibs_sim_params        - ibs_sim_params_struct: Parameters for IBS calculation
!   ix                    - integer: element index at which to perform calculation
!   s                     - real(rp): location in meters where to perform calculation. NOT CURRENTLY SUPPORTED.
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine kubo1_twiss_wrapper(lat, ibs_sim_params, rates, ix, s)
  ! Some parts of this subroutine are patterned from the SAD accelerator code.

  use mode3_mod
  use longitudinal_profile_mod

  implicit none

  ! Ps.sigma_mat rearranges the sigma matrix so that all xx terms are in top left block,
  ! all px terms are in top right block, and pp terms are in bottom left block.  Lower right
  ! block is px'.
  integer, parameter :: Ps(6,6) = reshape( [1,0,0,0,0,0, 0,0,1,0,0,0, 0,0,0,0,1,0, &
                                            0,1,0,0,0,0, 0,0,0,1,0,0, 0,0,0,0,0,1],[6,6] )

  type(lat_struct) :: lat
  type(ibs_sim_param_struct) ibs_sim_params
  type(ibs_struct) :: rates
  integer, OPTIONAL :: ix
  real(rp), OPTIONAL :: s
  real(rp) dt

  type(ele_struct), target :: stubele
  integer ix_use
  real(rp) sigma_mat(6,6)
  real(rp) sigma_mat_updated(6,6)
  real(rp) Nmat(6,6)
  real(rp) t6(6,6)
  real(rp) W(6,6)
  real(rp) normal(3)
  real(rp) err_rp

  integer i

  integer error
  logical lerror

  type(normal_modes_struct) mode
  type(ele_struct), pointer :: ele
  real(rp) energy, gamma, KE, rbeta

  logical lerr, ok, err_flag

  if( ibs_sim_params%set_dispersion ) then
    W = 0.0_rp
    do i=1,6
      W(i,i) = 1.0_rp
    enddo
    W(3,6) = -ibs_sim_params%eta_set
    W(4,6) = -ibs_sim_params%etap_set
    W(5,3) =  ibs_sim_params%etap_set
    W(5,4) = -ibs_sim_params%eta_set
  endif

  if(PRESENT(ix) .and. .not.PRESENT(s)) then
    ix_use = ix
    ele => lat%ele(ix_use)
  elseif(PRESENT(s) .and. .not.PRESENT(ix)) then
    ix_use = element_at_s(lat,s,.false.)
    call twiss_and_track_at_s(lat,s,stubele)
    ele => stubele
  else
    write(*,*) "FATAL ERROR in ibs_mod: Either ix or s (and not both) must be specified"
    STOP
  endif

  ! make sigma matrix in x,px,y,py,z,pz form
  if( ibs_sim_params%use_t6_cache ) then
    if( ASSOCIATED(lat%ele(ix_use)%r) ) then
      t6 = lat%ele(ix_use)%r(:,:,1)
    else
      write(*,'(A)') "use_t6_cache is true, but lat%ele(ix_use)%r not allocated.  Defaulting to transfer_matrix_calc."
      call transfer_matrix_calc (lat, .true., t6, ix1=ix_use, one_turn=.TRUE.)
    endif
  else
    call transfer_matrix_calc (lat, .true., t6, ix1=ix_use, one_turn=.TRUE.)
  endif

  if( PRESENT(s) ) then
    t6 = matmul(ele%mat6,  matmul(dagger(lat%ele(ix_use)%mat6), matmul(t6, matmul(lat%ele(ix_use)%mat6, dagger(ele%mat6)))))
  endif

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, -1, gamma, KE, rbeta)

  ! t6 = pwd_mat(lat, t6, ibs_sim_params%inductance, ele%z%sigma)
  ! call transfer_matrix_calc_special(lat, .true., t6, ix1=ix_use, one_turn=.TRUE., inductance=ibs_sim_params%inductance, sig_z=lat%ele(ix_use)%z%sigma)

  if( ibs_sim_params%set_dispersion ) then
    t6 = matmul(t6,W) 
  endif
  mode%a%emittance = ele%a%emit
  mode%b%emittance = ele%b%emit
  mode%z%emittance = ele%z%emit
  mode%a%tune = mod(lat%ele(lat%n_ele_track)%a%phi,twopi)
  mode%b%tune = mod(lat%ele(lat%n_ele_track)%b%phi,twopi)
  call make_smat_from_abc(t6, mode, sigma_mat, lerr, Nmat)

  if( lerr ) then
    write(*,'(A)') "BAD: make_smat_from_abc failed"
    rates%inv_Ta = 0.0d0
    rates%inv_Tb = 0.0d0
    rates%inv_Tz = 0.0d0
    return
  endif

  call kubo1(sigma_mat, ibs_sim_params, sigma_mat_updated, energy, lat%param%n_part)

  call get_emit_from_sigma_mat(sigma_mat_updated,normal, Nmat, err_flag)

  ! D_mat = matmul(dagger(N_mat),matmul(sigma_mat_updated,dagger(transpose(N_mat))))
  ! normal(1) = D_mat(1,1)
  ! normal(2) = D_mat(3,3)
  ! normal(3) = D_mat(5,5)
  ! dt = 1.0 / c_light
  ! rates%inv_Ta = ((normal(1)) / dt)/ele%a%emit/2.0d0
  ! rates%inv_Tb = ((normal(2)) / dt)/ele%b%emit/2.0d0
  ! rates%inv_Tz = ((normal(3)) / dt)/ele%z%emit/2.0d0
  
  dt = 1.0 / c_light
  rates%inv_Ta = ((normal(1) - ele%a%emit) / dt)/ele%a%emit/2.0d0
  rates%inv_Tb = ((normal(2) - ele%b%emit) / dt)/ele%b%emit/2.0d0
  rates%inv_Tz = ((normal(3) - ele%z%emit) / dt)/ele%z%emit/2.0d0

  !call deallocate_ele_pointers(ele)
end subroutine kubo1_twiss_wrapper

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! This is a sigma matrix based IBS calculation.
! It takes the beam sigma matrix and updates it for IBS effects applied over the length element_length. 
!
! Input:
!   sigma_mat(6,6)         - real(rp): beam sigma_matrix
!   ibs_sim_params         - ibs_sim_params_struct: parameters for IBS calculation
!   element_length         - real(rp): length of element
!   energy                 - real(rp): beam energy in eV
!   n_part                 - real(rp): number of particles in the bunch
!
! Output:
!   sigma_mat_updated(6,6) - real(rp): sigma matrix updated with effect of IBS
!-

subroutine kubo1(sigma_mat, ibs_sim_params, sigma_mat_updated, energy, n_part)
  ! Some parts of this subroutine are patterned from the SAD accelerator code.

  use eigen_mod
  use LA_PRECISION, only: WP => DP
  use f95_lapack

  implicit none

  ! Ps.sigma_mat rearranges the sigma matrix so that all xx terms are in top left block,
  ! all px terms are in top right block, and pp terms are in bottom left block.  Lower right
  ! block is px'.
  integer, parameter :: Ps(6,6) = reshape( [1,0,0,0,0,0, 0,0,1,0,0,0, 0,0,0,0,1,0, &
                                            0,1,0,0,0,0, 0,0,0,1,0,0, 0,0,0,0,0,1],[6,6] )

  real(rp) sigma_mat(6,6)
  type(ibs_sim_param_struct) ibs_sim_params
  real(rp) sigma_mat_updated(6,6)
  real(rp) n_part
  real(rp) energy

  real(rp) Tboost(6,6)
  real(rp) Tboost_inv(6,6)
  real(rp) Tspat(6,6)
  real(rp) Tspat_inv(6,6)
  real(rp) spatial_sigma_update(6,6)
  real(rp) sigma_update(6,6)
  real(rp) err_mat(6,6)
  real(rp) spatial_sigma_mat(6,6)
  real(rp) boosted_sigma_mat(6,6)
  real(rp) ar_sigma_mat(6,6)
  real(rp) sig_xx(3,3)
  real(rp) sig_xx_inv(3,3)
  real(rp) sig_xp(3,3)
  real(rp) sig_pp(3,3)
  real(rp) sig_pp_update(3,3)
  real(rp) sig_pl(3,3)
  real(rp) u(3)
  real(rp) R(3,3)
  integer etypes(3)
  real(rp) g1, g2, g3
  real(rp) clog
  real(rp) cI
  real(rp) Dw(3,3)
  real(rp) vol, vol1, pvol
  real(rp) ptrans, bn, bm, bmin, bmax
  real(rp) bmin1, bmin2

  integer error
  integer i, j, k

  real(rp), parameter :: pi_2 = pi/2.0d0

  !GSL for integrator
  type(fgsl_integration_workspace) :: integ_wk
  type(c_ptr) :: ptr
  type(fgsl_function) :: integrand_ready
  integer(fgsl_int) :: fgsl_status
  real(fgsl_double), target :: args(3)
  real(fgsl_double) :: abserr
  real(fgsl_double) :: integration_result

  real(rp) gamma, KE, rbeta

  logical ok

  call convert_total_energy_to(energy, -1, gamma, KE, rbeta)

  ! ! make transfer matrix from canonical to spatial coordinates
  ! Tspat = 0.0d0
  ! do i=1,6
  !   Tspat(i,i) = 1.0d0
  ! enddo
  ! !  Tspat(1,2) = element_length/4
  ! !  Tspat(3,4) = element_length/4
  ! !  Tspat(5,6) = element_length/4/gamma/gamma
  ! spatial_sigma_mat = matmul(Tspat,matmul(sigma_mat,transpose(Tspat)))
  spatial_sigma_mat = sigma_mat

  ! boost sigma matrix to COM frame of bunch
  Tboost = 0.0d0
  do i=1,6
    Tboost(i,i) = 1.0d0
  enddo
  Tboost(5,5) = 1.0d0 + gamma*gamma/(1.0d0+gamma)
  Tboost(6,6) = 1.0d0 / gamma
  boosted_sigma_mat = matmul(Tboost,matmul(spatial_sigma_mat,transpose(Tboost)))

  ! permute sigma matrix to x,y,z,px,py,pz form
  ar_sigma_mat = matmul(matmul(transpose(Ps),boosted_sigma_mat),Ps)
  sig_xx = ar_sigma_mat(1:3,1:3)
  sig_xp = ar_sigma_mat(1:3,4:6)
  sig_pp = ar_sigma_mat(4:6,4:6)
  !call eigensys(sig_xx, u, R, etypes, 3, error)
  R=sig_xx  ! LA_SYEV destroys the contents of R
  call LA_SYEV(R,u,JOBZ='V')  !evals and evecs of symmetric real matrix

  vol1 = sqrt(u(1)*u(2)*u(3))
  vol = sqrt(4.0d0*pi)**3 * vol1
  bm = sqrt(min( u(1), u(2), u(3) ))  !minimum beam dimension
  call mat_inverse(sig_xx,sig_xx_inv,ok)
  if( .not. ok ) then
    write(*,*) "BAD: Could not invert sig_xx"
    sigma_mat_updated = sigma_mat
    return
  endif

  !Get local momentum matrix
  sig_pl = sig_pp - matmul(transpose(sig_xp),matmul(sig_xx_inv,sig_xp))

  !Get eigen vectors of local momentum matrix
  call eigensys(sig_pl, u, R, etypes, 3, error)

  !R=sig_pl
  !call LA_SYEV(R,u,JOBZ='V',INFO=error)  !evals and evecs of symmetric real matrix
  !LA_SYEV seems to be less robust than eigensys.
  u(1) = max(u(1),1.0d-20)
  u(2) = max(u(2),1.0d-20)
  u(3) = max(u(3),1.0d-20)

  if( error .ne. 0 ) then
    write(*,'(A,I6," ",A)') "BAD: Eigenvectors of local momentum matrix not found."
    sigma_mat_updated = sigma_mat
    return
  endif
  R=transpose(R) !R is defined as inverse of eigenvector matrix, and tr(R) = inv(r)
  ptrans = sqrt(u(1)+u(2)+u(3))
  pvol = sqrt(u(1)*u(2)*u(3))

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

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)
  !- End integration with fgsl

  bn = (vol/n_part)**(1.0d0/3.0d0)
  bmax = min( bm, bn )  !minimum dimension or debye radius
  if( ibs_sim_params%clog_to_use .NE. 1 ) then
    !kubo's tail cut formula
    bmin1 = r_e/(ptrans*gamma)**2
    bmin2 = sqrt(abs(vol/n_part/pi/(ptrans*c_light)/ibs_sim_params%tau_a))
    bmin = max( bmin1, bmin2 )
  else
    !no tail cut
    bmin = r_e/(ptrans*gamma)**2
  endif
  clog = log(bmax/bmin)

  cI = (r_e**2)*n_part*clog/4.0d0/pi/(gamma**4)/vol1/pvol 

  Dw = 0.0d0
  Dw(1,1) = cI*(g2-g1+g3-g1)
  Dw(2,2) = cI*(g1-g2+g3-g2)
  Dw(3,3) = cI*(g1-g3+g2-g3)

  !- Build update matrix
  sig_pp_update = matmul(matmul(R,Dw),transpose(R))
  sigma_update = 0.0d0
  do i=1,3
    do j=1,3
      sigma_update(i*2,j*2) = sig_pp_update(i,j)
    enddo
  enddo

  ! boost updates to lab frame
  call mat_inverse(Tboost,Tboost_inv,ok)
  if( .not. ok ) then
    write(*,*) "BAD: Could not invert Tboost"
    sigma_mat_updated = sigma_mat
    return
  endif
  spatial_sigma_update = matmul(Tboost_inv,matmul(sigma_update,transpose(Tboost_inv)))

  ! ! back to canonical coordinates
  ! call mat_inverse(Tspat, Tspat_inv, ok)
  ! sigma_update = matmul(Tspat_inv,matmul(spatial_sigma_update,transpose(Tspat_inv)))
  sigma_update = spatial_sigma_update

  !- Update sigma matrix
   sigma_mat_updated = sigma_mat + sigma_update

end subroutine kubo1

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------

function kubo_integrand(s, params) bind(c)
    real(c_double), value :: s
    type(c_ptr), value :: params
    real(c_double) :: kubo_integrand

    real(c_double), pointer :: args(:)
    real(c_double) u1,u2,u3

    call c_f_pointer(params,args,[3])
    u1 = args(1)
    u2 = args(2)
    u3 = args(3)

    kubo_integrand = (2.0d0*u1*sin(s)**2*cos(s)) / &
                     sqrt((sin(s)**2 + u1/u2*cos(s)**2)*(sin(s)**2+u1/u3*cos(s)**2))
end function kubo_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine bane1(ele, ibs_sim_params, rates, n_part)
!
!  This is an implementation of equations 10-15 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a high energy approximation of the Bjorken-Mtingwa IBS
!  formulation.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    ibs_sim_params    - ibs_sim_params_struct: parameters for IBS calculation
!    n_part            - real(rp): number of particles in the bunch.
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine bane1(ele, ibs_sim_params, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct), intent(out) :: rates

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_b, sigma_b_beta
  real(rp) Da, Db, Dap, Dbp
  real(rp) alpha_a, alpha_b, coulomb_log
  real(rp) a, b, g_bane
  real(rp) n_part, big_A
  real(rp) sigma_H, Ha, Hb
  real(rp) inv_Tz, inv_Ta, inv_Tb

  real(fgsl_double), target :: Elpha
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, -1, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*n_part/16.0/(gamma**3)/(emit_a**(3./4.))/(emit_b**(3./4.))/sigma_z/(sigma_p**3)

  beta_a = ele%a%beta
  beta_b = ele%b%beta
  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  Dap = ele%a%etap
  Dbp = ele%b%etap
  Da = ele%a%eta
  Db = ele%b%eta
  sigma_b_beta = sqrt(beta_b * emit_b)
  sigma_b = sqrt(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  !coulomb_log = log( (gamma**2)*sigma_b*emit_a/r_e/beta_a )

  call multi_coulomb_log(ibs_sim_params, ele, coulomb_log, n_part)

  Ha = ( (Da**2) + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( (Db**2) + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b
  sigma_H = 1.0/sqrt(1.0/(sigma_p**2) + Ha/emit_a + Hb/emit_b)

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)

  Elpha = a/b
  ptr = c_loc(Elpha)
  integ_wk = fgsl_integration_workspace_alloc(limit)
  integrand_ready = fgsl_function_init(integrand, ptr)
  fgsl_status = fgsl_integration_qagiu(integrand_ready, 0.0d0, eps7, eps7, &
                                       limit, integ_wk, integration_result, abserr)
  g_bane = 2.0d0 * sqrt(Elpha) / pi * integration_result

  call fgsl_function_free(integrand_ready)
  call fgsl_integration_workspace_free(integ_wk)

  inv_Tz = big_A*coulomb_log*sigma_H*g_bane*((beta_a*beta_b)**(-1./4.))
  inv_Ta = (sigma_p**2)*Ha/emit_a*inv_Tz
  inv_Tb = (sigma_p**2)*Hb/emit_b*inv_Tz

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb

end subroutine bane1

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------

function integrand(x, params) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: integrand

    real(c_double), pointer :: alpha
    call c_f_pointer(params, alpha)

    integrand = 1.0_c_double / ( sqrt(1+x*x)*sqrt(alpha*alpha+x*x) )
end function integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine mpxx1(ele, ibs_sim_params, rates, n_part)
!
!  Modified Piwinski, further modified to treat Coulomb Log
!  in the same manner as Bjorken-Mtingwa, CIMP, Bane, Kubo & Oide, etc.
!  This formula is derived in Section 2.8.4 of Michael Ehrlichman's Graduate Thesis.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    ibs_sim_params    - ibs_sim_params_struct: parameters for IBS calculation
!    n_part            - real(rp): number of particles in the bunch.
!
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine mpxx1(ele, ibs_sim_params, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct), intent(out) :: rates

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_a, sigma_b, sigma_a_beta, sigma_b_beta
  real(rp) Da, Db, Dap, Dbp
  real(rp) alpha_a, alpha_b, coulomb_log
  real(rp) a,b,q
  real(rp) n_part, big_A
  real(rp) sigma_H, Ha, Hb
  real(rp) inv_Tz, inv_Ta, inv_Tb
  real(rp) fab, f1b, f1a

  real(fgsl_double), target :: args(1:2)
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, -1, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_a_beta = sqrt(beta_a * emit_a)
  sigma_b_beta = sqrt(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = sqrt(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = sqrt(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( Db**2 + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b

  sigma_H = 1.0/sqrt( 1.0/(sigma_p**2)+ Ha/emit_a + Hb/emit_b )

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)

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

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)
  !------------------------End calls to GSL integrator

  call multi_coulomb_log(ibs_sim_params, ele, coulomb_log, n_part)

  inv_Tz = coulomb_log * big_A * sigma_H**2 / sigma_p**2 * fab
  inv_Ta = coulomb_log * big_A * (f1b + Ha*sigma_H**2/emit_a*fab)
  inv_Tb = coulomb_log * big_A * (f1a + Hb*sigma_H**2/emit_b*fab)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb
end subroutine mpxx1

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------

function mpxx_integrand(x, params) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: mpxx_integrand

    real(c_double), pointer :: args(:)
    real(c_double) av, bv

    real(c_double) u

    call c_f_pointer(params,args,[2])
    av = args(1)
    bv = args(2)

    mpxx_integrand = 8.0d0*pi*(1-3.0d0*x*x) / &
                     sqrt(av*av+(1-av*av)*x*x) / &
                     sqrt(bv*bv+(1-bv*bv)*x*x)

end function mpxx_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine mpzt1(ele, ibs_sim_params, rates, n_part)
!
!  Modified Piwinski with Zotter's integral.  This is Piwinski's original derivation,
!  generalized to take the derivatives of the optics functions.  Also, Piwinski's
!  original cumbersome triple integral is reaplaced by Zotter's single integral.  Zotter's
!  integral is exact, and not an approximation.
!
!  rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    ibs_sim_params    - ibs_sim_params_struct: parameters for IBS calculation
!    n_part            - real(rp): number of particles in the bunch.
!
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine mpzt1(ele, ibs_sim_params, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct), intent(out) :: rates

  real(rp) sigma_p, emit_a, emit_b, sigma_z, energy
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_a, sigma_b, sigma_a_beta, sigma_b_beta
  real(rp) Da, Db, Dap, Dbp
  real(rp) alpha_a, alpha_b, coulomb_log
  real(rp) a,b,q
  real(rp) n_part, big_A
  real(rp) sigma_H, Ha, Hb
  real(rp) inv_Tz, inv_Ta, inv_Tb
  real(rp) fabq, f1bq, f1aq

  real(fgsl_double), target :: args(1:3)
  type(fgsl_function) :: integrand_ready
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr
  type(c_ptr) :: ptr
  type(fgsl_integration_workspace) :: integ_wk
  integer(fgsl_int) :: fgsl_status

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, -1, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_a_beta = sqrt(beta_a * emit_a)
  sigma_b_beta = sqrt(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = sqrt(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = sqrt(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  Hb = ( Db**2 + (beta_b*Dbp + alpha_b*Db)**2 ) / beta_b

  sigma_H = 1.0/sqrt( 1.0/(sigma_p**2)+ Ha/emit_a + Hb/emit_b )

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)
  q = sigma_H*rbeta*sqrt(2.0_rp*sigma_b/r_e)
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

  call fgsl_integration_workspace_free(integ_wk)
  call fgsl_function_free(integrand_ready)
  !------------------------End calls to GSL integrator

  inv_Tz = big_A * sigma_H**2 / sigma_p**2 * fabq
  inv_Ta = big_A * (f1bq + Ha*sigma_H**2/emit_a*fabq)
  inv_Tb = big_A * (f1aq + Hb*sigma_H**2/emit_b*fabq)

  rates%inv_Tz = inv_Tz
  rates%inv_Ta = inv_Ta
  rates%inv_Tb = inv_Tb
end subroutine mpzt1

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------

function zot_integrand(x, params) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: zot_integrand

    real(c_double), pointer :: args(:)
    real(c_double) av, bv, qv, P, Q

    real(c_double) u

    call c_f_pointer(params,args,[3])
    av = args(1)
    bv = args(2)
    qv = args(3)

    !COV to remove singularity at endpoints
    u=30.*((x**5.)/5. - (x**4.)/2. + (x**3.)/3.)

    P = sqrt(av*av + (1.-av*av)*u*u)
    Q = sqrt(bv*bv + (1.-bv*bv)*u*u)
    zot_integrand = 8.0_rp*pi * (1.0_rp-3.0_rp*u*u)/P/Q * &
                (2.0_rp*log(qv/2.0_rp*(1/P+1/Q)) - 0.577215665)

    zot_integrand = zot_integrand * 30.*(x**4. - 2.*(x**3.) + x**2.)  !COV to remove singularity

    !without COV
    !P = sqrt(av*av + (1-av*av)*x*x)
    !Q = sqrt(bv*bv + (1-bv*bv)*x*x)
    !zot_integrand = 8.0_rp*pi * (1.0_rp-3.0_rp*x*x)/P/Q * &
    !            (2.0_rp*log(qv/2.0_rp*(1/P+1/Q)) - 0.577215665)

end function zot_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
!  subroutine cimp1(ele, ibs_sim_params, rates, n_part)
!
!  This is an implementation of equations 34,38-40 from "Intrabeam
!  scattering formulas for high energy beams" Kubo,Mtingwa,Wolski.
!  It is a modified version of the Piwinski IBS formulation.
!  The integral (34) is handled with a piecewise interpolation generated
!  in mathematica.  The interpolation is accurate beyond 1% through it's
!  effective range (.0001 - 3000).
!
!  This is the quickest of the three IBS formuations in this module. 
!
!  rates returns betatron growth rates.  Multiply by two to get transverse emittance growth rates.
!
!  Input:
!    ele               - ele_struct: contains Twiss parameters used in IBS formula
!    ibs_sim_params    - ibs_sim_params_struct: parameters for IBS calculation
!    n_part            - real(rp): number of particles in the bunch.
!
!  Output:
!    rates             - ibs_struct
!         %inv_Ta      - real(rp): a-mode betatron growth rate.
!         %inv_Tb      - real(rp): b-mode betatron growth rate.
!         %inv_Tz      - real(rp): energy spread growth rate.
!-

subroutine cimp1(ele, ibs_sim_params, rates, n_part)

  implicit none

  type(ele_struct) :: ele
  type(ibs_sim_param_struct) :: ibs_sim_params
  type(ibs_struct), intent(out) :: rates
  real(rp) element_length, E_TOT, n_part

  real(rp) sigma_p, emit_a, emit_b, sigma_z
  real(rp) gamma, KE, rbeta, beta_a, beta_b
  real(rp) sigma_x, sigma_y, sigma_x_beta, sigma_y_beta
  real(rp) Dx, Dy, Dxp, Dyp
  real(rp) alpha_a, alpha_b, coulomb_log
  real(rp) a, b
  real(rp) big_A
  real(rp) sigma_H, Hx, Hy
  real(rp) inv_Tz, inv_Ta, inv_Tb
  real(rp) g_ab,g_ba
  real(rp) bminstar, bmax
  real(rp) energy

  !- Code specific to alternitive log representation
  real(rp) q, lnqa, lnqb
  !-

  energy = ele%value(E_TOT$)

  call convert_total_energy_to(energy, -1, gamma, KE, rbeta)

  sigma_p = ele%z%sigma_p
  sigma_z = ele%z%sigma
  emit_a = ele%a%emit
  emit_b = ele%b%emit

  big_A=(r_e**2)*c_light*n_part/64.0/(pi**2)/(rbeta**3)/(gamma**4)/emit_a/emit_b/sigma_z/sigma_p

  alpha_a = ele%a%alpha
  alpha_b = ele%b%alpha
  beta_a = ele%a%beta
  beta_b = ele%b%beta
  sigma_x_beta = sqrt(beta_a * emit_a)
  sigma_y_beta = sqrt(beta_b * emit_b)
  Dx = ele%a%eta
  Dy = ele%b%eta
  Dxp = ele%a%etap
  Dyp = ele%b%etap
  sigma_x = sqrt(sigma_x_beta**2 + (Dx**2)*(sigma_p**2))
  sigma_y = sqrt(sigma_y_beta**2 + (Dy**2)*(sigma_p**2))

  Hx = ( Dx**2 + (beta_a*Dxp + alpha_a*Dx)**2 ) / beta_a
  Hy = ( Dy**2 + (beta_b*Dyp + alpha_b*Dy)**2 ) / beta_b

  sigma_H = 1.0/sqrt( 1.0/(sigma_p**2)+ Hx/emit_a + Hy/emit_b )

  call multi_coulomb_log(ibs_sim_params, ele, coulomb_log, n_part)

  a = sigma_H/gamma*sqrt(beta_a/emit_a)
  b = sigma_H/gamma*sqrt(beta_b/emit_b)

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
!  rates%inv_Tz = 0.0
!  rates%inv_Ta = 0.0
!  rates%inv_Tb = 0.0

end subroutine cimp1

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
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

function g(u)
                                                                                   
  real(rp), intent(in) :: u
  real(rp) :: g
  real(rp), dimension(0:13) ::  o,p,pa,pb,qb,qa,q,ra,rb,r

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

  if    (u .gt. 3000.0) then
    write(*,*) "CRITICAL WARNING: interpolation range exceeded"
    g = 0.
  elseif(u .gt. 300.0) then
    g = o(0)+o(1)*u+o(2)*(u**2)+o(3)*(u**3)+o(4)*(u**4)+o(5)*(u**5)+ &
        o(6)*(u**6)+o(7)*(u**7)+o(8)*(u**8)+o(9)*(u**9)+o(10)*(u**10)+ &
        o(11)*(u**11)+o(12)*(1.0E-10)*(u**12)+o(13)*(1.0E-10)*(u**13)
  elseif(u .gt. 30.0)  then
    g = p(0)+p(1)*u+p(2)*(u**2)+p(3)*(u**3)+p(4)*(u**4)+p(5)*(u**5)+ &
        p(6)*(u**6)+p(7)*(u**7)+p(8)*(u**8)+p(9)*(u**9)+p(10)*(u**10)+ &
        p(11)*(u**11)+p(12)*(u**12)+p(13)*(u**13)
  elseif(u .gt. 10.0)  then
    g = pa(0)+pa(1)*u+pa(2)*(u**2)+pa(3)*(u**3)+pa(4)*(u**4)+pa(5)*(u**5)+ &
        pa(6)*(u**6)+pa(7)*(u**7)+pa(8)*(u**8)+pa(9)*(u**9)+pa(10)*(u**10)+ &
        pa(11)*(u**11)+pa(12)*(u**12)+pa(13)*(u**13)
  elseif(u .gt. 1.6)   then
    g = pb(0)+pb(1)*u+pb(2)*(u**2)+pb(3)*(u**3)+pb(4)*(u**4)+pb(5)*(u**5)+ &
        pb(6)*(u**6)+pb(7)*(u**7)+pb(8)*(u**8)+pb(9)*(u**9)+pb(10)*(u**10)+ &
        pb(11)*(u**11)+pb(12)*(u**12)+pb(13)*(u**13)
  elseif(u .gt. .20)   then
    g = qb(0)+qb(1)*u+qb(2)*(u**2)+qb(3)*(u**3)+qb(4)*(u**4)+qb(5)*(u**5)+ &
        qb(6)*(u**6)+qb(7)*(u**7)+qb(8)*(u**8)+qb(9)*(u**9)+qb(10)*(u**10)+ &
        qb(11)*(u**11)+qb(12)*(u**12)+qb(13)*(u**13)
  elseif(u .gt. .10)   then
    g = qa(0)+qa(1)*u+qa(2)*(u**2)+qa(3)*(u**3)+qa(4)*(u**4)+qa(5)*(u**5)+ &
        qa(6)*(u**6)+qa(7)*(u**7)+qa(8)*(u**8)+qa(9)*(u**9)+qa(10)*(u**10)+ &
        qa(11)*(u**11)+qa(12)*(u**12)+qa(13)*(u**13)
  elseif(u .gt. .02)  then
    g = q(0)+q(1)*u+q(2)*(u**2)+q(3)*(u**3)+q(4)*(u**4)+q(5)*(u**5)+ &
        q(6)*(u**6)+q(7)*(u**7)+q(8)*(u**8)+q(9)*(u**9)+q(10)*(u**10)+ &
        q(11)*(u**11)+q(12)*(u**12)+q(13)*(u**13)
  elseif(u .gt. .01)  then
    g = ra(0)+ra(1)*u+ra(2)*(u**2)+ra(3)*(u**3)+ra(4)*(u**4)+ra(5)*(u**5)+ &
        ra(6)*(u**6)+ra(7)*(u**7)+ra(8)*(u**8)+ra(9)*(u**9)+ra(10)*(u**10)+ &
        ra(11)*(u**11)+ra(12)*(u**12)+ra(13)*(u**13)
  elseif(u .gt. .001)  then
    g = rb(0)+rb(1)*u+rb(2)*(u**2)+rb(3)*(u**3)+rb(4)*(u**4)+rb(5)*(u**5)+ &
        rb(6)*(u**6)+rb(7)*(u**7)+rb(8)*(u**8)+rb(9)*(u**9)+rb(10)*(u**10)+ &
        rb(11)*(u**11)+rb(12)*(u**12)+rb(13)*(u**13)
  elseif(u .gt. .0001) then
    g = r(0)+r(1)*u+r(2)*(u**2)+r(3)*(u**3)+r(4)*(u**4)+r(5)*(u**5)+ &
        r(6)*(u**6)+r(7)*(u**7)+r(8)*(u**8)+r(9)*(u**9)+r(10)*(u**10)+ &
        r(11)*(u**11)+r(12)*(1.0E10)*(u**12)+r(13)*(1.0E10)*(u**13)
  else
    write(*,*) "CRITICAL WARNING: interpolation range exceeded"
    g = 0.
  endif
end function g

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! Subroutine multi_coulomb_log(ibs_sim_params, ele, coulomb_log, n_part)
!
! Calculates the value of the Coulomb log using various methods.
!
! ibs_sim_params%clog_to_use == 1   Classic coulomb log (pi/2 max scattering angle)
! ibs_sim_params%clog_to_use == 2   Integral based tail-cut prescribed by Raubenheimer.
! ibs_sim_params%clog_to_use == 3   Bane tail cut. 1 event/part/damping period.  
! ibs_sim_params%clog_to_use == 4   Kubo and Oide tail cut. Used by CesrTA publications.
!
! Input:
!   ibs_sim_params        - ibs_sim_params_struct: parameters for IBS calculation
!   ele                   - ele_struct: populated with Twiss parameters
!   n_part                - real(rp): number of particles in the bunch
! Output:
!   coulomb_log           - real(rp): Value of the Coulomb Logarithm
!-

subroutine multi_coulomb_log(ibs_sim_params, ele, coulomb_log, n_part)
  type(ibs_sim_param_struct) :: ibs_sim_params
  real(rp) coulomb_log

  type(ele_struct) ele
  real(rp) gamma, energy, g2
  real(rp) sigma_a, sigma_b, sigma_z, sigma_p, sp2
  real(rp) sigma_a_beta, sigma_b_beta
  real(rp) emit_a, emit_b
  real(rp) beta_a, beta_b
  real(rp) alpha_a, alpha_b
  real(rp) Da, Db, Dap, Dbp
  real(rp) n_part
  real(rp) bminstar
  real(rp) bmax
  real(rp) gamma_a, gamma_b
  real(rp) Ha, Hb
  real(rp) Bbar
  real(rp) u, v, w  !used for Raubenheimer's calculation
  real(rp) qmin, qmax
  real(rp) vol, debye

  !fgsl variables
  type(c_ptr) :: ptr
  real(fgsl_double), target :: args(1:3)
  type(fgsl_integration_workspace) :: integ_wk
  type(fgsl_function) :: integrand_ready
  integer(fgsl_int) :: fgsl_status
  real(fgsl_double) :: integration_result
  real(fgsl_double) :: abserr

  energy = ele%value(E_TOT$)
  call convert_total_energy_to(energy, -1, gamma=gamma)

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
  sigma_a_beta = sqrt(beta_a * emit_a)
  sigma_b_beta = sqrt(beta_b * emit_b)
  Da = ele%a%eta
  Db = ele%b%eta
  Dap = ele%a%etap
  Dbp = ele%b%etap
  sigma_a = sqrt(sigma_a_beta**2 + (Da**2)*(sigma_p**2))
  sigma_b = sqrt(sigma_b_beta**2 + (Db**2)*(sigma_p**2))

  if( ibs_sim_params%clog_to_use == 1 ) then
    !Classic Coulomb Log.
    coulomb_log = log( (gamma**2)*sigma_b*emit_a/r_e/beta_a )
  elseif( ibs_sim_params%clog_to_use == 2 ) then
    !Tail cut Raubenheimer gstar (ignores vertical dispersion !)
    Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
    g2 = gamma*gamma
    sp2 = sigma_p*sigma_p
    u = g2*( Ha/emit_a + 1.0/sp2 + beta_a/g2/emit_a + beta_b/g2/emit_b )
    v = g2*( Ha*beta_b/emit_a/emit_b + beta_a*beta_b/g2/emit_a/emit_b + beta_a/sp2/emit_a + beta_b/sp2/emit_b + Da*Da/emit_a/emit_a ) 
    w = g2*( Da*Da*beta_b/emit_a/emit_a/emit_b + beta_a*beta_b/sp2/emit_a/emit_b )

    Bbar = gamma*sqrt( gamma_a*emit_a + gamma_b*emit_b + sp2/g2)
    qmin = 2.0 * r_e / sigma_b / Bbar

    !FGSL integration
    ptr = c_loc(args)
    args = (/u,v,w/)
    integ_wk = fgsl_integration_workspace_alloc(limit)
    integrand_ready = fgsl_function_init(rclog_integrand, ptr)
    fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 40.0d0, eps7, eps7, &
                                           limit, 3, integ_wk, integration_result, abserr)

    call fgsl_function_free(integrand_ready)
    call fgsl_integration_workspace_free(integ_wk)

    qmax = sqrt( ibs_sim_params%tau_a*n_part*c_light*r_e*r_e/4.0/pi/g2/emit_a/emit_b/sigma_z/sigma_p * integration_result )

    coulomb_log = log(qmax/qmin)
  elseif( ibs_sim_params%clog_to_use == 3 ) then
    !Tail cut Bane form: SLAC-PUB-9227
    vol = (4.0d0*pi)**(3.0d0/2.0d0)*sigma_a*sigma_b*sigma_z*gamma
    bminstar = sqrt(vol/n_part/pi/(ibs_sim_params%tau_a/gamma)) / sqrt(c_light*gamma*sqrt(emit_a/beta_a))
    debye = (vol/n_part)**(1.0d0/3.0d0)
    bmax = min(min(sigma_a,sigma_b),debye)
    coulomb_log = log(bmax/bminstar)
  elseif( ibs_sim_params%clog_to_use == 4) then
    !Tail cut Kubo, but using relative velocity in all 3 dims, rather than just horizontal
    vol = (4.0d0*pi)**(3.0d0/2.0d0)*sigma_a*sigma_b*sigma_z*gamma
    Bbar = sqrt( gamma_a*emit_a + gamma_b*emit_b + sigma_p*sigma_p/gamma/gamma)  !beta in the rest frame
    bminstar = sqrt(vol/n_part/pi/(ibs_sim_params%tau_a/gamma)) / sqrt(c_light*gamma*Bbar)
    debye = (vol/n_part)**(1.0d0/3.0d0)
    bmax = min(min(sigma_a,sigma_b),debye)
    coulomb_log = log(bmax/bminstar)
  endif

end subroutine multi_coulomb_log

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------

function rclog_integrand(x, params) bind(c)
    real(c_double), value :: x
    type(c_ptr), value :: params
    real(c_double) :: rclog_integrand

    real(c_double), pointer :: args(:)
    real(c_double) u, v, w

    call c_f_pointer(params, args, [3])
    u = args(1)
    v = args(2)
    w = args(3)

    rclog_integrand = 2.0_rp/sqrt(exp(6.0_rp*x) + u*exp(4.0_rp*x) + v*exp(2.0_rp*x) + w)*exp(x)
end function rclog_integrand

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! Subroutine bl_via_vlassov(current,alpha,Energy,sigma_p,Vrf,omega,U0,circ,R,L,sigma_z)
!
! This is a frontend for get_bl_from_fwhm from longitudinal_profile_mod.
! See longitudinal_profile_mod for details.  In short, this implements a model of potential well distortion
! based on the Vlassov equation which uses an effective Resistive, Inductive, and Capacitive impedance.
!
! Input:
!   current       -- real(rp): Beam current in amps
!   alpha         -- real(rp): Momentum compaction
!   Energy        -- real(rp): beam energy
!   sigma_p       -- real(rp): energy spread
!   Vrf           -- real(rp): total RF voltage in Volts
!   omega         -- real(rp): rf frequency in radians/s
!   U0            -- real(rp): energy loss per turn (eV)
!   circ          -- real(rp): circumpherence
!   R             -- real(rp): Resistive part of effective impedance
!   L             -- real(rp): Inductive part of effective impedance
! Output:
!   sigma_z       -- real(rp): Bunch length. FWHM/TwoRootTwoLogTwo from bunch profile
!-
subroutine bl_via_vlassov(current,alpha,Energy,sigma_p,Vrf,omega,U0,circ,R,L,sigma_z)
  use longitudinal_profile_mod

  implicit none

  real(rp) current
  real(rp) alpha
  real(rp) Energy
  real(rp) sigma_p
  real(rp) Vrf
  real(rp) omega
  real(rp) U0
  real(rp) circ
  real(rp) R
  real(rp) L
  real(rp) sigma_z

  real(rp) delta_e
  real(rp) A
  real(rp) Q
  real(rp) T0
  real(rp) phi
  real(rp), parameter :: bound = 500.0d-12

  real(rp) args(1:8)

  delta_e = sigma_p * Energy
  T0 = circ/c_light
  Q = current * T0
  phi = -acos(U0/Vrf)
  A = Energy/delta_e/delta_e/alpha/T0

  args(1) = A
  args(2) = Vrf
  args(3) = Q
  args(4) = omega
  args(5) = phi
  args(6) = R
  args(7) = L
  args(8) = U0

  call get_bl_from_fwhm(bound,args,sigma_z)
end subroutine bl_via_vlassov

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! Subroutine bl_via_mat(lat, ibs_sim_params, mode, sig_z)
!
! Calculates bunch length while taking PWD effects into account.  
! PWD is approximated as a defocusing rf voltage.
!
! Input:
!   lat             - lat_struct
!   ibs_sim_params  - ibs_sim_params_struct: parameters for IBS calculation
!   mode            - normal_modes_modes_struct: energy spread, bunch length, and z emittance.
!
! Output:
!   sig_z           - real(rp): bunch length after taking PWD into account.
!-

subroutine bl_via_mat(lat, ibs_sim_params, mode, sig_z)

  use nr
  use mode3_mod
  use longitudinal_profile_mod

  implicit none

  type(lat_struct) lat
  type(ibs_sim_param_struct) ibs_sim_params
  type(normal_modes_struct) mode
  real(rp) sig_z
  real(rp) lb, mb, ub, min_val

  lb = 0.90 * mode%sig_z
  mb = mode%sig_z
  ub = 1.50 * mode%sig_z

  min_val = brent(lb,mb,ub,residual_pwd_sig_z, 0.00001d0, sig_z)

  contains

    function residual_pwd_sig_z(zz)
      real(rp) residual_pwd_sig_z
      real(rp) sigma_mat(6,6)
      real(rp) t6(6,6)
      real(rp), intent(in) :: zz
      logical error

      mode%z%emittance = zz * mode%sigE_E

      call transfer_matrix_calc (lat, .true., t6, ix1=0, one_turn=.TRUE.)
      t6 = pwd_mat(lat, t6, ibs_sim_params%inductance, zz)
      ! call transfer_matrix_calc_special (lat, .true., t6, ix1=0, one_turn=.true., inductance=ibs_sim_params%inductance, sig_z=zz)
      call make_smat_from_abc(t6, mode, sigma_mat, error)

      residual_pwd_sig_z = (abs(sqrt(sigma_mat(5,5)) - zz))/zz
    end function residual_pwd_sig_z

end subroutine bl_via_mat

end module ibs_mod






