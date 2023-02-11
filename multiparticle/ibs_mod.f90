module ibs_mod

use mode3_mod
use ibs_rates_mod
use twiss_and_track_mod
use longitudinal_profile_mod
use fgsl
use, intrinsic :: iso_c_binding

type ibs_sim_param_struct
  real(rp) :: tau_a = 0.0d0          ! horizontal damping rate (needed for coulomb log tail cut)
  integer :: clog_to_use = 1      ! see multi_coulomb_log subroutine for valid settings.  Set to 1 to disable tail-cut.  Set to 1 for linacs.
  logical :: set_dispersion = .false.   ! True: add vertical dispersion to transfer matrix.  Valid for kubo method.
  real(rp) :: eta_set = 0.0d0         ! If set_dispersion, then this value is used to add y-z coupling to the transfer matrix.
  real(rp) :: etap_set = 0.0d0        ! If set_dispersion, then this value is used to add y-z coupling to the transfer matrix.
  logical :: do_pwd = .false.           ! If true, then use potential well distortion to calculate bunch lengths.  If false,
                           ! bunch length is proportional to energy spread.
  real(rp) :: inductance = 0.0d0      ! Inductive part of impedance for pwd calc.  
  character(4) :: formula = 'bjmt'     ! Which IBS formulation to use.  See subroutine ibs1 for a list.
  !real(rp) :: fake_3HC = -1   ! If greater than zero, divide growth rates by this factor.
end type

type ibs_lifetime_struct  ! Beam lifetime based on IBS.  Useful for linacs.  These quantities are populated with time
                          ! required for beam size to increase by some amount.
  real(rp) Tlx
  real(rp) Tly
  real(rp) Tlp
end type

type ibs_maxratio_struct  ! Parameters for IBS lifetime calculation
  real(rp) rx
  real(rp) ry
  real(rp) r_p
end type

real(fgsl_double), parameter :: eps_7 = 1.0d-7
integer(fgsl_size_t), parameter :: space_limit = 1000_fgsl_size_t

contains

!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------
!+
! Subroutine ibs_equib_rlx(lat,ibs_sim_params,inmode,ibsmode,ratio,initial_blow_up,granularity)
! Iterates to equilibrium beam conditions using relaxation method
!
! This method requires that the initial beam size be larger than the equilibrium beam size.
! An initial_blow_up of 3 to 5 is a good place to start.
!
! See ibs_rates subroutine for available IBS rate formulas.
!
! Input:
!   lat             -- lat_struct: lattice for tracking
!     %param%n_part  -- Real: number of particles in bunch
!   ibs_sim_params   -- ibs_sim_param_struct: parameters for IBS calculation
!   inmode           -- normal_modes_struct: natural beam parameters 
!   ratio            -- Real: Ratio of vert_emit_coupling / vert_emit_total
!   initial_blow_up  -- Real: Factor multiplied to all thre bunch dimensions prior to starting iteration.
!   granularity      -- Real: Step size for slicing lattice.  i.e. set to 1 to calculate IBS rates every 1 meter.
!
! Output:
!   ibsmode          -- normal_modes_struct: beam parameters after IBS effects
!-

subroutine ibs_equib_rlx(lat,ibs_sim_params,inmode,ibsmode,ratio,initial_blow_up,granularity)

implicit none
                                                                                                       
type(lat_struct) :: lat
type(ibs_sim_param_struct) ibs_sim_params
type(normal_modes_struct) :: inmode
type(normal_modes_struct) :: ibsmode
real(rp) :: granularity
real(rp) :: ratio
real(rp) :: initial_blow_up(3)
type(ibs_struct) rates

real(rp) time_for_one_turn
real(rp) tau_a, tau_b, tau_z
real(rp) Ta, Tb, Tz
real(rp) afactor,bfactor,zfactor
real(rp) emit_a, emit_b
real(rp) emit_a0, emit_b0
real(rp) advance, threshold
real(rp) sigE_E, sigE_E0, sigma_z0, L_ratio
real(rp) ratio_a, ratio_b, ratio_z
real(rp) emit_a_prev, emit_b_prev, sigE_E_prev
real(rp) running_emit_x(1:40), running_emit_y(1:40), running_sigE_E(1:40)
integer half
real(rp) runavg_emit_x_A, runavg_emit_y_A, runavg_sigE_E_A
real(rp) runavg_emit_x_B, runavg_emit_y_B, runavg_sigE_E_B
real(rp) energy, gamma, KE, rbeta
logical converged
integer counter, i
integer ibs_mod_lun

character(20) :: r_name = 'ibs_equib_rlx'

emit_a0 = inmode%a%emittance
emit_b0 = inmode%b%emittance
sigE_E0 = inmode%sigE_E

!compute the SR betatron damping times
time_for_one_turn = lat%param%total_length / c_light
tau_a = time_for_one_turn / inmode%a%alpha_damp
tau_b = time_for_one_turn / inmode%b%alpha_damp
tau_z = time_for_one_turn / inmode%z%alpha_damp 

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

ibsmode = inmode
L_ratio = inmode%sig_z / inmode%sigE_E
ibsmode%a%emittance = inmode%a%emittance * initial_blow_up(1)
ibsmode%b%emittance = inmode%b%emittance * initial_blow_up(2)
ibsmode%sigE_E = inmode%sigE_E * sqrt(initial_blow_up(3))
ibsmode%sig_z = inmode%sig_z * sqrt(initial_blow_up(3))
ibsmode%z%emittance = ibsmode%sig_z * ibsmode%sigE_E

do i=1, SIZE(running_emit_x)
  running_emit_x(i) = 0.0
  running_emit_y(i) = 0.0
  running_sigE_E(i) = 0.0
enddo
half = SIZE(running_emit_x) / 2

ibs_mod_lun = lunget()
open(ibs_mod_lun,file="ibs_mod.progress")
do WHILE(.not.converged)
  counter = counter + 1
  do i=1, lat%n_ele_track
    lat%ele(i)%a%emit = ibsmode%a%emittance
    lat%ele(i)%b%emit = ibsmode%b%emittance
    lat%ele(i)%z%sigma = ibsmode%sig_z
    lat%ele(i)%z%sigma_p = ibsmode%sigE_E
    lat%ele(i)%z%emit = ibsmode%sig_z * ibsmode%sigE_E
  enddo

  call ibs_rates1turn(lat,ibs_sim_params,rates,granularity)

  ratio_a = tau_a*rates%inv_Ta
  ratio_b = tau_b*rates%inv_Tb
  ratio_z = tau_z*rates%inv_Tz

  !Initial blow up must be such that the first time through this loop, the following conditional
  !will be false (i.e. the else will be taken).

  if( (ratio_a .ge. 1) .or. (ratio_b .ge. 1) .or. (ratio_z .ge. 1) ) then
    write(*,*) "Negative emittance encountered.  Reducing step size by half."
    advance = advance / 2.0d0
    ibsmode%a%emittance = emit_a_prev + advance*( afactor*emit_a0 - emit_a_prev ) 
    ibsmode%b%emittance = emit_b_prev + advance*( ((1.0-ratio)*bfactor+ratio*afactor)*emit_b0 - emit_b_prev )
    ibsmode%sigE_E      = sigE_E_prev + advance*( zfactor*sigE_E0 - sigE_E_prev )
  else
    !advance = advance + (advance0-advance)*0.1
    if(mod(counter,500)==0) then
      write(*,'(a,i10)') "Still working.  See ibs_mod.progress to monitor progress.  Steps: ", counter 
    endif

    emit_a_prev = ibsmode%a%emittance
    emit_b_prev = ibsmode%b%emittance
    sigE_E_prev = ibsmode%sigE_E

    afactor = 1.0/(1.0-ratio_a)
    bfactor = 1.0/(1.0-ratio_b)
    zfactor = 1.0/(1.0-ratio_z)

    emit_a = ibsmode%a%emittance + advance*( afactor*emit_a0 - ibsmode%a%emittance ) 
    emit_b = ibsmode%b%emittance + advance*( ((1.0-ratio)*bfactor+ratio*afactor)*emit_b0 - ibsmode%b%emittance )
    sigE_E = ibsmode%sigE_E      + advance*( zfactor*sigE_E0 - ibsmode%sigE_E )

    write(ibs_mod_lun,'(3es14.5,es14.5)') emit_a, emit_b, sigE_E

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
      ibsmode%sig_z = L_ratio * ibsmode%sigE_E
      ibsmode%z%emittance = ibsmode%sig_z * ibsmode%sigE_E 
    endif
  endif
enddo
close(ibs_mod_lun)
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
type(normal_modes_struct) :: inmode
type(normal_modes_struct) :: ibsmode
real(rp) :: granularity
type(ibs_struct) rates

real(rp) time_for_one_turn
real(rp) tau_a, tau_b, tau_z
real(rp) Ta, Tb, Tz
real(rp) emit_a, emit_b, sigE_E
real(rp) emit_a0, emit_b0, sigE_E0
real(rp) dadt, dbdt, dsEdt, dT
real(rp) threshold
real(rp) sigma_z0, L_ratio
real(rp) sig_z

logical converged
integer i

!compute the SR betatron damping times
time_for_one_turn = lat%param%total_length / c_light
tau_a = time_for_one_turn / inmode%a%alpha_damp
tau_b = time_for_one_turn / inmode%b%alpha_damp
tau_z = time_for_one_turn / inmode%z%alpha_damp

sigma_z0 = inmode%sig_z
sigE_E0 = inmode%sigE_E
emit_a0 = inmode%a%emittance
emit_b0 = inmode%b%emittance
L_ratio = sigma_z0 / sigE_E0

threshold = .0000001 !fractional changes in emittance smaller than this
                   !indicate convergence
converged = .false.
dT = tau_a / 10.0 !Time to advance per iteration

ibsmode = inmode
ibsmode%sig_z = L_ratio * ibsmode%sigE_E
ibsmode%z%emittance = ibsmode%sig_z * ibsmode%sigE_E

do WHILE(.not.converged)
  do i=1, lat%n_ele_track
    lat%ele(i)%a%emit = ibsmode%a%emittance
    lat%ele(i)%b%emit = ibsmode%b%emittance
    lat%ele(i)%z%emit = ibsmode%z%emittance
    lat%ele(i)%z%sigma = ibsmode%sig_z
    lat%ele(i)%z%sigma_p = ibsmode%sigE_E
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
  ! tau_a,b and Ta,b are betatron growth times ... the factors of 2 convert them to emittance times, so that dadt and dbdt are
  ! emittance growth rates.  tau_z is a time rate of change in sigE_E.
  dadt = -(emit_a-emit_a0)*2./tau_a + emit_a*2./Ta
  dbdt = -(emit_b-emit_b0)*2./tau_b + emit_b*2./Tb 
  dsEdt= -(sigE_E-sigE_E0)/tau_z + sigE_E/Tz

  if( abs(dadt*dT)/emit_a .lt. threshold ) then
    if( abs(dbdt*dT)/emit_b .lt. threshold ) then
      if( abs(dsEdt*dT)/ibsmode%sigE_E .lt. threshold ) then
        converged = .true.
      endif
    endif
  endif

  ibsmode%a%emittance = emit_a + dadt*dT
  ibsmode%b%emittance = emit_b + dbdt*dT
  ibsmode%sigE_E = ibsmode%sigE_E + dsEdt*dT
  ibsmode%sig_z = L_ratio * ibsmode%sigE_E
  ibsmode%z%emittance = ibsmode%sig_z * ibsmode%sigE_E 

  !write(*,'(3es15.6)') ibsmode%a%emittance, ibsmode%b%emittance, ibsmode%z%emittance
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

call ibs1(lat, ibs_sim_params, rates1ele, i=ix)
if (present(delta_sigma_energy)) delta_sigma_energy = ele%value(l$)/c_light*rates1ele%inv_Tz*ele%z%sigma_p*ele%value(E_TOT$)
if (present(delta_emit_a))       delta_emit_a     = 2*ele%value(l$)/c_light*rates1ele%inv_Ta*ele%a%emit
if (present(delta_emit_b))       delta_emit_b     = 2*ele%value(l$)/c_light*rates1ele%inv_Tb*ele%b%emit
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

    call convert_total_energy_to(lat%ele(i)%value(E_TOT$), lat%param%particle, gamma1)
    call convert_total_energy_to(lat%ele(i+1)%value(E_TOT$), lat%param%particle, gamma2)
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
real(rp) coulomb_log
real(rp) n_part

n_part = lat%param%n_part

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

call multi_coulomb_log(ibs_sim_params, ele, coulomb_log, n_part)

if(ibs_sim_params%formula == 'cimp') then
  call cimp1(ele, coulomb_log, rates, n_part)
elseif(ibs_sim_params%formula == 'bjmt') then
  call bjmt1(ele, coulomb_log, rates, n_part)
elseif(ibs_sim_params%formula == 'bane') then
  call bane1(ele, coulomb_log, rates, n_part)
elseif(ibs_sim_params%formula == 'mpzt') then
  call mpzt1(ele, coulomb_log, rates, n_part)
elseif(ibs_sim_params%formula == 'mpxx') then
  call mpxx1(ele, coulomb_log, rates, n_part)
else
  write(*,*) "Invalid IBS formula selected ... terminating"
  write(*,*) "Valid IBS formulas are: cimp, bjmt, bane, mpzt, and mpxx"
  STOP
endif

!replaced by fake_3HC in driver ! IBS theory gives growth rates that go exactly as 1/sigma_z.  Simulate effect of 3rd harmonic
!replaced by fake_3HC in driver ! cavity by dividing rates by some factor.
!replaced by fake_3HC in driver if( ibs_sim_params%fake_3HC .gt. 0 ) then
!replaced by fake_3HC in driver   rates%inv_Tz = rates%inv_Tz / ibs_sim_params%fake_3HC
!replaced by fake_3HC in driver   rates%inv_Ta = rates%inv_Ta / ibs_sim_params%fake_3HC
!replaced by fake_3HC in driver   rates%inv_Tb = rates%inv_Tb / ibs_sim_params%fake_3HC
!replaced by fake_3HC in driver endif
end subroutine ibs1

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
real(rp) classical_radius

!fgsl variables
type(c_ptr) :: ptr
real(fgsl_double), target :: args(1:3)
type(fgsl_integration_workspace) :: integ_wk
type(fgsl_function) :: integrand_ready
integer(fgsl_int) :: fgsl_status
real(fgsl_double) :: integration_result
real(fgsl_double) :: abserr

energy = ele%value(E_TOT$)
call convert_total_energy_to(energy, ele%branch%lat%param%particle, gamma=gamma)

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

classical_radius = c_light*c_light*e_charge*1.0d-7*abs(charge_of(ele%branch%lat%param%particle))/mass_of(ele%branch%lat%param%particle)

if( ibs_sim_params%clog_to_use == 1 ) then
  !Classic Coulomb Log.
  coulomb_log = log( (gamma**2)*sigma_b*emit_a/classical_radius/beta_a )
elseif( ibs_sim_params%clog_to_use == 2 ) then
  !Tail cut Raubenheimer gstar (ignores vertical dispersion !)
  Ha = ( Da**2 + (beta_a*Dap + alpha_a*Da)**2 ) / beta_a
  g2 = gamma*gamma
  sp2 = sigma_p*sigma_p
  u = g2*( Ha/emit_a + 1.0/sp2 + beta_a/g2/emit_a + beta_b/g2/emit_b )
  v = g2*( Ha*beta_b/emit_a/emit_b + beta_a*beta_b/g2/emit_a/emit_b + beta_a/sp2/emit_a + beta_b/sp2/emit_b + Da*Da/emit_a/emit_a ) 
  w = g2*( Da*Da*beta_b/emit_a/emit_a/emit_b + beta_a*beta_b/sp2/emit_a/emit_b )

  Bbar = gamma*sqrt( gamma_a*emit_a + gamma_b*emit_b + sp2/g2)
  qmin = 2.0 * classical_radius / sigma_b / Bbar

  !FGSL integration
  ptr = c_loc(args)
  args = (/u,v,w/)
  integ_wk = fgsl_integration_workspace_alloc(space_limit)
  integrand_ready = fgsl_function_init(rclog_integrand, ptr)
  fgsl_status = fgsl_integration_qag(integrand_ready, 0.0d0, 40.0d0, eps_7, eps_7, &
                                         space_limit, 3, integ_wk, integration_result, abserr)

  call fgsl_function_free(integrand_ready)
  call fgsl_integration_workspace_free(integ_wk)

  qmax = sqrt( ibs_sim_params%tau_a*n_part*c_light*classical_radius**2/4.0d0/pi/g2/emit_a/emit_b/sigma_z/sigma_p * integration_result )

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
!
! Output:
!   sigma_z       -- real(rp): Bunch length. FWHM/TwoRootTwoLogTwo from bunch profile
!-
subroutine bl_via_vlassov(current,alpha,Energy,sigma_p,Vrf,omega,U0,circ,R,L,sigma_z)

implicit none

real(rp) current, alpha, Energy, sigma_p, Vrf, omega, U0, circ, R, L
real(rp) sigma_z, delta_e, A, Q, T0, phi
real(rp), parameter :: bound = 500.0d-12
real(rp) args(1:8)

delta_e = sigma_p * Energy
T0 = circ/c_light
Q = current * T0
phi = -acos(U0/Vrf)
A = Energy/delta_e/delta_e/alpha/T0

args = [A, Vrf, Q, omega, phi, R, L, U0]

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

use super_recipes_mod, only: super_brent

implicit none

type(lat_struct) lat
type(ibs_sim_param_struct) ibs_sim_params
type(normal_modes_struct) mode
real(rp) sig_z
real(rp) lb, mb, ub, min_val
integer status

!

lb = 0.90 * mode%sig_z
mb = mode%sig_z
ub = 1.50 * mode%sig_z

min_val = super_brent(lb,mb,ub, residual_pwd_sig_z, 1.0e-5_rp, 1e-18_rp*abs(lb), sig_z, status)

!------------------------------------------------------------------
contains

function residual_pwd_sig_z(zz, status)
real(rp) residual_pwd_sig_z
real(rp) sigma_mat(6,6)
real(rp) t6(6,6)
real(rp), intent(in) :: zz
integer, optional :: status
logical error

!

mode%z%emittance = zz * mode%sigE_E

call transfer_matrix_calc (lat, t6, ix1=0, one_turn=.true.)
t6 = pwd_mat(lat, t6, ibs_sim_params%inductance, zz)
! call transfer_matrix_calc_special (lat, t6, ix1=0, one_turn=.true., inductance=ibs_sim_params%inductance, sig_z=zz)
call make_smat_from_abc(t6, mode, sigma_mat, error)

residual_pwd_sig_z = (abs(sqrt(sigma_mat(5,5)) - zz))/zz
end function residual_pwd_sig_z

end subroutine bl_via_mat

end module ibs_mod






