program ibs_ring

use bmad
use mode3_mod
use ibs_mod
use ibs_rates_mod, only: ibs_struct
use ptc_layout_mod
use sim_utils_interface
use mode3_mod
use longitudinal_profile_mod

implicit none

type ibs_data_struct
  real(rp) current
  real(rp) a_emittance
  real(rp) b_emittance
  real(rp) sigE_E
  real(rp) mode_sig_z
  real(rp) sigma_x
  real(rp) sigma_y
  real(rp) sigma_z
end type ibs_data_struct

type(ibs_data_struct), allocatable :: ibs_data(:)
type(ele_struct) ele_at_s
type(ele_struct) pwd_ele

integer, parameter :: N_MAX_CURRENTS = 100000
real(rp) currents(N_MAX_CURRENTS)

real(rp) :: high_current
real(rp) :: a_emit, b_emit
real(rp) :: energy_spread
real(rp) :: ratio
real(rp) view_sigma_x, view_sigma_y, view_sigma_z
real(rp) :: delta_current, low_current
real(rp) :: granularity
real(rp) inductance
real(rp) eta_set, etap_set
real(rp) t6(6,6)
real(rp) sigma_mat(6,6)
real(rp) inv_Ta_int, inv_Tb_int, inv_Tz_int
real(rp) s, delta_s
real(rp) Ha, Hb
real(rp) L_ratio
real(rp) :: fake_3HC
real(rp) initial_blow_up(3)

logical error, do_pwd
logical :: ptc_calc
logical set_dispersion

character(50) in_file
character(4) :: ibs_formula
character(3)::  eqb_method
character(130) lat_file

integer :: x_view, y_view, z_view
integer i, n_steps
integer radcache
integer stdoutlun, dotinlun, scalinglun
integer emitlun
integer rateslun
integer int_rateslun
integer clog_to_use
integer status
integer ix_cache

type(ibs_struct) rates
type(normal_modes_struct) mode
type(normal_modes_struct) mode0
type(lat_struct) :: lat, lat0
type(ibs_sim_param_struct) ibs_sim_params
type(coord_struct), target, allocatable :: orb(:)
type(coord_struct) ptc_co

namelist /parameters/ &
  lat_file, &        ! Lattice file in BMAD format.
  granularity, &     ! Step size along lattice in meters.  Set to -1 for element-by-element.
  ptc_calc, &        ! Use PTC for emittance calculation.
  b_emit, &          ! Zero current vertical emittance.  Set to -1 for rad int calc.
  a_emit, &          ! Zero current horizontal emittance.  Set to -1 for rad int calc.
  energy_spread, &   ! Zero current energy spread.  Set to -1 for rad int calc.
  fake_3HC, &        ! If greater than zero, reduce rates by this factor.  IBS rates scale with 1/sigma_z
  high_current, &    ! Largest current per bunch in mA.
  delta_current, &   ! mA step size.
  low_current, &     ! Smallest current per bunch in mA.
  ibs_formula, &     ! 'cimp', 'bjmt', 'bane', 'mpzt', 'mpxx'
  clog_to_use, &     ! 1=classic, no tail cut.  2=Raubenheimer.  3=Oide, 4=Bane.  See multi_coulomb_log in ibs_mod.f90
  eqb_method, &      ! 'der' for derivatives.  'rlx' for relaxation approach.  Use 'der'.
  initial_blow_up, & ! initial starting point for 'rlx' method.
  ratio, &           ! "Coupling parameter r" hack for including coupling.
  x_view, &          ! index of element where projection is taken for horizontal beam size calculation.
  y_view, &          ! index of element where projection is taken for vertical beam size calculation.
  z_view, &          ! index of element where projection is taken for longitudinal beam size calculation.
  do_pwd, &          ! .true. or .false.: Simulate PWD by inserting element at start of lattice with 
                     ! some 
  inductance, &      ! Longitudinal inductance for PWD calc.  Effects bunch length vs. current.
  set_dispersion, &  ! If true, then apply eta_set and etap_set.  If false, then do not.
                     ! Note: this adjusts the vertical dispersion that is plugged into the IBS formulas.
                     ! This vertical dispersion is included in the sigma_x, sigma_y, sigma_z values that
                     ! emittance.dat is populated with.  See set_t6_eta for details.
  eta_set, &         ! Sets vertical dispersion at every element to this fixed value.
  etap_set           ! Sets vertical dispersion prime at every element to this fixes value.

call load_parameters_file()
if(do_pwd) then
  write(*,*) "PWD is currently disabled pending calculation improvements."
  do_pwd = .false.
endif

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
if(do_pwd) then
  call init_ele(pwd_ele)
  pwd_ele%key = taylor$
  do i=1,6
    call init_taylor_series (pwd_ele%taylor(i), 1)
    pwd_ele%taylor(i)%term(1)%coef = 1.0
    pwd_ele%taylor(i)%term(1)%expn = 0
    pwd_ele%taylor(i)%term(1)%expn(i) = 1
  enddo
  call init_taylor_series (pwd_ele%taylor(6), 2)
  pwd_ele%taylor(6)%term(1)%coef = 1.0
  pwd_ele%taylor(6)%term(1)%expn = 0
  pwd_ele%taylor(6)%term(1)%expn(6) = 1
  pwd_ele%taylor(6)%term(2)%coef = 0.0
  pwd_ele%taylor(6)%term(2)%expn = 0
  pwd_ele%taylor(6)%term(2)%expn(5) = 1
  call insert_element(lat,pwd_ele,1)
endif
bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.
IF(ptc_calc)  call lat_to_ptc_layout(lat)
call twiss_and_track(lat,orb,status)

ix_cache = -1
call radiation_integrals(lat, orb, mode, ix_cache)


if( ptc_calc ) then
  call ptc_emit_calc(lat%ele(0), mode, sigma_mat, ptc_co)
  mode%sig_z = sqrt(sigma_mat(5,5))
  mode%sige_e = sqrt(sigma_mat(6,6))
  mode%a%tune = mode%a%tune * twopi
  mode%b%tune = mode%b%tune * twopi
  mode%z%tune = mode%z%tune * twopi
  lat%a%tune = mode%a%tune
  lat%b%tune = mode%b%tune
  lat%z%tune = mode%z%tune
else
  ! associate tunes with modes
  lat%a%tune = mod(lat%ele(lat%n_ele_track)%a%phi,twopi)
  lat%b%tune = mod(lat%ele(lat%n_ele_track)%b%phi,twopi)
  call calc_z_tune(lat%branch(0))
  mode%a%tune = lat%a%tune
  mode%b%tune = lat%b%tune
  if(lat%z%tune .lt. 0) then
    mode%z%tune = twopi+lat%z%tune
  else
    mode%z%tune = lat%z%tune
  endif
endif
write(*,*) "Lattice preparation complete..."

stdoutlun = lunget()
open(stdoutlun,file='properties.out')

do i=6,stdoutlun,stdoutlun-6
  write(i,*) "Beam parameters from radiation calculation:"
  write(i,*) "   emit_a (m*rad)  : ", mode%a%emittance
  write(i,*) "   emit_b (m*rad)  : ", mode%b%emittance
  write(i,*) "   sigmaE_E (rel.) : ", mode%sigE_E
  write(i,*) "   sigma_z (m)     : ", mode%sig_z
  write(i,*)
enddo

if( b_emit .gt. 0.0 ) then
  mode%b%emittance = b_emit
endif
if( a_emit .gt. 0.0 ) then
  mode%a%emittance = a_emit
endif
if( energy_spread .gt. 0.0 ) then
  mode%sigE_E = energy_spread
endif
if( fake_3HC .gt. 0.0 ) then
  mode%sig_z = mode%sig_z * fake_3HC
endif

lat%param%n_part = 0.0d0
! without current, set_pwd_ele effectively adjusts bunch length for sigE_E
if(do_pwd) call set_pwd_ele(lat, mode, inductance)

do i=6,stdoutlun,stdoutlun-6
  write(i,*) "User-adjusted beam parameters at zero current:"
  write(i,*) "   emit_a (m*rad)  : ", mode%a%emittance
  write(i,*) "   emit_b (m*rad)  : ", mode%b%emittance
  write(i,*) "   sigmaE_E (rel.) : ", mode%sigE_E
  write(i,*) "   sigma_z (m)     : ", mode%sig_z
  write(i,*)
enddo

mode0=mode

lat%param%n_part = high_current*0.001_rp*(lat%param%total_length/c_light)/e_charge/abs(charge_of(lat%param%particle))
if(do_pwd) call set_pwd_ele(lat, mode, inductance)
mode%z%emittance = mode%sigE_E * mode%sig_z
do i=6,stdoutlun,stdoutlun-6
  write(i,*) "Beam parameters at full current with PWD, but without IBS:"
  write(i,*) "   emit_a (m*rad)  : ", mode%a%emittance
  write(i,*) "   emit_b (m*rad)  : ", mode%b%emittance
  write(i,*) "   sigmaE_E (rel.) : ", mode%sigE_E
  write(i,*) "   sigma_z (m)     : ", mode%sig_z
  write(i,*)
enddo

!populate IBS parameters structure to prepare for IBS calculation
ibs_sim_params%tau_a = lat%param%total_length / c_light / mode%a%alpha_damp
ibs_sim_params%clog_to_use = clog_to_use
ibs_sim_params%set_dispersion = set_dispersion
ibs_sim_params%eta_set = eta_set
ibs_sim_params%etap_set = etap_set
ibs_sim_params%do_pwd = do_pwd
ibs_sim_params%inductance = inductance
ibs_sim_params%formula = ibs_formula
ibs_sim_params%tau_a = lat%param%total_length / c_light / mode%a%alpha_damp  !needed for tail cut calculation

if(eqb_method == 'rlx') then
  call ibs_equib_rlx(lat,ibs_sim_params,mode0,mode,ratio,initial_blow_up,granularity)  !relaxation method
elseif(eqb_method == 'der') then
  call ibs_equib_der(lat,ibs_sim_params,mode0,mode,granularity)  !derivatives method
else
  write(*,*) "ERROR: Unrecognized setting for eqb_method: ", eqb_method
  write(*,*) "TERMINATING EXECUTION"
  stop
endif

do i=6,stdoutlun,stdoutlun-6
  write(i,*) "Beam parameters at full current with PWD and IBS:"
  write(i,*) "   emit_a (m*rad)  : ", mode%a%emittance
  write(i,*) "   emit_b (m*rad)  : ", mode%b%emittance
  write(i,*) "   sigmaE_E (rel.) : ", mode%sigE_E
  write(i,*) "   sigma_z (m)     : ", mode%sig_z
  write(i,*)
enddo

close(stdoutlun)

!Create array of currents for beam size versus current simulation
currents(:) = 0.0d0
n_steps = ceiling( (high_current-low_current) / delta_current) 
if(n_steps+1 .gt. N_MAX_CURRENTS) then
  write(*,*) "Number of currents to calculate exceeds hard-coded limit of ", N_MAX_CURRENTS
  stop
endif
do i=1,n_steps
  currents(i) = high_current - delta_current*(i-1)
enddo
! add one more current, to ensure that current at precisely low_current is calculated
n_steps = n_steps + 1
currents(n_steps) = low_current
currents = 0.001 * currents ! mA to A
ALLOCATE(ibs_data(1:n_steps))

!
!Growth rates at full current at natural emittances
!
rateslun = lunget()
open(rateslun,file='ibs_rates0.out')
write(rateslun,'(4a14)') "#s", "invTz", "invTa", "invTb"
write(rateslun,'(4a14)') "# (m)", "(1/s)", "(1/s)", "(1/s)" 
write(rateslun,'(a)') "# Element-by-element rates at start of simulation at full current."
lat%param%n_part = high_current * lat%param%total_length / e_charge / abs(charge_of(lat%param%particle)) / c_light
call ibs_rates1turn(lat,ibs_sim_params,rates,granularity)
write(rateslun,'(a,3f15.7)') "# Initial average rates: ", rates%inv_Tz, rates%inv_Ta, rates%inv_Tb
do i=1,lat%n_ele_track
  call ibs1(lat,ibs_sim_params,rates,i=i)
  write(rateslun,'(F14.5,3f16.6)') lat%ele(i)%s, rates%inv_Tz, rates%inv_Ta, rates%inv_Tb
enddo
close(rateslun)

!
!Find equilibrium emittance at each current.
!
lat0 = lat
do i=1,n_steps
  lat = lat0
  lat%param%n_part = currents(i) * lat%param%total_length / e_charge / abs(charge_of(lat%param%particle)) / c_light

  if(eqb_method == 'rlx') then
    call ibs_equib_rlx(lat,ibs_sim_params,mode0,mode,ratio,initial_blow_up,granularity)  !relaxation method
  elseif(eqb_method == 'der') then
    call ibs_equib_der(lat,ibs_sim_params,mode0,mode,granularity)  !derivatives method
  else
    write(*,*) "ERROR: Unrecognized setting for eqb_method: ", eqb_method
    write(*,*) "TERMINATING EXECUTION"
    stop
  endif

  call transfer_matrix_calc (lat, t6, ix1=x_view, one_turn=.TRUE.)
  if(set_dispersion) call set_t6_eta(t6)
  call make_smat_from_abc(t6, mode, sigma_mat, error)
  view_sigma_x = SQRT(sigma_mat(1,1))

  call transfer_matrix_calc (lat, t6, ix1=y_view, one_turn=.TRUE.)
  if(set_dispersion) call set_t6_eta(t6)
  call make_smat_from_abc(t6, mode, sigma_mat, error)
  view_sigma_y = SQRT(sigma_mat(3,3))

  call transfer_matrix_calc (lat, t6, ix1=z_view, one_turn=.TRUE.)
  if(set_dispersion) call set_t6_eta(t6)
  call make_smat_from_abc(t6, mode, sigma_mat, error)
  view_sigma_z = sqrt(sigma_mat(5,5))

  ibs_data(i) = ibs_data_struct(currents(i),mode%a%emittance,mode%b%emittance,mode%sigE_E,mode%sig_z, view_sigma_x, view_sigma_y, view_sigma_z)

  write(*,'(A,I10,A,I10,A)') "Step ", i, " of ", n_steps, " complete!"
enddo

emitlun = LUNGET()
open(emitlun, file='emittance.dat',status='replace')
write(emitlun,'(A,I0,"   ",A)') "# sigma_x calculated at ", x_view, lat%ele(x_view)%name
write(emitlun,'(A,I0,"   ",A)') "# sigma_y calculated at ", y_view, lat%ele(y_view)%name
write(emitlun,'(A,I0,"   ",A)') "# sigma_z calculated at ", z_view, lat%ele(z_view)%name
write(emitlun,'(94X,A)') "-------------Values from sigma matrix-------------"
write(emitlun,'(A14,7A18)') "# current", "emit_a", "emit_b", "sigE/E", "sig z", "sigma_x", "sigma_y", "sigmz_z"
write(emitlun,'(A14,7A18)') "#     (A)", "(m*rad)", "(m*rad)", "(rel.)", "(m)", "(m)", "(m)", "(m)"
do i=1,n_steps
  write(emitlun,"(es18.8,'   ',es15.8,'   ',es15.8,'   ',es15.8,'   'es15.8,'   ',es15.8,'   ',es15.8,'   ',es15.8)") ibs_data(i)
enddo
close(emitlun)
deallocate(ibs_data)

!
! Calculate rates for one turn at full current.
!
lat = lat0
lat%param%n_part = 0.001 * high_current * lat%param%total_length / e_charge / abs(charge_of(lat%param%particle)) / c_light
call ibs_equib_der(lat,ibs_sim_params,mode0,mode,granularity)
do i=1,lat%n_ele_track
  lat%ele(i)%a%emit = mode%a%emittance
  lat%ele(i)%b%emit = mode%b%emittance
  lat%ele(i)%z%sigma = mode%sig_z
  lat%ele(i)%z%sigma_p = mode%sigE_E
  lat%ele(i)%z%emit = mode%sig_z * mode%sigE_E
enddo
rateslun = lunget()
open(rateslun,file='ibs_rates.out')
int_rateslun = lunget()
open(int_rateslun,file='ibs_rates_integrated.out')
scalinglun = lunget()
open(scalinglun,file='ibs_scaling.out')
write(rateslun,'(a)') "# Element-by-element rates at equilibrium at full current"
write(rateslun,'(4a14)') "# s", "inv_Ta", "inv_Tb", "inv_Tz"
write(rateslun,'(4a14)') "# (m)", "(1/s)", "(1/s)", "(1/s)" 
write(int_rateslun,'(4a14)') "# s", "inv_Ta", "inv_Tb", "inv_Tz"
write(int_rateslun,'(4a14)') "# (m)", "(1/s)", "(1/s)", "(1/s)" 
write(scalinglun,'(5a14)') "# s", "inv_Ta", "inv_Tb", "Ha", "Hb"
write(scalinglun,'(5a14)') "# (m)", "(1/s)", "(1/s)", "(m)", "(m)"


inv_Ta_int = 0.0d0
inv_Tb_int = 0.0d0
inv_Tz_int = 0.0d0
delta_s = 0.1
s=delta_s
write(rateslun,*) "# delta_s = ", delta_s
write(int_rateslun,*) "# delta_s = ", delta_s

do while(s .lt. lat%param%total_length)
  call ibs1(lat,ibs_sim_params,rates,s=s)
  write(rateslun,'(F14.5,3ES14.4)') s, rates%inv_Ta, rates%inv_Tb, rates%inv_Tz

  inv_Ta_int = inv_Ta_int + rates%inv_Ta * delta_s
  inv_Tb_int = inv_Tb_int + rates%inv_Tb * delta_s
  inv_Tz_int = inv_Tz_int + rates%inv_Tz * delta_s
  write(int_rateslun,'(F14.5,3ES14.4)') s, inv_Ta_int, inv_Tb_int, inv_Tz_int

  call twiss_and_track_at_s(lat, s, ele_at_s)
  Ha = ( (ele_at_s%a%eta**2) + ( ele_at_s%a%beta*ele_at_s%a%etap + ele_at_s%a%alpha*ele_at_s%a%eta )**2)/ele_at_s%a%beta
  Hb = ( (ele_at_s%b%eta**2) + ( ele_at_s%b%beta*ele_at_s%b%etap + ele_at_s%b%alpha*ele_at_s%b%eta )**2)/ele_at_s%b%beta
  write(scalinglun,'(f14.5,4es14.4)') s, Ha/ele_at_s%b%beta/sqrt(ele_at_s%a%beta*ele_at_s%b%beta), Hb/ele_at_s%a%beta/sqrt(ele_at_s%a%beta*ele_at_s%b%beta), Ha, Hb

  s = s + delta_s
enddo
close(rateslun)
close(int_rateslun)
close(scalinglun)

!-----------------------------------------------------------
contains
  subroutine set_t6_eta(t6)
    real(rp) t6(6,6), W(6,6)
    integer i

    W = 0.0_rp
    do i=1,6
      W(i,i) = 1.0_rp
    enddo
    W(3,6) = -ibs_sim_params%eta_set
    W(4,6) = -ibs_sim_params%etap_set
    W(5,3) =  ibs_sim_params%etap_set
    W(5,4) = -ibs_sim_params%eta_set

    t6 = matmul(t6,W)
  end subroutine set_t6_eta

!-----------------------------------------------------------
! contains

  subroutine load_parameters_file()
    !-Set bogus values for namelist parameters, so we can check that they were
    !-set by the .in file.
    lat_file      = ''
    b_emit        = -99.0
    a_emit        = -99.0
    energy_spread = -99.0
    fake_3HC      = -99.0
    high_current  = -99.0
    ibs_formula   = ''
    ratio         = -99.0
    delta_current      = -99.0
    low_current       = -99.0
    x_view        = 0
    y_view        = 0
    z_view        = 0
    granularity   = -1.0
    inductance    = -99.0
    clog_to_use   = -99
    eqb_method    = ''
    initial_blow_up = (/-1.0d0,-1.0d0,-1.0d0/)

    eta_set = -99.0
    etap_set = -99.0
    ptc_calc      = .false.
    do_pwd = .false.

    dotinlun = LUNGET()
    call getarg(1,in_file)
    open(dotinlun,file=in_file,status='old')
    read(dotinlun,nml=parameters)
    close(dotinlun)

    !-Check if any parameters were missing from the .in file.
    if( lat_file == '' ) call param_bomb('lat_file')
    if( b_emit .lt. -90 ) call param_bomb('b_emit')
    if( a_emit .lt. -90 ) call param_bomb('a_emit')
    if( fake_3HC .lt. -90 ) call param_bomb('fake_3HC')
    if( energy_spread .lt. -90 ) call param_bomb('energy_spread')
    if( high_current .lt. -90 ) call param_bomb('high_current')
    if( ibs_formula == '' ) call param_bomb('ibs_formula')
    if( ratio .lt. -90 ) call param_bomb('ratio')
    if( delta_current .lt. -90 ) call param_bomb('delta_current')
    if( low_current .lt. -90 ) call param_bomb('low_current')
    if( x_view .lt. -90 ) call param_bomb('x_view')
    if( y_view .lt. -90 ) call param_bomb('y_view')
    if( z_view .lt. -90 ) call param_bomb('z_view')
    if( eta_set .lt. -90 ) call param_bomb('eta_set')
    if( etap_set .lt. -90 ) call param_bomb('etap_set')
    if( granularity .lt. -90 ) call param_bomb('granularity')
    if( inductance .lt. -90 ) call param_bomb('inductance')
    if( clog_to_use .lt. -90 ) call param_bomb('clog_to_use')
    if( eqb_method == '' ) call param_bomb('eqb_method')
    if( any(initial_blow_up .lt. 0) ) call param_bomb('initial_blow_up')
  end subroutine load_parameters_file

end program ibs_ring

subroutine param_bomb(parameter)
  character(*) parameter
  write(*,*) "ERROR: parameter '", parameter, "' is missing from the .in file."
  write(*,*) "TERMINATING EXECUTION"
  stop
end subroutine param_bomb


