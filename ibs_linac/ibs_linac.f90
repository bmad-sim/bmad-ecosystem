!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!+
! program ibs_linac
!
! Simple program to calculate IBS in a linac.  
!
! Usage: ibs_linac (looks for ibs_linac.in file)
!        ibs_linac file.in
!
!-


program ibs_linac

use beam_mod
use twiss_and_track_mod
use ibs_mod
use random_mod

implicit none

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (coord_struct) :: orb0
type (coord_struct), allocatable :: orbit(:)
type (rad_int_all_ele_struct), target :: rad_int_by_ele
type (normal_modes_struct) :: normal_mode
type (rad_int1_struct), pointer :: rad_int1
type (beam_init_struct) :: beam_init
type (beam_struct), target :: beam
type (bunch_struct), pointer :: bunch
type (bunch_params_struct) :: bunch_params
type (ibs_sim_param_struct) :: ibs_sim_params

real(rp) :: delta_sigma_energy, delta_emit_a, delta_emit_b, rad_delta_eV2, gamma
real(rp), parameter :: c_q = 3.84e-13
real(rp) :: initial_slice_energy_spread_eV, energy_spread_eV, sigma_z

integer :: ix, status, ix_branch
integer :: namelist_file, n_char

character(100) :: lat_name, lat_path, base_name, in_file
character(30), parameter :: r_name = 'ibs_linac'
character(4) :: ibs_formula

logical :: radiation_damping_on, radiation_fluctuations_on, ISR_energy_spread_on
logical :: err, use_beam, verbose, ibs_affects_bunch, ibs_on

namelist / ibs_linac_params / &
    lat_name, ibs_formula, ix_branch, &
    use_beam, beam_init, radiation_damping_on, radiation_fluctuations_on, &
    ISR_energy_spread_on, ibs_affects_bunch, ibs_on, &
    initial_slice_energy_spread_eV, verbose

!------------------------------------------
!Defaults for namelist
lat_name = 'lat.bmad'
ibs_formula = 'bjmt'   ! See ibs_mod
use_beam = .true.      ! Beam is tracked to set sigma_z, emittances for every element 
initial_slice_energy_spread_eV = 5e3 ! Slice energy spread is kept track of separate from the bunch
ISR_energy_spread_on = .true.      ! For energy spread calc
radiation_damping_on = .false.     ! For bunch tracking
radiation_fluctuations_on = .true. ! For bunch tracking
verbose = .false. 
beam_init%n_bunch = 1
ibs_affects_bunch = .false.
ibs_on = .true.
ix_branch = 0

!Read namelist
in_file = 'ibs_linac.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)

namelist_file = lunget()
print *, 'Opening: ', trim(in_file)
open (namelist_file, file = in_file, status = "old")
read (namelist_file, nml = ibs_linac_params)
close (namelist_file)

! Set params
ibs_sim_params%tau_a = 0            ! No damping time in a linac
ibs_sim_params%clog_to_use = 1      ! This disables the tail-cut in 'kubo'
ibs_sim_params%formula = ibs_formula

!Trim filename
n_char= SplitFileName(lat_name, lat_path, base_name) 

!Parse Lattice
call bmad_parser (lat_name, lat)
branch => lat%branch(ix_branch)

!--------------------
! Radiation 
bmad_com%radiation_damping_on = radiation_damping_on
bmad_com%radiation_fluctuations_on = radiation_fluctuations_on
print *, 'bmad_com%radiation_damping_on: ', bmad_com%radiation_damping_on
print *, 'bmad_com%radiation_fluctuations_on: ', bmad_com%radiation_fluctuations_on

call twiss_and_track(lat, orbit, status)
if (status /= ok$) then
  print *, 'problem with twiss_and_track'
  stop
endif

if (verbose) print *, 'radiation_integrals'
call radiation_integrals (lat, orbit, normal_mode, rad_int_by_ele = rad_int_by_ele)

! Set element 0
ele => branch%ele(0)

call init_beam_distribution (ele, lat%param, beam_init, beam)

if (upcase(beam_init%distribution_type(1)) == 'FILE') then
  write (*, '(2a)') 'Using particles from beam_file: ', beam_init%position_file
endif

bunch => beam%bunch(1)

call calc_bunch_params (bunch, bunch_params, err, print_err = .true.)

! Set running energy spread and sigma_z
sigma_z = sqrt(bunch_params%sigma(5,5))
energy_spread_eV = initial_slice_energy_spread_eV !Old: sqrt(bunch_params%sigma(s66$))*ele%value(e_tot$)

if (verbose) print *, 'Bunch initialized'

! Set first element
call set_ele(lat%ele(0), bunch_params, energy_spread_eV)

! Set 
lat%param%n_part = bunch_params%charge_live / e_charge
write (*, '(a)')           'Beginning bunch:'
write (*, '(a, f15.7, a)') '      energy            : ', 1e-6_rp*ele%value(e_tot$), ' MeV'
write (*, '(a, f15.7, a)') '      charge            : ', 1e12_rp*bunch_params%charge_live, ' pC'
write (*, '(a, f15.7, a)') '      norm_emit_a       : ', 1e6_rp*bunch_params%a%norm_emit, ' mm-mrad'
write (*, '(a, f15.7, a)') '      norm_emit_b       : ', 1e6_rp*bunch_params%b%norm_emit, ' mm-mrad'
write (*, '(a, f15.7, a)') '      sigma_z/c         : ', 1e12_rp*sqrt(bunch_params%sigma(5,5))/c_light, ' ps'
write (*, '(a, f15.7, a)') '      slice sigma_E     : ', energy_spread_eV, ' eV'
if (use_beam) then
  write (*,'(a)')  'Using beam for sigma_z calc'
else
  write(*,'(a)')   'Not tracking beam'
endif
if (ibs_on) then
  write (*, '(2a)') 'Using IBS formula: ', ibs_sim_params%formula
else
  delta_sigma_energy = 0
  delta_emit_a = 0
  delta_emit_b = 0
endif
if (ISR_energy_spread_on) write (*, '(a)') 'ISR included in energy spread calculation'
if(ibs_affects_bunch) write (*, '(a)') 'IBS will affect the bunch energy spread'


! IBS loop

write (*, '(a)') 'BEGIN_DATA'
write (*, '(9a15)') 's', 'gamma', 'norm_emit_a', 'norm_emit_b ', 'sigma_z/c ', 'sigma_E', 'demit_a', 'demit_b', 'dsigma_E'
write (*, '(9a15)') 'm', '1', 'mm-mrad', 'mm-mrad', 'fs', 'eV', 'mm-mrad', 'mm-mrad', 'eV'

do ix=1, lat%n_ele_track
  ele => lat%ele(ix)
  rad_int1 => rad_int_by_ele%branch(ix_branch)%ele(ix)
  call convert_total_energy_to(ele%value(e_tot$), lat%param%particle, gamma)
  
  if (use_beam) then
    call track1_bunch (bunch, ele, err)
    if (err) then
      print *, 'Bunch tracking error in ele: ', trim(ele%name)
      stop
    endif
    ! Don't worry about calc_bunc_params errors
    call calc_bunch_params (bunch, bunch_params, err, print_err = .false.)
  endif  
  
  ! Take into account bunch compression by magnifying the previous energy spread by the compression ratio
  !  (using previous sigma_z, then set sigma_z
  energy_spread_eV = energy_spread_eV * sigma_z /  sqrt(bunch_params%sigma(5,5)) 
  sigma_z = sqrt(bunch_params%sigma(5,5)) 
  
  ! Set ele's emittances and sigma_z from the bunch. The energy spread is set separately
  call set_ele(ele, bunch_params, energy_spread_eV)
  
  ! Radiation integral contribution
  if (ISR_energy_spread_on) then
    rad_delta_eV2 = (4./3.)*c_q*r_e*(m_electron**2)*rad_int1%lin_i3_E7
  else
    rad_delta_eV2 = 0
  endif
  
  ! IBS deltas
  if (ibs_on) call ibs_delta_calc(lat, ix, ibs_sim_params, sigma_mat = bunch_params%sigma, & ! sigma_mat only used for 'kubo' formula
  	delta_emit_a = delta_emit_a, &
  	delta_emit_b = delta_emit_b, &
  	delta_sigma_energy = delta_sigma_energy)
  !if (verbose) write (*, '(a20, 3es15.7)') trim(ele%name), delta_emit_a, delta_emit_b, delta_sigma_energy
  
  ! Spread bunch energy
  if (ibs_affects_bunch .and. ibs_on .and. (delta_sigma_energy > 0) ) then
    call spread_energy( sqrt(2*(delta_sigma_energy/ele%value(e_tot$)) * ele%z%sigma_p) )
  endif
  
  energy_spread_eV = sqrt( (energy_spread_eV + delta_sigma_energy)**2 + rad_delta_ev2)
  ele%a%emit = ele%a%emit + delta_emit_a
  ele%b%emit = ele%b%emit + delta_emit_b
  ele%z%sigma_p = energy_spread_eV/ele%value(e_tot$)
  
   write (*, '(5f15.7, 10es15.7)') ele%s, gamma, 1e6_rp*ele%a%emit*gamma, 1e6_rp*ele%b%emit*gamma, 1e15_rp*ele%z%sigma/c_light, energy_spread_eV, &
     1e6_rp*gamma*delta_emit_a, 1e6_rp*gamma*delta_emit_b, delta_sigma_energy
  
enddo
write (*, '(a)') 'END_DATA'


contains


subroutine set_ele(ele, bunch_params, energy_spread_eV)
implicit none
type (ele_struct) :: ele
type (bunch_params_struct) :: bunch_params
real (rp) :: gamma, energy_spread_eV
!
call convert_total_energy_to(ele%value(e_tot$), lat%param%particle, gamma)
ele%z%sigma   = sqrt(bunch_params%sigma(5,5))
ele%z%sigma_p = energy_spread_eV / ele%value(e_tot$)
ele%a%emit    = bunch_params%a%norm_emit / gamma
ele%b%emit    = bunch_params%b%norm_emit / gamma
end subroutine

subroutine spread_energy(sigma_delta)
implicit none 
real(rp) ::this_ran, sigma_delta
integer :: i
do i = 1, size(bunch%particle)
  call ran_gauss (this_ran)
  bunch%particle(i)%vec(6) = bunch%particle(i)%vec(6) + sigma_delta*this_ran
enddo 

end subroutine

end program
