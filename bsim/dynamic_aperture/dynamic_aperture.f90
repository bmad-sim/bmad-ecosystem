!+
! Program dynamic_aperture
!
! Program for maping out the dynamic aperture.
! Not to be confused with the CESR centric dynamic_aperture_cesr program.
!-

program dynamic_aperture_program

use da_program_mod

!$ use omp_lib

implicit none

type (lat_struct), target :: lat
type (aperture_param_struct) :: da_param

type (ltt_params_struct), target :: ltt
type (ltt_com_struct), target :: ltt_com

type (aperture_point_struct), pointer :: da_point
type (aperture_scan_struct), allocatable, target :: aperture_scan(:)
type (aperture_scan_struct), pointer :: da
type (coord_struct), allocatable :: closed_orb(:)
type (ele_pointer_struct), allocatable :: eles(:)
type (ele_struct), pointer :: ele
type (branch_struct), pointer :: branch

real(rp) dpz(20)
real(rp) :: ramping_start_time = 0
integer nargs, ios, i, j, n_dpz, nt

logical :: ramping_on = .false., set_rf_off = .false.
logical err

character(160) :: in_file, lat_file = '', dat_file, gnu_command

namelist / params / bmad_com, ltt, da_param, set_rf_off, dpz, dat_file, &
            ramping_start_time, lat_file, ramping_on

!

track1_preprocess_ptr => ltt_track1_preprocess
track1_bunch_hook_ptr => ltt_track1_bunch_hook
track_many_hook_ptr   => track_many_hook

! Set inits

bmad_com%auto_bookkeeper = .false.   ! Makes tracking faster
bmad_com%absolute_time_tracking = .true.
dpz = real_garbage$

! Read parameters
! Read the master input file again after bmad_parser is called so that bmad_com parameters
! set in the file take precedence over bmad_com parameters set in the lattice file.

in_file = 'dynamic_aperture.init'

nargs = command_argument_count()
if (nargs == 1) then
  call get_command_argument(1, in_file)
elseif (nargs > 1) then
  print *, 'ERROR: EXTRA STUFF ON COMMAND LINE'
  stop
endif

print *, 'Parameter input file: ', trim(in_file)

open (unit = 1, file = in_file, status = 'old', iostat = ios)

if (ios /= 0) then
  print *, 'ERROR: CANNOT OPEN FILE: ', trim(in_file)
  stop
endif

read (1, nml = params)

ltt_params_global => ltt
ltt_com_global => ltt_com
if (ramping_on)              ltt%ramping_on = ramping_on                   ! Old style
if (ramping_start_time /= 0) ltt%ramping_start_time = ramping_start_time   ! Old style
if (lat_file /= '')          ltt%lat_file = lat_file                       ! Old style
if (set_rf_off)              ltt%rfcavity_on = .false.                     ! Old style
ltt%tracking_method = 'OLD'

call bmad_parser (ltt%lat_file, ltt_com%lat)

rewind(1)
read (1, nml = params)

close (1)

if (ltt%ele_start /= '' .and. da_param%start_ele /= '') then
  print *, 'SETTING BOTH LTT%ELE_START *AND* DA_PARAM%START_ELE NOT PERMITTED!'
  stop
endif

if (ltt%ele_start /= '') da_param%start_ele = ltt%ele_start
if (da_param%start_ele /= '') ltt%ele_start = da_param%start_ele

!

if (abs(da_param%max_angle) > 7 .or. abs(da_param%min_angle) > 7) then
  print *, '!!! da_param%min_angle and da_param%max_angle are not in radians instead of '
  print *, "!!! degrees to be consistent with Tao's implementation of dynamic aperture."
  print *, "!!! Please modify the init file appropriately!!!"
  stop
endif

!da_param%min_angle = da_param%min_angle * pi / 180
!da_param%max_angle = da_param%max_angle * pi / 180

! Print number of threads

nt = 0
!$ nt = omp_get_max_threads()
!$ print '(a, i0)', 'Note: Maximum number of OpenMP threads available: ', nt
if (nt == 0) print '(a)', 'Note: Program has been compiled without OpenMP threading.'

ltt%simulation_mode = 'SINGLE'
ltt_com%track_bypass = .true.
call ltt_init_params(ltt, ltt_com)
call ltt_init_tracking(ltt, ltt_com)
ltt_com%track_bypass = .false.
ltt%ramp_update_each_particle = .true.    ! Needed since doing single particle tracking
branch => ltt_com%tracking_lat%branch(ltt_com%ix_branch)

! Read in lattice

n_dpz = count(dpz /= real_garbage$)
print *, 'Note: Number of dpz points: ', n_dpz

if (.not. bmad_com%absolute_time_tracking) then
  print *, 'Note: absolute time tracking is OFF!'
endif

print *, 'Data file: ', trim(dat_file)

write (gnu_command, '(a, i0, 3a)') 'plot for [IDX=1:', n_dpz, '] "', &
                  trim(dat_file), '" index (IDX-1) u 1:2 w lines title columnheader(1)'

! Scan

open (1, file = dat_file)
write (1, '(2a)') '# master_input_file                  = ' // quote(ltt_com%master_input_file)
write (1, '(2a)') '# ltt%lat_file                       = ' // quote(ltt%lat_file)
write (1, '(2a)') '# ltt%tracking_method                = ' // quote(ltt%tracking_method)
write (1, '(2a)') '# ltt%ele_start                      = ' // quote(ltt%ele_start)
if (ltt%tracking_method == 'MAP' .or. ltt%simulation_mode == 'CHECK') then
  write (1, '(2a)') '# ltt%map_order                      = ' // int_str(ltt%map_order)
  write (1, '(2a)') '# ltt%exclude_from_maps              = ' // quote(ltt%exclude_from_maps)
  write (1, '(2a)') '# ltt%symplectic_map_tracking        = ' // logic_str(ltt%symplectic_map_tracking)
  write (1, '(2a)') '# Number_of_maps                     = ' // int_str(ltt_com%num_maps)
endif
write (1, '(2a)') '# ltt%split_bends_for_stochastic_rad = ' // logic_str(ltt%split_bends_for_stochastic_rad)
write (1, '(2a)') '# ltt%use_rf_clock                   = ' // logic_str(ltt%use_rf_clock)
write (1, '(2a)') '# ltt%ramping_on                     = ' // logic_str(ltt%ramping_on)
write (1, '(2a)') '# ltt%ramping_start_time             = ' // real_str(ltt%ramping_start_time, 6)
write (1, '(2a)') '# ltt%set_beambeam_z_crossing        = ' // logic_str(ltt%set_beambeam_z_crossing)
write (1, '(2a)') '# ltt%random_seed                    = ' // int_str(ltt%random_seed)
if (ltt%random_seed == 0) then
  write (1, '(2a)') '# random_seed_actual                 = ' // int_str(ltt_com%random_seed_actual)
endif
write (1, '(2a)') '# ltt%rfcavity_on                    = ' // logic_str(ltt%rfcavity_on)
write (1, '(2a)') '# is_RF_on                           = ' // logic_str(rf_is_on(branch)) // '  #  M65 /= 0 ?'
write (1, '(a, f10.6)')  '# da_param%min_angle      =', da_param%min_angle
write (1, '(a, f10.6)')  '# da_param%max_angle      =', da_param%max_angle
write (1, '(a, es10.2)') '# da_param%rel_accuracy   =', da_param%rel_accuracy
write (1, '(a, es10.2)') '# da_param%abs_accuracy   =', da_param%abs_accuracy
write (1, '(a, f10.5)')  '# da_param%x_init         =', da_param%x_init
write (1, '(a, f10.5)')  '# da_param%y_init         =', da_param%y_init
write (1, '(a, i0)')     '# da_param%n_turn         = ', da_param%n_turn
write (1, '(a, i0)')     '# da_param%n_angle        = ', da_param%n_angle
write (1, '(2a)')        '# bmad_com%radiation_damping_on      = ' // logic_str(bmad_com%radiation_damping_on)
write (1, '(2a)')        '# bmad_com%radiation_fluctuations_on = ' // logic_str(bmad_com%radiation_fluctuations_on)
write (1, '(2a)')        '## gnuplot plotting command:'
write (1, '(2a)')        '##   ', trim(gnu_command)

call dynamic_aperture_scan (aperture_scan, da_param, dpz(1:n_dpz), ltt_com%tracking_lat)

do i = 1, n_dpz
  da => aperture_scan(i)

  write (1, *)
  write (1, *)
  write (1, '(a, f10.6, a)') '"dpz =', dpz(i), '"'
  write (1, '(a, f10.6, a)') '"x_ref_orb =', da%ref_orb%vec(1), '"'
  write (1, '(a, f10.6, a)') '"y_ref_orb =', da%ref_orb%vec(3), '"'
  do j = 1, da_param%n_angle
    da_point => da%point(j)
    write (1, '(2f11.6, i7, 6x, a, 3x, a)') da_point%x, da_point%y, da_point%i_turn, &
                              coord_state_name(da_point%plane), trim(branch%ele(da_point%ix_ele)%name)
  enddo
enddo

close (1)

print '(2a)', 'gnuplot plotting command: ', trim(gnu_command)

end program
