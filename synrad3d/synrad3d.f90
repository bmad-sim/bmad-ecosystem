!+
! Program synrad3d
!
! Program to calculate photoelectron distributions in a lattice.
! See the synrad3d manual for more details.
!-

program synrad3d

use synrad3d_plot_mod
use synrad3d_output_mod
use synrad3d_test_mod
use synrad3d_parse_wall
use bookkeeper_mod

implicit none

type (ele_struct) ele_here
type (ele_struct), pointer :: ele
type (lat_struct), target :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) orbit_here, orbit, orb1
type (rad_int_all_ele_struct) rad_int_ele
type (normal_modes_struct) modes
type (branch_struct), pointer :: branch
type (wall3d_struct), pointer :: wall3d
type (sr3d_photon_track_struct), allocatable, target :: photons(:)
type (sr3d_photon_track_struct), pointer :: photon
type (sr3d_plot_param_struct) plot_param
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)
type (random_state_struct) ran_state

real(rp) ds_step_min, d_i0, i0_tot, ds, gx, gy, s_offset
real(rp) emit_a, emit_b, sig_e, g, gamma, r, dtrack, photon_number_factor
real(rp) e_filter_min, e_filter_max, s_filter_min, s_filter_max
real(rp) e_init_filter_min, e_init_filter_max, timer_time, old_time
real(rp) surface_roughness_rms, roughness_correlation_len, rms_set, correlation_set

integer i, n, iu, iu2, ix, random_seed, iu_start, j_photon, ix_ele, status, ix_branch
integer n_photon_generated, n_photon_array, i0_ele, n_photon_ele, n_photon_here
integer ix_ele_track_start, ix_ele_track_end, iu_hit_file, iu_lat_file
integer photon_direction, num_photons, num_photons_per_pass, n_phot, ios
integer n_photons_per_pass, num_ignore_generated_outside_wall, ix_photon_out

character(200) lattice_file, wall_hit_file, reflect_file, lat_ele_file, photon_track_file
character(200) photon_start_input_file, photon_start_output_file, surface_reflection_file
character(300) line3
character(100) dat_file, dat2_file, wall_file, param_file, arg, line
character(40) plotting, test, who
character(16) chamber_end_geometry
character(16) :: r_name = 'synrad3d'

logical ok, s_wrap_on, filter_this, err, filter_phantom_photons
logical is_inside, turn_off_kickers_in_lattice

namelist / synrad3d_parameters / ix_ele_track_start, ix_ele_track_end, chamber_end_geometry, &
            photon_direction, num_photons, lattice_file, ds_step_min, num_photons_per_pass, &
            emit_a, emit_b, sig_e, sr3d_params, wall_file, dat_file, random_seed, &
            e_filter_min, e_filter_max, s_filter_min, s_filter_max, wall_hit_file, &
            photon_start_input_file, photon_start_output_file, reflect_file, lat_ele_file, &
            num_ignore_generated_outside_wall, turn_off_kickers_in_lattice, &
            e_init_filter_min, e_init_filter_max, plot_param, surface_reflection_file, &
            surface_roughness_rms, roughness_correlation_len, photon_track_file, filter_phantom_photons

namelist / start / orbit, ix_branch, ran_state, random_seed

! Parse command line args

ok = .true.
plotting = ''
test = ''
param_file = ''
surface_reflection_file = ''
photon_start_input_file = '' 
ix_photon_out = -1

i = 0
do while (i < cesr_iargc())
  i = i + 1
  call cesr_getarg(i, arg)
  select case (arg)
  case ('-plot')
    i = i + 1
    call cesr_getarg(i, plotting)
    if (plotting == '') ok = .false.
  case ('-test')
    i = i + 1
    call cesr_getarg(i, test)
  case ('-in')
    i = i + 1
    call cesr_getarg(i, photon_start_input_file)
  case ('-out')
    i = i + 1
    call cesr_getarg(i, arg)
    read (arg, *) ix_photon_out
  case default
    if (arg(1:1) == '-') then
      print *, 'I DO NOT UNDERSTAND: ', trim(arg)
      print *
      ok = .false.
    endif
    param_file = arg
  end select
enddo

if (param_file == '') param_file = 'synrad3d.init'

if (.not. ok) then
  print '(a)', 'Usage:'
  print '(a)', '  synrad3d {options} {<init_file>}'
  print '(a)', 'Default:'
  print '(a)', '  <init_file> = synrad3d.init'
  print '(a)', 'Options: [Note: Standard photon tracking not done with -plot nor -test]'
  print '(a)', '  -plot <type>  ! <type> = "reflect", "xy", "xs", or "ys".'
  print '(a)', '  -test <who>   ! <who> = "diffuse_reflection" .or "specular_reflection".'
  print '(a)', '  -in <file>    ! Use <file> to initialize photon(s). '
  print '(a)', '                !   This option is used for debugging synrad3d.'
  print '(a)', '  -out <n>      ! Create error_photon_start file using the n^th generated photon.'
  print '(a)', '                !   This option is used for debugging synrad3d.'
  stop
endif

! test 

if (test /= '') then
  call match_word (trim(test), ['diffuse_reflection ', 'specular_reflection'], i, .true., .true., who)
  select case (who)
  case ('diffuse_reflection ')
    call sr3d_diffuse_reflection_test (param_file)
  case ('specular_reflection')
    call sr3d_specular_reflection_test (param_file)
  case default
    print '(a)', 'I DO NOT UNDERSTAND THIS TEST: ' // trim(test)
  end select
  stop
endif

! Get parameters.
! Radiation is produced from the end of ix_ele_track_start to the end of ix_ele_track_end.

random_seed = 0
ix_ele_track_start = 0    ! defaults
ix_ele_track_end = -1
ds_step_min = 0.01
emit_a = -1
emit_b = -1
sig_e  = -1
dat_file = 'synrad3d.dat'
wall_file = 'synrad3d.wall'
photon_track_file = ''
photon_direction = 1
e_init_filter_min = -1
e_init_filter_max = -1
e_filter_min = -1
e_filter_max = -1
s_filter_min = -1
s_filter_max = -1
filter_phantom_photons = .true.
wall_hit_file = ''
reflect_file = ''
lat_ele_file = ''
iu_lat_file = 0
photon_start_output_file = ''
num_photons = -1
num_photons_per_pass = -1
num_ignore_generated_outside_wall = 0
turn_off_kickers_in_lattice = .false.
surface_roughness_rms = -1; roughness_correlation_len = -1
chamber_end_geometry = ''

sr3d_params%debug_on = .false.
sr3d_params%ix_generated_warn = -1
sr3d_params%iu_wall_hit = 0

print *, 'Input parameter file: ', trim(param_file)
open (1, file = param_file, status = 'old')
read (1, nml = synrad3d_parameters)
close (1)

print *,'Lattice file: ',trim(lattice_file)
print *,'Wall file: ',trim(wall_file)

if (reflect_file /= '') wall_hit_file = reflect_file  ! Accept old syntax.
sr3d_params%photon_track_file = photon_track_file
sr3d_params%wall_hit_file = wall_hit_file

call ran_seed_put (random_seed)

! When a filter parameter is set, only photons that satisfy the filter criteria are kept

s_wrap_on = (s_filter_min >= 0) .and. (s_filter_max >= 0) .and. (s_filter_min > s_filter_max)

! Get lattice

call bmad_and_xsif_parser(lattice_file, lat)
branch => lat%branch(0)

select case (chamber_end_geometry)
case ('open')
  sr3d_params%chamber_end_geometry = open$
case ('closed')
  sr3d_params%chamber_end_geometry = closed$
case ('')
  sr3d_params%chamber_end_geometry = lat%branch(0)%param%geometry
case default
  print *, 'Bad "chamber_end_geometry" setting: ', chamber_end_geometry
  stop
end select

if (ix_ele_track_end < 0) ix_ele_track_end = branch%n_ele_track

! Wall init

call sr3d_read_wall_file (wall_file, lat)

! Load different surface reflection parameters if wanted

if (surface_reflection_file /= '') call read_surface_reflection_file (surface_reflection_file, lat%surface(1))
if (surface_roughness_rms > 0) lat%surface(1)%surface_roughness_rms = surface_roughness_rms
if (roughness_correlation_len > 0) lat%surface(1)%roughness_correlation_len = roughness_correlation_len

! Plot wall cross-sections or reflections. 
! The plotting routines never return back to the main program.

if (plotting /= '') then
  if (plotting == 'xy') then
    call sr3d_plot_wall_cross_sections (plot_param, lat)
  elseif (plotting == 'xs' .or. plotting == 'ys') then
    call sr3d_plot_wall_vs_s (plot_param, lat, plotting)
  elseif (index('reflect', trim(plotting)) == 1) then
    call sr3d_plot_reflection_probability(plot_param, lat)
  else
    call out_io (s_fatal$, r_name, 'I DO NOT UNDERSTAND WHAT TO PLOT: ' // plotting)
    stop
  endif
endif

! Lattice setup

if (turn_off_kickers_in_lattice) then
  do i = 1, branch%n_ele_max
    ele => branch%ele(i)
    if (attribute_index(ele, 'HKICK') > 0) then
      ele%value(hkick$) = 0
      ele%value(vkick$) = 0
      ele%value(bl_hkick$) = 0
      ele%value(bl_vkick$) = 0
    endif
    if (attribute_index(ele, 'KICK') > 0) then  
      ele%value(kick$) = 0
      ele%value(bl_kick$) = 0
    endif
    if (associated(ele%a_pole)) then
      ele%a_pole(0) = 0
      ele%b_pole(0) = 0
    endif
  enddo
  call lattice_bookkeeper (lat)
endif

call twiss_and_track (lat, orb, status, branch%ix_branch, use_beam_start = .true.)
if (status /= ok$) stop
  
! Find out much radiation is produced

call radiation_integrals (lat, orb, modes, rad_int_by_ele = rad_int_ele, ix_branch = branch%ix_branch)

if (ix_ele_track_end > ix_ele_track_start) then
  i0_tot = sum(rad_int_ele%ele(ix_ele_track_start+1:ix_ele_track_end)%i0)
else
  i0_tot = sum(rad_int_ele%ele(ix_ele_track_start+1:branch%n_ele_track)%i0) + &
           sum(rad_int_ele%ele(1:ix_ele_track_end)%i0)
endif

! To generate photons we either need bends or wigglers or a photon init file.

if (i0_tot == 0 .and. photon_start_input_file == '') then
  call out_io (s_fatal$, r_name, 'NO BENDS OR OTHER ELEMENTS TO GENERATE RADIATION IN REGION OF INTEREST!')
  call err_exit
endif

! Print some info

print *, 'I0 Radiation Integral of entire lattice:', modes%synch_int(0)
print *, 'I0 Radiation Integral over emission region:', i0_tot
print *, 'Closed orbit max X amplitude (meters):', maxval(abs(orb(:)%vec(1)))
print *, 'Closed orbit max Y amplitude (meters):', maxval(abs(orb(:)%vec(3)))

! d_i0 determines the number of photons to generatie per unit i0 integral.

n_photons_per_pass = num_photons_per_pass
if (n_photons_per_pass < 1) n_photons_per_pass = ceiling(0.2 * num_photons)
d_i0 = i0_tot / n_photons_per_pass

! Determine the emittance

if (emit_a < 0) then
  emit_a = modes%a%emittance
  print *, 'Using emit_a =', emit_a
endif

if (emit_b < 0) then
  emit_b = modes%b%emittance
  print *, 'Using emit_b =', emit_b
endif

if (sig_e < 0) then
  sig_e  = modes%sige_e
  print *, 'Using sig_e =', sig_e
endif

! Open data files.

if (sr3d_params%wall_hit_file /= '') then
  sr3d_params%iu_wall_hit = lunget()
  open (sr3d_params%iu_wall_hit, file = sr3d_params%wall_hit_file, recl = 300)
  print *, 'Creating photon hit point output file: ', trim(sr3d_params%wall_hit_file)
endif

if (lat_ele_file /= '') then
  iu_lat_file = lunget()
  open (iu_lat_file, file = lat_ele_file, recl = 120)
  print *, 'Creating lattice element output file: ', trim(lat_ele_file)
    write (iu_lat_file, *) 'I0 Radiation Integral of entire lattice:', modes%synch_int(0)
    write (iu_lat_file, *) 'I0 Radiation Integral over emission region:', i0_tot
    write (iu_lat_file, *) ''
    write (iu_lat_file, *) &
        'Index  Name                Type                  S       L          I0    N_phot  ds_step'
endif


sr3d_params%iu_dat_file = lunget()
open (sr3d_params%iu_dat_file, recl = 300, file = dat_file)
print *, 'Data file is: ', trim(dat_file)

if (sr3d_params%photon_track_file /= '') then
  open (5, file = sr3d_params%photon_track_file)
  sr3d_params%iu_photon_track = 5
endif

! Write header info

call write_this_header (sr3d_params%iu_dat_file)

! Track through the elements and generate photons.

bmad_com%auto_bookkeeper = .false.  ! Since we are not changing any element params.

n_photon_generated = 0
n_photon_array = 0

allocate (wall_hit(0:10))
print *, 'Initialization done. Tracking beginning...'
call run_timer ('START')
old_time = -1d3

!--------------------------------------------------------------------------
! If the photon_start input file exists then use that

if (photon_start_input_file /= '') then

  ! Open photon start input file and count the number of photons

  print *, 'Opening photon starting position input file: ', trim(photon_start_input_file)
  open (1, file = photon_start_input_file, status = 'old')

  allocate (photons(1000))

  ! Now read the file and initialize photons and the random number generator.
  ! Only set the random number generator if ran_state is set in the file.

  photon_loop: do

    n_photon_array = n_photon_array + 1
    if (n_photon_array > size(photons)) call reallocate_photon_array (photons, 2*size(photons))
    photon => photons(n_photon_array)

    do
      ran_state%iy = -1  ! To see if ran_state is set by the read.
      random_seed = -1
      ix_branch = 0   ! Default
      read (1, nml = start, iostat = ios)
      if (ios < 0) exit photon_loop
      if (ios > 0) then
        print *, 'Error reading photon starting position at photon index:', n_photon_generated
        rewind (1)
        do
          read (1, nml = start) ! will generate error message
        enddo
      endif
      ix_ele = element_at_s(lat, orbit%s, .true., branch%ix_branch)
      call init_coord (photon%start%orb, orbit%vec, branch%ele(ix_ele), inside$, 0, photon$, orbit%p0c)
      photon%start%ix_branch = ix_branch
      photon%start%orb%s = orbit%s

      call sr3d_check_if_photon_init_coords_outside_wall (photon%start, lat, is_inside, num_ignore_generated_outside_wall)
      if (is_inside) exit
    enddo

    n_photon_generated = n_photon_generated + 1
    photon%ix_photon_generated = n_photon_generated
    photon%n_wall_hit = 0

    if (ran_state%iy > 0) call ran_default_state (set_state = ran_state)
    if (random_seed > -1) call ran_seed_put (seed = random_seed)
    call check_filter_restrictions(ok, .true.)
    if (.not. ok) cycle

    if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('START_RECORDING')
    call sr3d_track_photon (photon, lat, wall_hit, err)

    ! ix_photon_out is used for generating a file of the photon starting position.
    ! This is used for diagnostic purposes.

    if (photon%ix_photon_generated == ix_photon_out) then
      call sr3d_print_hit_points (-1, photon, wall_hit, branch)
      stop
     endif

    !

    if (err) then
      n_photon_array = n_photon_array - 1  ! Delete photon from the array.
      cycle
    endif

    call check_filter_restrictions(ok, .false.)
    if (ok) then
      if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('MOVE_TRACK_TO_FILE')
    else
      if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('ERASE_RECORDING')
      cycle
    endif

    call sr3d_print_hit_points (sr3d_params%iu_wall_hit, photon, wall_hit, branch)
    call write_photon_data (n_photon_array, photon)

  enddo photon_loop

  close (1)

!--------------------------------------------------------------------------
! Regular photon generation

else

  ! Open photon start output file

  if (photon_start_output_file /= '') then
    iu_start = lunget()
    open (iu_start, file = photon_start_output_file, recl = 140)
    print *, 'Creating photon start output file: ', trim(photon_start_output_file)
  endif

  ! Track loop

  allocate (photons(num_photons + n_photons_per_pass + 1))

  ix_ele = ix_ele_track_start
  do 

    if (ix_ele == ix_ele_track_end) then
      if (n_photon_array >= num_photons) exit
      ix_ele = ix_ele_track_start
      if (iu_lat_file > 0) close (iu_lat_file)
      iu_lat_file = 0 ! To stop further output
      if (n_photon_generated == 0) then
        print *, 'NO PHOTONS GENERATED. N_PHOTONS OR N_PHOTONS_PER_PASS IS TOO SMALL!'
        call err_exit
      endif
    endif

    ix_ele = ix_ele + 1
    if (ix_ele == branch%n_ele_track+1) ix_ele = 0

    ele => branch%ele(ix_ele)

    n_phot = nint(rad_int_ele%ele(ix_ele)%i0 / d_i0)
    if (n_phot == 0) cycle

    ds = ele%value(l$) / n_phot
    if (ds < ds_step_min) ds = (1+int(ds_step_min/ds)) * ds

    ! Write info to lat_ele_file

    if (iu_lat_file > 0) then
      write (iu_lat_file, '(i6, 2x, a20, a16, f10.3, f8.3, f10.1, i9, f8.3)') ix_ele, ele%name, &
              key_name(ele%key), ele%s, ele%value(l$), rad_int_ele%ele(ix_ele)%i0, n_phot, ds
    endif

    ! Loop over all photon generating points.
    ! First point is random in the range [0, ds] to avoid correlations between passes when
    ! there are multiple passes.

    call ran_uniform (r)
    s_offset = r * ds 
    i0_ele = 0         ! Integrated i0 for this element
    n_photon_ele = 0   

    do

      call run_timer ('READ', timer_time)
      if (old_time > old_time + 60) then
        print '(a, f9.1)', 'Time from start (min):', timer_time/60
        print *, '    Num photons generated:          ', n_photon_generated
        print *, '    Num photons passed filter tests:', n_photon_array
        old_time = timer_time
      endif

      call sr3d_get_emission_pt_params (branch, orb, ix_ele, s_offset, ele_here, orbit_here, gx, gy)
      g = sqrt(gx**2 + gy**2) 
      call convert_total_energy_to (ele%value(e_tot$),  branch%param%particle, gamma)
      ! Generate photons, track to the wall 

      n_photon_here = nint(g * gamma * ds / d_i0)
      do j_photon = 1, n_photon_here
        n_photon_generated = n_photon_generated + 1
        n_photon_array = n_photon_array + 1
        if (n_photon_array > size(photons)) then
          print *, 'INTERNAL ERROR: NUMBER OF PHOTONS GENERATED TOO LARGE!'
          call err_exit
        endif
        photon => photons(n_photon_array)
        photon%ix_photon = n_photon_array
        photon%ix_photon_generated = n_photon_generated

        if (n_photon_generated == sr3d_params%ix_generated_warn) then
          print *, 'Note: At sr3d_params%ix_generated_warn:', n_photon_generated ! For debug purposes.
        endif

        photon%n_wall_hit = 0

        do
          call sr3d_emit_photon (ele_here, orbit_here, gx, gy, emit_a, emit_b, sig_e, photon_direction, photon%start)
          call sr3d_check_if_photon_init_coords_outside_wall (photon%start, lat, is_inside, num_ignore_generated_outside_wall)
          if (is_inside) exit
        enddo

        if (photon_start_output_file /= '') then
          call ran_default_state (get_state = ran_state)
          write (iu_start, '(a)')           '&start'
          write (iu_start, '(a, 6es20.12)') '  p%vec     =', photon%start%orb%vec
          write (iu_start, '(a, es20.12)')  '  p%s       =', photon%start%orb%s
          write (iu_start, '(a, es20.12)')  '  p%p0c     =', photon%start%orb%p0c
          write (iu_start, '(a, i4)')       '  ix_branch =', ix_branch
          write (iu_start, *)               '  ran_state = ', ran_state
          write (iu_start, '(a)')           '/'
        endif

        call check_filter_restrictions(ok, .true.)
        if (.not. ok) cycle

        if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('START_RECORDING')
        call sr3d_track_photon (photon, lat, wall_hit, err)

        ! ix_photon_out is used for generating a file of the photon starting position.
        ! This is used for diagnostic purposes.

        if (photon%ix_photon_generated == ix_photon_out) then
          call sr3d_print_hit_points (-1, photon, wall_hit, branch)
          stop
        endif

        !

        if (err) then
          n_photon_array = n_photon_array - 1  ! Delete photon from the array.
          cycle
        endif

        call check_filter_restrictions (ok, .false.)
        if (ok) then
          if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('MOVE_TRACK_TO_FILE')
        else
          if (sr3d_params%iu_photon_track > 0) call sr3d_record_photon_position('ERASE_RECORDING')
          cycle
        endif

        call sr3d_print_hit_points (sr3d_params%iu_wall_hit, photon, wall_hit, branch)
        call write_photon_data (n_photon_array, photon)

      enddo

      s_offset = s_offset + ds
      if (s_offset > ele%value(l$)) exit

    enddo

  enddo

endif

! Write photon_number_factor = (Num actual photons emitted per beam particle) / (Num macro photons generated in simulation)

iu = sr3d_params%iu_dat_file
iu2 = lunget()
open (iu2, status = 'scratch')

rewind(iu)

do
  read (iu, '(a)', iostat = ios) line3
  write (iu2, '(a)') trim(line3)
  if (ios /= 0) exit
enddo

rewind(iu)
rewind(iu2)

photon_number_factor = 5 * sqrt(3.0) * classical_radius_factor * i0_tot / (6 * h_bar_planck * c_light * n_photon_generated)

write (iu, '(a, es11.3, a)') '# photon_number_factor       =', photon_number_factor
write (iu, '(a, i11, a)')    '# num_photons_generated      =', n_photon_generated, '   ! Total including filtered photons.'
write (iu, '(a, i11, a)')    '# num_photons_passed_test    =', n_photon_array,     '   ! Num that passed filter tests'
write (iu, '(a, es11.3, a)') '# I0_tot                     =', i0_tot, '   ! I0 radiation integral for the entire ring'

do
  read (iu2, '(a)', iostat = ios) line3
  write (iu, '(a)') trim(line3)
  if (ios /= 0) exit
enddo

close (iu)
close (iu2)

!--------------------------------------------------------------------------------------------
contains

subroutine write_photon_data (n_photon, photon)

type (sr3d_photon_track_struct) :: photon
type (branch_struct), pointer :: branch
real(rp) start_vec(6), now_vec(6), dtrack
integer n_photon, j, iu
character(40) wall_name

!

branch => lat%branch(photon%now%ix_branch)
start_vec = [photon%start%orb%vec(1:4), photon%start%orb%s, photon%start%orb%vec(6)]
now_vec = [photon%now%orb%vec(1:4), photon%now%orb%s, photon%now%orb%vec(6)]
dtrack = photon%now%orb%path_len - photon_direction * modulo2((photon%now%orb%s - photon%start%orb%s), branch%param%total_length/2)
j = photon%now%orb%ix_ele
wall_name = branch%wall3d(photon%now%ix_wall3d)%name
if (wall_name == '') wall_name = '<default_subchamber>'

iu = sr3d_params%iu_dat_file
write (iu, '(i8, i6, es14.6, (4f12.6, f12.3, f12.6), i4, (4f12.6, f12.3, f12.6), 2f12.6, i4, i7, 3x, a16, 3x, a)') &
        n_photon, photon%n_wall_hit, photon%start%orb%p0c, start_vec, photon%start%ix_branch, now_vec, &
        photon%now%orb%path_len, dtrack, branch%ix_branch, j, key_name(branch%ele(j)%key), trim(wall_name)

end subroutine write_photon_data

!------------------------------------------------------------------------------------------
! contains
!+
! Subroutine check_filter_restrictions (ok, check_init_filters)
!
! Routine to check if a photon has passed the filter requirements.
!
! Input:
!   check_init_filters -- logical: If true then photon coordinates represent the photon
!                           at generation time. This involves different filter cuts.
!   
! Output:
!   ok -- logical: Set True if passed. False otherwise.
!-

subroutine check_filter_restrictions (ok, check_init_filters)

type (branch_struct), pointer :: branch
type (sr3d_coord_struct), pointer :: now
type (photon_reflect_surface_struct), pointer :: surface
logical ok, check_init_filters

! Check filter restrictions

ok = .true.

filter_this = .false.

if (check_init_filters) then
  now => photon%start   ! Photon is just starting so use %start and not %now
  if (e_init_filter_min > 0 .and. now%orb%p0c < e_init_filter_min) filter_this = .true.
  if (e_init_filter_max > 0 .and. now%orb%p0c > e_init_filter_max) filter_this = .true.

else
  now => photon%now
  if (e_filter_min > 0 .and. now%orb%p0c < e_filter_min) filter_this = .true.
  if (e_filter_max > 0 .and. now%orb%p0c > e_filter_max) filter_this = .true.
  if (s_wrap_on) then
    if (now%orb%s > s_filter_max .and. now%orb%s < s_filter_min) filter_this = .true.
  else
    if (s_filter_min > 0 .and. now%orb%s < s_filter_min) filter_this = .true.
    if (s_filter_max > 0 .and. now%orb%s > s_filter_max) filter_this = .true.
  endif

  if (filter_phantom_photons) then
    branch => lat%branch(now%ix_branch)
    surface => branch%wall3d(now%ix_wall3d)%section(now%ix_wall_section+1)%surface
    if (surface%name == 'PHANTOM') filter_this = .true.
  endif
endif

!

if (filter_this) then
  n_photon_array = n_photon_array - 1  ! Delete photon from the array.
  ok = .false.
endif

end subroutine

!--------------------------------------------------------------------------------------------
! contains

!+
! Subroutine reallocate_photon_array (photon_array, n_size)
! 
! Routine to enlarge an array of photons while saving the original information.
! The final size will be at least n_size.
!
! Input:
!   photon_array(:) -- sr3d_photon_track_struct, allocatable: Array of photons
!   n_size          -- Integer: Minimum size
!
! Output:
!   photon_array(:) -- sr3d_photon_track_struct, allocatable: Array with size at least n_size.
!-

subroutine reallocate_photon_array (photon_array, n_size)

implicit none

type (sr3d_photon_track_struct), allocatable :: photon_array(:), temp(:)
integer n, n_size, n_old, j

!

if (.not. allocated(photon_array)) then
  allocate (photon_array(n_size))
  return
endif

n_old = size(photon_array)
if (n_old >= n_size) return

allocate(temp(n_old))
temp = photon_array

deallocate(photon_array)
allocate(photon_array(n_size))

photon_array(1:n_old) = temp

deallocate (temp)

end subroutine

!--------------------------------------------------------------------------------------------
! contains

!+
! Subroutine write_this_header (iu)
!
! Routine to write the header into to the output data file.
!-

subroutine write_this_header (iu)

integer iu
character(40) line
character(20) date_and_time

call date_and_time_stamp (date_and_time)
line = ''
write (iu, '(2a)')           '# date                       = ', date_and_time
write (iu, '(a, i0, a)')     '# ix_ele_track_start         = ', ix_ele_track_start
write (iu, '(a, i0, a)')     '# ix_ele_track_end           = ', ix_ele_track_end
write (iu, '(a, i0, a)')     '# photon_direction           = ', photon_direction
write (iu, '(a, i0, a)')     '# num_photons                = ', num_photons, '   ! Input target number to generate'
write (iu, '(a, i0, a)')     '# num_photons_per_pass       = ', num_photons_per_pass
write (iu, '(a, i0, a)')     '# random_seed                = ', random_seed
write (iu, '(a, 3a)')        '# lattice_file               = ', '"', trim(lattice_file), '"'
write (iu, '(a, 3a)')        '# photon_start_input_file    = ', '"', trim(photon_start_input_file), '"'
write (iu, '(a, 3a)')        '# wall_file                  = ', '"', trim(wall_file), '"'
write (iu, '(a, 3a)')        '# dat_file                   = ', '"', trim(dat_file), '"'
write (iu, '(a, 3a)')        '# chamber_end_geometry       = ', '"', trim(chamber_end_geometry), '"'
write (iu, '(a, es10.3)')    '# ds_step_min                = ', ds_step_min
write (iu, '(a, es10.3)')    '# emit_a                     = ', emit_a
write (iu, '(a, es10.3)')    '# emit_b                     = ', emit_b
write (iu, '(a, es10.3)')    '# sig_e                      = ', sig_e
write (iu, '(a, es10.3)')    '# e_init_filter_min          = ', e_init_filter_min
write (iu, '(a, es10.3)')    '# e_init_filter_max          = ', e_init_filter_max
write (iu, '(a, es10.3)')    '# e_filter_min               = ', e_filter_min
write (iu, '(a, es10.3)')    '# e_filter_max               = ', e_filter_max
write (iu, '(a, f11.4)')     '# s_filter_min               = ', s_filter_min
write (iu, '(a, f11.4)')     '# s_filter_max               = ', s_filter_max
write (iu, '(a, l1)')        '# filter_phantom_photons     = ', filter_phantom_photons
write (iu, '(a, es10.3)')    '# sr3d_params%ds_track_step_max         = ', sr3d_params%ds_track_step_max
write (iu, '(a, es10.3)')    '# sr3d_params%dr_track_step_max         = ', sr3d_params%dr_track_step_max
write (iu, '(a, es10.3)')    '# surface_roughness_rms (input)         = ', surface_roughness_rms
write (iu, '(a, es10.3)')    '# surface_roughness_rms (set value)     = ', lat%surface(1)%surface_roughness_rms
write (iu, '(a, es10.3)')    '# roughness_correlation_len (input)     = ', roughness_correlation_len
write (iu, '(a, es10.3)')    '# roughness_correlation_len (set value) = ', lat%surface(1)%roughness_correlation_len
write (iu, '(a, l1)')        '# sr3d_params%allow_reflections         = ', sr3d_params%allow_reflections
write (iu, '(a, l1)')        '# sr3d_params%specular_reflection_only  = ', sr3d_params%specular_reflection_only
write (iu, '(a)') '#                                                                                                                                                                                                         Ix_Branch'
write (iu, '(a)') '# Photon   Num               |                        Init Postion                                   Ix_   |                        Final Position                                 |   Path     ds_photon    Ix_Ele   Ele Type'
write (iu, '(a)') '#  index   Hit   E_photon    |    x           Vx          y           Vy          s           Vz    Branch |    x           Vx          y          Vy          s           Vz      |   Length   - ds_beam     @ Hit   @ Hit              Chamber_name'

end subroutine

end program
