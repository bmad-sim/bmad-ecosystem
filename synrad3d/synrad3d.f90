!+
! Program synrad3d
!
! Program to calculate photoelectron distributions in a lattice
!-

program synrad3d

use synrad3d_plot_mod
use synrad3d_output_mod
use synrad3d_test_mod
use bookkeeper_mod

implicit none

type (ele_struct) ele_here
type (ele_struct), pointer :: ele
type (lat_struct), target :: lat
type (coord_struct), allocatable :: orb(:)
type (coord_struct) orbit_here
type (rad_int_all_ele_struct) rad_int_ele
type (normal_modes_struct) modes
type (sr3d_photon_track_struct), allocatable, target :: photons(:)
type (sr3d_photon_track_struct), pointer :: photon
type (sr3d_wall_struct), target :: wall
type (sr3d_photon_coord_struct) p
type (sr3d_plot_param_struct) plot_param
type (random_state_struct) ran_state
type (sr3d_photon_wall_hit_struct), allocatable :: wall_hit(:)

real(rp) ds_step_min, d_i0, i0_tot, ds, gx, gy, s_offset
real(rp) emit_a, emit_b, sig_e, g, gamma, r, dtrack, photon_number_factor
real(rp) e_filter_min, e_filter_max, s_filter_min, s_filter_max
real(rp) e_init_filter_min, e_init_filter_max, timer_time
real(rp) surface_roughness_rms, roughness_correlation_len, rms_set, correlation_set

integer i, j, n, iu, ix, random_seed, iu_start
integer ix_ele, n_photon_generated, n_photon_array, i0_ele, n_photon_ele, n_photon_here
integer ix_ele_track_start, ix_ele_track_end, iu_hit_file, iu_lat_file
integer photon_direction, num_photons, num_photons_per_pass, n_phot, ios
integer n_photons_per_pass, num_ignore_generated_outside_wall, ix_photon_out

character(200) lattice_file, wall_hit_file, reflect_file, lat_ele_file
character(200) photon_start_input_file, photon_start_output_file, surface_reflection_file

character(100) dat_file, dat2_file, wall_file, param_file, arg, line
character(40) plotting, test, who
character(16) :: r_name = 'synrad3d'

logical ok, filter_on, s_wrap_on, filter_this, err
logical is_inside, turn_off_kickers_in_lattice

namelist / synrad3d_parameters / ix_ele_track_start, ix_ele_track_end, &
            photon_direction, num_photons, lattice_file, ds_step_min, num_photons_per_pass, &
            emit_a, emit_b, sig_e, sr3d_params, wall_file, dat_file, random_seed, &
            e_filter_min, e_filter_max, s_filter_min, s_filter_max, wall_hit_file, &
            photon_start_input_file, photon_start_output_file, reflect_file, lat_ele_file, &
            num_ignore_generated_outside_wall, turn_off_kickers_in_lattice, &
            e_init_filter_min, e_init_filter_max, plot_param, surface_reflection_file, &
            surface_roughness_rms, roughness_correlation_len

namelist / start / p, ran_state, random_seed

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
photon_direction = 1
e_init_filter_min = -1
e_init_filter_max = -1
e_filter_min = -1
e_filter_max = -1
s_filter_min = -1
s_filter_max = -1
wall_hit_file = ''
reflect_file = ''
lat_ele_file = ''
iu_lat_file = 0
iu_hit_file = 0
photon_start_output_file = ''
num_photons = -1
num_photons_per_pass = -1
num_ignore_generated_outside_wall = 0
turn_off_kickers_in_lattice = .false.
surface_roughness_rms = -1; roughness_correlation_len = -1

sr3d_params%debug_on = .false.
sr3d_params%ix_generated_warn = -1

print *, 'Input parameter file: ', trim(param_file)
open (1, file = param_file, status = 'old')
read (1, nml = synrad3d_parameters)
close (1)

if (reflect_file /= '') wall_hit_file = reflect_file  ! Accept old syntax.

! When a filter parameter is set, only photons that satisfy the filter criteria are kept

filter_on = (e_init_filter_min > 0) .or. (e_init_filter_max > 0) .or. &
            (e_filter_min > 0) .or. (e_filter_max > 0) .or. (s_filter_min >= 0) .or. (s_filter_max >= 0)
s_wrap_on = (s_filter_min >= 0) .and. (s_filter_max >= 0) .and. (s_filter_min > s_filter_max)

! Get lattice

if (lattice_file(1:6) == 'xsif::') then
  call xsif_parser(lattice_file(7:), lat)
else
  call bmad_parser (lattice_file, lat)
endif

if (turn_off_kickers_in_lattice) then
  do i = 1, lat%n_ele_max
    ele => lat%ele(i)
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

call twiss_and_track (lat, orb, ok)
if (.not. ok) stop
  
if (ix_ele_track_end < 0) ix_ele_track_end = lat%n_ele_track

call ran_seed_put (random_seed)

! Wall init

call sr3d_read_wall_file (wall_file, lat%ele(lat%n_ele_track)%s, lat%param%geometry, wall)

! Load different surface reflection parameters if wanted

if (surface_reflection_file /= '') call read_surface_reflection_file (surface_reflection_file, wall%surface(1)%info)
if (surface_roughness_rms > 0) wall%surface(1)%info%surface_roughness_rms = surface_roughness_rms
if (roughness_correlation_len > 0) wall%surface(1)%info%roughness_correlation_len = roughness_correlation_len

! Plot wall cross-sections or reflections. 
! The plotting routines never return back to the main program.

if (plotting /= '') then
  if (plotting == 'xy') then
    call sr3d_plot_wall_cross_sections (plot_param, wall, lat)
  elseif (plotting == 'xs' .or. plotting == 'ys') then
    call sr3d_plot_wall_vs_s (plot_param, wall, lat, plotting)
  elseif (index('reflect', trim(plotting)) == 1) then
    call sr3d_plot_reflection_probability(plot_param, wall)
  else
    call out_io (s_fatal$, r_name, 'I DO NOT UNDERSTAND WHAT TO PLOT: ' // plotting)
    call err_exit
  endif
endif

! Find out much radiation is produced

call radiation_integrals (lat, orb, modes, rad_int_by_ele = rad_int_ele)

if (ix_ele_track_end > ix_ele_track_start) then
  i0_tot = sum(rad_int_ele%ele(ix_ele_track_start+1:ix_ele_track_end)%i0)
else
  i0_tot = sum(rad_int_ele%ele(ix_ele_track_start+1:lat%n_ele_track)%i0) + &
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
if (n_photons_per_pass < 1) n_photons_per_pass = num_photons
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

if (wall_hit_file /= '') then
  iu_hit_file = lunget()
  open (iu_hit_file, file = wall_hit_file)
  print *, 'Creating photon hit point output file: ', trim(wall_hit_file)
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

open (1, file = dat_file)

open (3, file = trim(dat_file) // '_table', recl = 240)
print *, 'Data file is: ', trim(dat_file)
print *, 'Data file in table format is: ', trim(dat_file) // '_table'

if (sr3d_params%stop_if_hit_antechamber) then
  dat2_file = trim(dat_file) // '.antechamber'
  open (2, file = dat2_file)
  print *, 'Data file for photons hitting the antechamber: ', trim(dat2_file)
endif

! Write header info

call write_this_header (1)

! Track through the elements and generate photons.

bmad_com%auto_bookkeeper = .false.  ! Since we are not changing any element params.

n_photon_generated = 0
n_photon_array = 0

allocate (wall_hit(0:10))
call run_timer ('START')

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

    do
      ran_state%iy = -1  ! To see if ran_state is set by the read.
      random_seed = -1
      read (1, nml = start, iostat = ios)
      if (ios < 0) exit photon_loop
      if (ios > 0) then
        print *, 'Error reading photon starting position at photon index:', n_photon_generated
        call err_exit
      endif
      call sr3d_check_if_photon_init_coords_outside_wall (p, wall, is_inside, num_ignore_generated_outside_wall)
      if (is_inside) exit
    enddo

    p%ix_ele = element_at_s(lat, p%s, .true.)

    n_photon_generated = n_photon_generated + 1
    n_photon_array = n_photon_array + 1
    if (n_photon_array > size(photons)) call reallocate_photon_array (photons, 2*size(photons))
    photon => photons(n_photon_array)
    photon%start = p
    photon%n_wall_hit = 0
    photon%ix_photon_generated = n_photon_generated
    if (ran_state%iy > 0) call ran_default_state (set_state = ran_state)
    if (random_seed > -1) call ran_seed_put (seed = random_seed)
    call check_filter_restrictions(ok, .true.)
    if (.not. ok) cycle
    call sr3d_track_photon (photon, lat, wall, wall_hit, err)

    ! ix_photon_out is used for generating a file of the photon starting position.
    ! This is used for diagnostic purposes.

    if (photon%ix_photon_generated == ix_photon_out) then
      call print_hit_points (-1, photon, wall_hit)
       stop
     endif

    !

    if (err) then
      n_photon_array = n_photon_array - 1  ! Delete photon from the array.
      cycle
    endif

    call check_filter_restrictions(ok, .false.)
    if (.not. ok) cycle

    call print_hit_points (iu_hit_file, photon, wall_hit)
    call write_photon_data (n_photon_array, photon)

  enddo photon_loop

  close (1)

! Regular photon generation

else

  ! Open photon start output file

  if (photon_start_output_file /= '') then
    iu_start = lunget()
    open (iu_start, file = photon_start_output_file, recl = 140)
    print *, 'Creating photon start output file: ', trim(photon_start_output_file)
  endif

  ! Allocate photons array

  n = max(num_photons, n_photons_per_pass)
  if (filter_on) then
    allocate (photons(nint(2.1*n)))   
  else
    allocate (photons(nint(1.1*n)))   ! Allow for some slop
  endif

  ix_ele = ix_ele_track_start
  do 

    if (ix_ele == ix_ele_track_end) then
      if (.not. filter_on .or. n_photon_array > 0.9 * num_photons) exit
      ix_ele = ix_ele_track_start
      if (iu_lat_file > 0) close (iu_lat_file)
      iu_lat_file = 0 ! To stop further output
      call run_timer ('READ', timer_time)
      print *, 'Time from start (min):', nint(timer_time/60)
      print *, '    Num photons generated:          ', n_photon_generated
      print *, '    Num photons passed filter tests:', n_photon_array
      if (n_photon_generated == 0) then
        print *, 'NO PHOTONS GENERATED. N_PHOTONS OR N_PHOTONS_PER_PASS IS TOO SMALL!'
        call err_exit
      endif
    endif

    ix_ele = ix_ele + 1
    if (ix_ele == lat%n_ele_track+1) ix_ele = 0

    ele => lat%ele(ix_ele)

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
      call sr3d_get_emission_pt_params (lat, orb, ix_ele, s_offset, ele_here, orbit_here, gx, gy)
      g = sqrt(gx**2 + gy**2) 
      call convert_total_energy_to (ele%value(e_tot$),  lat%param%particle, gamma)
      ! Generate photons, track to the wall 

      n_photon_here = nint(g * gamma * ds / d_i0)
      do j = 1, n_photon_here
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
        photon%start%ix_ele = ix_ele

        do
          call sr3d_emit_photon (ele_here, orbit_here, gx, gy, &
                               emit_a, emit_b, sig_e, photon_direction, photon%start)
          call sr3d_check_if_photon_init_coords_outside_wall (photon%start, wall, is_inside, num_ignore_generated_outside_wall)
          if (is_inside) exit
        enddo

        if (photon_start_output_file /= '') then
          call ran_default_state (get_state = ran_state)
          write (iu_start, '(a)')           '&start'
          write (iu_start, '(a, 6es20.12)') '  p%vec     = ', photon%start%vec
          write (iu_start, '(a, es20.12)')  '  p%s       =', photon%start%s
          write (iu_start, '(a, es20.12)')  '  p%energy  =', photon%start%energy
          write (iu_start, *)               '  ran_state = ', ran_state
          write (iu_start, '(a)')           '/'
        endif

        call check_filter_restrictions(ok, .true.)
        if (.not. ok) cycle

        call sr3d_track_photon (photon, lat, wall, wall_hit, err)

        ! ix_photon_out is used for generating a file of the photon starting position.
        ! This is used for diagnostic purposes.

        if (photon%ix_photon_generated == ix_photon_out) then
          call print_hit_points (-1, photon, wall_hit)
          stop
        endif

        !

        if (err) then
          n_photon_array = n_photon_array - 1  ! Delete photon from the array.
          cycle
        endif

        call check_filter_restrictions (ok, .false.)
        if (.not. ok) cycle

        call print_hit_points (iu_hit_file, photon, wall_hit)
        call write_photon_data (n_photon_array, photon)

      enddo

      s_offset = s_offset + ds
      if (s_offset > ele%value(l$)) exit

    enddo

  enddo

endif

! Write photon_number_factor = (Num actual photons emitted per beam particle) / (Num macro photons generated in simulation)
! Using direct access to be able to modify the first line in the file without touching the rest of the file.

close (1)

open (1, file = dat_file, access = 'direct', recl = 36, form = 'formatted')
photon_number_factor = 5 * sqrt(3.0) * classical_radius_factor * i0_tot / &
                                             (6 * h_bar_planck * c_light * n_photon_generated)

write (line, '(a, es11.3)') 'photon_number_factor    =', photon_number_factor
write (1, '(a36)', rec = 1) line

close (1)

if (sr3d_params%stop_if_hit_antechamber) close (2)

!--------------------------------------------------------------------------------------------
contains

subroutine write_photon_data (n_photon, photon)

type (sr3d_photon_track_struct) :: photon
real(rp) start_vec(6), now_vec(6)
integer n_photon

!

start_vec = [photon%start%vec(1:4), photon%start%s, photon%start%vec(6)]
now_vec = [photon%now%vec(1:4), photon%now%s, photon%now%vec(6)]

iu = 1
if (sr3d_params%stop_if_hit_antechamber .and. photon%hit_antechamber) iu = 2
write (iu, '(2i8, f12.4, 2x, a)') n_photon, photon%n_wall_hit, photon%start%energy, '! index, n_wall_hit, eV'
write (iu, '(4f12.6, f12.3, f12.6, a)') start_vec, '  ! Start position'
write (iu, '(4f12.6, f12.3, f12.6, a)') now_vec,   '  ! End position'
write (iu, '(f12.6, a)') photon%now%track_len, '  ! photon_track_len' 
dtrack = photon%now%track_len - photon_direction * &
    modulo2((photon%now%s - photon%start%s), lat%param%total_length/2)
write (iu, '(f12.6, a)') dtrack, '  ! photon_track_len - ds_beam'
j = photon%now%ix_ele
write (iu, '(i8, 3x, 2a)') j, key_name(lat%ele(j)%key), '  ! Lat ele index and class'

if (iu == 1) then
  write (3, '(2i8, es14.6, 2(4f12.6, f12.3, f12.6), 2f12.6, i8, 3x, a)') &
        n_photon, photon%n_wall_hit, photon%start%energy, start_vec, now_vec, &
        photon%now%track_len, dtrack, j, trim(key_name(lat%ele(j)%key)) 
endif

end subroutine write_photon_data

!------------------------------------------------------------------------------------------
! contains
!+
! Subroutine check_filter_restrictions (ok, init_filter)
!
! Routine to check if a photon has passed the filter requirements.
!
! Input:
!   init_filter -- logical: If true then photon coordinates represent the photon
!                    at generation time. This involves different filter cuts.
!   
! Output:
!   ok -- logical: Set True if passed. False otherwise.
!-

subroutine check_filter_restrictions (ok, init_filter)

logical ok, init_filter

! Check filter restrictions

ok = .true.

if (filter_on) then
  filter_this = .false.
  if (init_filter) then
    if (e_init_filter_min > 0 .and. photon%now%energy < e_init_filter_min) filter_this = .true.
    if (e_init_filter_max > 0 .and. photon%now%energy > e_init_filter_max) filter_this = .true.
  else
    if (e_filter_min > 0 .and. photon%now%energy < e_filter_min) filter_this = .true.
    if (e_filter_max > 0 .and. photon%now%energy > e_filter_max) filter_this = .true.
    if (s_wrap_on) then
      if (photon%now%s > s_filter_max .and. photon%now%s < s_filter_min) filter_this = .true.
    else
      if (s_filter_min > 0 .and. photon%now%s < s_filter_min) filter_this = .true.
      if (s_filter_max > 0 .and. photon%now%s > s_filter_max) filter_this = .true.
    endif
  endif

  if (filter_this) then
    n_photon_array = n_photon_array - 1  ! Delete photon from the array.
    ok = .false.
  endif

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

write (iu, '(a36)') 'photon_number_factor    = 0.000E+00  '
write (iu, '(a, i0)') 'ix_ele_track_start      = ', ix_ele_track_start
write (iu, '(a, i0)') 'ix_ele_track_end        = ', ix_ele_track_end
write (iu, '(a, i0)') 'photon_direction        = ', photon_direction
write (iu, '(a, i0)') 'num_photons             = ', num_photons
write (iu, '(a, i0)') 'num_photons_per_pass    = ', num_photons_per_pass
write (iu, '(a, i0)') 'random_seed             = ', random_seed
write (iu, '(a, a)') 'lattice_file            = ', trim(lattice_file)
write (iu, '(a, a)') 'photon_start_input_file = ', trim(photon_start_input_file)
write (iu, '(a, a)') 'wall_file               = ', trim(wall_file)
write (iu, '(a, a)') 'dat_file                = ', trim(dat_file)
write (iu, '(a, es10.3)') 'ds_step_min             = ', ds_step_min
write (iu, '(a, es10.3)') 'emit_a                  = ', emit_a
write (iu, '(a, es10.3)') 'emit_b                  = ', emit_b
write (iu, '(a, es10.3)') 'sig_e                   = ', sig_e
write (iu, '(a, es10.3)') 'e_init_filter_min       = ', e_init_filter_min
write (iu, '(a, es10.3)') 'e_init_filter_max       = ', e_init_filter_max
write (iu, '(a, es10.3)') 'e_filter_min            = ', e_filter_min
write (iu, '(a, es10.3)') 'e_filter_max            = ', e_filter_max
write (iu, '(a,  f11.4)') 's_filter_min            = ', s_filter_min
write (iu, '(a,  f11.4)') 's_filter_max            = ', s_filter_max
write (iu, '(a, es10.3)') 'sr3d_params%ds_track_step_max         = ', sr3d_params%ds_track_step_max
write (iu, '(a, es10.3)') 'sr3d_params%dr_track_step_max         = ', sr3d_params%dr_track_step_max
write (iu, '(a, es10.3)') 'surface_roughness_rms (input)         = ', surface_roughness_rms
write (iu, '(a, es10.3)') 'surface_roughness_rms (set value)     = ', wall%surface(1)%info%surface_roughness_rms
write (iu, '(a, es10.3)') 'roughness_correlation_len (input)     = ', roughness_correlation_len
write (iu, '(a, es10.3)') 'roughness_correlation_len (set value) = ', wall%surface(1)%info%roughness_correlation_len
write (iu, '(a, l1)') 'sr3d_params%allow_reflections         = ', sr3d_params%allow_reflections
write (iu, '(a, l1)') 'sr3d_params%specular_reflection_only  = ', sr3d_params%specular_reflection_only
write (iu, '(a, l1)') 'sr3d_params%stop_if_hit_antechamber   = ', sr3d_params%stop_if_hit_antechamber
write (iu, *)

end subroutine

end program
