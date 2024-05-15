!+ 
! Program dark_current_tracker
!
! Program to track multiple particles through an arbitrary lattice
! (backwards or forwards through elements) using a custom time-
! based tracking routine.
!
! Note: Partcles are assumed to be electrons
!
! Modules Needed:
!   use bmad
!
! Input (namelist)
!   dark_current_tracker.in
!
! Output
!   [particle file name].out
!-

program dark_current_tracker

use dark_current_mod
use quick_plot
use random_mod
use time_tracker_mod
use runge_kutta_mod ! for common struct only
!$ use omp_lib


implicit none

type (lat_struct), target :: lat
type (lat_struct), allocatable :: omp_lat(:)
type (branch_struct), pointer :: branch
type (ele_struct), pointer :: ele
type (track_struct) :: single_track
type (track_struct), allocatable :: omp_single_track(:)
type (coord_struct),  allocatable :: start_particles(:), omp_particles(:)
type (coord_struct) :: start,  end_particle, monitor_orb, dummy_particle
type (ele_pointer_struct), allocatable :: monitor_eles(:)
type (dark_current_param_struct) :: dc_param
type (dark_current_param_struct), allocatable :: omp_dc_param(:)
type (dark_current_tally_struct), allocatable :: tally(:)

character(2000) :: in_file, lat_name, lat2_name, lat_path, base_name, particle_file_name, monitor_file_name, outfile_name
character(400) :: element_monitor_list
real(rp) :: e_tot, max_charge, min_charge_ratio, dt_save
real(rp) :: color
real(rp) :: r1, plot_y_min, plot_y_max, plot_color_scale, plot_x_min, plot_x_max, plot_x_size, plot_y_size
real(rp) :: rel_tol_adaptive_tracking, abs_tol_adaptive_tracking
!real(rp) :: ref_p0c, ref_time
integer :: outfile, namelist_file, monitor_file
integer :: plot_every_n, species
integer, allocatable :: omp_counter(:), omp_particle_id(:)


integer :: ele_id, particle_id, iteration, point_id, plot_id, particle_set
integer :: ix_ele, n_monitor_eles
integer :: n_char
integer :: omp_n, omp_i, i, n_sets
logical :: save_tracks, verbose, err, save_field, global_frame, plot_on
logical :: parallel, error, do_tally
character(1), parameter :: carriage_return = achar(13) 
character(20) :: particle_species

character(30), parameter :: r_name = 'dark_current_tracker'

namelist / dark_current_tracker_params / lat_name, lat2_name, particle_file_name, save_tracks, &
		save_field, global_frame, verbose, min_charge_ratio, dt_save, element_monitor_list, plot_on, parallel, &
    plot_every_n, plot_y_min, plot_y_max, plot_color_scale, plot_x_min, plot_x_max, plot_x_size, plot_y_size, &
    rel_tol_adaptive_tracking, abs_tol_adaptive_tracking, do_tally, particle_species


!------
print '(a)',  ''
print '(a)',  ' __        __           __        __   __   ___      ___    ___  __        __        ___  __  ' 
print '(a)',  '|  \  /\  |__) |__/    /  ` |  | |__) |__) |__  |\ |  |      |  |__)  /\  /  ` |__/ |__  |__) '
print '(a)',  '|__/ /~~\ |  \ |  \    \__, \__/ |  \ |  \ |___ | \|  |      |  |  \ /~~\ \__, |  \ |___ |  \ '
print '(a)',  ''


!------------------------------------------
!Defaults for namelist
lat_name = 'lat.bmad'
lat2_name = ''
particle_file_name = 'ReferenceParticles.dat'
save_tracks = .false.
save_field  = .false. 
global_frame = .true. 
verbose = .false.
min_charge_ratio = 0
dt_save = 0.001/c_light
element_monitor_list = ''
parallel = .false.
do_tally = .true.
plot_on = .false.
plot_x_min = 0
plot_x_max = 0
plot_y_min = -.05_rp
plot_y_max =  .05_rp
plot_color_scale = .01_rp
plot_x_size = 800.0_rp
plot_y_size = 400.0_rp
rel_tol_adaptive_tracking = bmad_com%rel_tol_adaptive_tracking
abs_tol_adaptive_tracking = bmad_com%abs_tol_adaptive_tracking
particle_species = 'electron'

!Read namelist
in_file = 'dark_current_tracker.in'
! TODO: add this to the -init option. Need to restructure init
!if (command_argument_count() > 0) call get_command_argument(1, in_file)

namelist_file = lunget()
print *, 'Opening: ', trim(in_file)
open (namelist_file, file = in_file, status = "old")
read (namelist_file, nml = dark_current_tracker_params)
close (namelist_file)

!TEMP
dc_param%particle_file_name = particle_file_name
dc_param%save_field = save_field
dc_param%save_tracks = save_tracks
dc_param%verbose = verbose
dc_param%dt_save = dt_save
dc_param%global_frame = global_frame
dc_param%plot_on = plot_on
dc_param%id = 1




species = species_id(particle_species)
print *, 'Tracking particle type: '//trim(species_name(species))
print *, '  Mass (u)     = ', mass_of(species)/atomic_mass_unit
print *, '  Charge (|e|) = ', charge_of(species)

bmad_com%debug = verbose

! Read command line arguments
call dark_current_tracker_parse_command_args (dc_param, error)
if (error) then
  call out_io (s_fatal$, r_name, 'problem reading command line argument')
  stop
endif

print *, 'particle_file_name: ', trim(dc_param%particle_file_name)

! Force saving tracks with plot
if (dc_param%plot_on) dc_param%save_tracks = .true.


if ((.not. save_tracks) .and. save_field) then
  call out_io (s_warn$, r_name, 'WARNING: save_field only allowed with save_tracks. Disabling save_field')
  !print *, 'WARNING: save_field only allowed with save_tracks. Disabling save_field'
   dc_param%save_field = .false.
   save_field = .false.
endif


if (rel_tol_adaptive_tracking /= bmad_com%rel_tol_adaptive_tracking) then
  bmad_com%rel_tol_adaptive_tracking = rel_tol_adaptive_tracking
  print *, 'setting custom rel_tol_adaptive_tracking: ', bmad_com%rel_tol_adaptive_tracking
else
  print *, 'using default rel_tol_adaptive_tracking: ', bmad_com%rel_tol_adaptive_tracking
endif

if (abs_tol_adaptive_tracking /= bmad_com%abs_tol_adaptive_tracking) then
  bmad_com%abs_tol_adaptive_tracking = abs_tol_adaptive_tracking
  print *, 'setting custom abs_tol_adaptive_tracking: ', bmad_com%abs_tol_adaptive_tracking 
else
  print *, 'using default abs_tol_adaptive_tracking: ', bmad_com%abs_tol_adaptive_tracking  
endif

if(global_frame) then
  print *, 'writing in global frame'
else
  print *, 'writing in s-coordinates'
endif

if (save_tracks) then
   print *, "save_tracks ON; all particle tracks will be saved"
else
   print *, "save_tracks OFF; only final coordinates will be saved"
endif

if (save_field) print *, "save_field ON"


!Trim filename
n_char= SplitFileName(lat_name, lat_path, base_name) 

!Prepare output file name
call file_suffixer (dc_param%particle_file_name, outfile_name, '.out', .true.)

!Parse Lattice
call bmad_parser (lat_name, lat)
branch => lat%branch(0)

!Parse additional settings
if (lat2_name /= '') then
  print *, 'Parsing: '//trim(lat2_name)
  call bmad_parser2 (lat2_name, lat)
endif

! Check for absolute time tracking. If not, abort!
if (.not. bmad_com%absolute_time_tracking ) then
  call out_io (s_error$, r_name, 'absolute time tracking must be set to True.')
  stop
endif

!------------------------------------------
!Track through multiple elements

!Import  BMAD-T style particles
call import_time_distribution(dc_param%particle_file_name, start_particles)

!Switch all elements to time_runge_kutta$ tracking
do ele_id = 1, branch%n_ele_track
   if (branch%ele(ele_id)%key == taylor$) cycle ! Skip
   if (branch%ele(ele_id)%key == patch$) cycle ! Skip
   if (branch%ele(ele_id)%value(L$) == 0) cycle ! Skip
   branch%ele(ele_id)%tracking_method = time_runge_kutta$
   
   !Also set logic default for monitor use
   branch%ele(ele_id)%logic = .false.
end do

if (do_tally) then
  allocate(tally(branch%n_ele_track))
endif

!Monitor elements
if (element_monitor_list /= '') then  
  call lat_ele_locator (element_monitor_list, lat, monitor_eles, n_monitor_eles, err)
  if (err) then
  	call out_io (s_fatal$, r_name, 'PROBLEM WITH MONITOR ELES' // element_monitor_list)
    call err_exit
  end if
  
  print *, 'Monitoring: '
  do ix_ele = 1, n_monitor_eles
    monitor_eles(ix_ele)%ele%logic = .true.
    print *, '            ', trim( monitor_eles(ix_ele)%ele%name ), ' at s = ', monitor_eles(ix_ele)%ele%s, ' m'
  end do

  monitor_file = lunget()
  call file_suffixer (dc_param%particle_file_name, monitor_file_name, '.monitor', .true.)
  open(monitor_file, file = monitor_file_name)
  if (save_field) then
    write (monitor_file, '(10a19)') 'x', 'cp_x', 'y', 'cp_y', 's', 'cp_s', 't', 'charge', 'hit_angle', 's_origin'
    write (monitor_file, '(10a19)') 'm', 'eV', 	'm', 'eV',  'm', 'eV',  's', 'C',      'rad',  'm'  
  end if

  !Coordinates will be monitored at the exit end of elements in the loop below
end if



!Open output file
outfile = lunget()
open(outfile, file = outfile_name)

!Write header
if (save_tracks .or. save_field) then
  call write_particle_track(outfile, lat, write_field = save_field, write_header = .true., global_frame = global_frame)
else if (global_frame) then
  write (outfile, '(11a19)') 'x', 'cp_x', 'y', 'cp_y', 'z', 'cp_z', 't', 'charge', 'hit_angle', 's_origin', 'ix_ele'
  write (outfile, '(11a19)') 'm', 'eV', 	'm', 'eV',  'm', 'eV',   's', 'C',      'rad',       'm',        '1'
else
  write (outfile, '(11a19)') 'x', 'cp_x', 'y', 'cp_y', 's', 'cp_s', 't', 'charge', 'hit_angle', 's_origin', 'ix_ele'
  write (outfile, '(11a19)') 'm', 'eV', 	'm', 'eV',  'm', 'eV',   's', 'C',      'rad',       'm'  ,      '1'
end if



!--------------------
!Iterate through particles to find maximum charge
max_charge = maxval(start_particles%charge)

!Track electrons
branch%param%default_tracking_species = electron$

!write (*, '(a)', advance='no') 'Starting tracking'

! Set ds_save
single_track%ds_save = dc_param%dt_save*c_light

!--------------------
! Plotting
if (plot_on) then
  if (plot_x_max == 0) plot_x_max = lat%ele(lat%n_ele_track)%s

  call qp_open_page ('X', plot_id, plot_x_size, plot_y_size, 'POINTS')
  call qp_set_page_border (0.05_rp, 0.05_rp, 0.05_rp, 0.05_rp, '%PAGE')
  call qp_set_margin (0.1_rp, 0.1_rp, 0.1_rp, 0.1_rp, '%PAGE')
  call qp_set_line_width_basic(0)
  call qp_calc_and_set_axis ('X', plot_x_min , plot_x_max, 4, 8, 'GENERAL')
  call qp_calc_and_set_axis ('Y', plot_y_min, plot_y_max, 4, 8, 'GENERAL')
  !call plscmap1n(256)
  !call plspal1('cmap1_blue_yellow.pal',1)
  call qp_draw_axes ("Z (m)", "X (m)")
  call qp_set_line_attrib ('PLOT', color ='black')
  call qp_set_line_width_basic(0)
  call plot_wall(lat)
endif



if (.not. parallel) then
!--------------------
!Track list of particles
particle_loop: do particle_id = 1, size(start_particles)
  if (verbose) print *, "Tracking particle: ", particle_id 

  !If particle has charge less than max_charge * min_charge_ratio, move on
  if (start_particles(particle_id)%charge < max_charge * min_charge_ratio) &
    cycle

  !Set start coords and starting element; reset flags
  start = start_particles(particle_id)

  call prepare_particle (start, lat) 

  !----------------------------
  !      Track particle 
  !----------------------------
  end_particle = start
  
  if (dc_param%save_tracks .or. dc_param%plot_on) then
    call track_until_dead (start, lat, end_particle, single_track)
  else
    call track_until_dead (start, lat, end_particle)
  endif
  
  if (plot_on) call plot_particle_track( lat,single_track, plot_color_scale )

  ! Write to file
  ! Do not write tracks when making plots TEMPORARY
  if (dc_param%save_tracks .and. (.not. dc_param%plot_on)) then
    write (outfile, '(a, i8)')   'particle ', particle_id
    call write_particle_track (outfile, lat, track = single_track, write_field = dc_param%save_field, global_frame = dc_param%global_frame)   
  end if
  
  ! Clean track for next particle
  ! depreciated: call init_saved_orbit (single_track, 0)
  single_track%n_pt = -1

  if (.not. dc_param%save_tracks)  then
  ! Write to file: final point
    !Convert to t-coordinates for file output
    if (global_frame) then
      end_particle = particle_in_global_frame (end_particle, branch)
    else
       call convert_particle_coordinates_s_to_t(end_particle, end_particle%s, branch%ele(end_particle%ix_ele)%orientation)
    endif
    call add_to_tally(end_particle)
    write (outfile, '(10es19.10E3, i19)')  end_particle%vec, end_particle%t,  end_particle%charge, end_particle%phase(2), start%s,  end_particle%ix_ele
  end if

enddo particle_loop
endif



if (parallel) then
  print *, 'Parallel mode'
  omp_n = 1
  !$ omp_n = omp_get_max_threads()
  !$ print *, 'omp_get_max_threads(): ', omp_n

  !Set up independent lattices 
  allocate(omp_single_track(plot_every_n))

  allocate(omp_particle_id(plot_every_n))

  allocate(omp_dc_param(omp_n))
  allocate(omp_lat(omp_n))
  allocate(omp_counter(omp_n))
  do i=1, omp_n
    omp_lat(i) = lat
    !print *, 'parsing lattice: ', i
    !call bmad_parser (lat_name, omp_lat(i))
    !call bmad_parser2 (lat2_name, omp_lat(i))
    omp_dc_param(i) = dc_param
    omp_dc_param(i)%id = i
    omp_counter(i) = 0
  end do

  !Set up particle temp array
  allocate(omp_particles(plot_every_n))

  n_sets = ceiling(size(start_particles) / (1.0*plot_every_n) )
  print *, 'particle sets: ', n_sets


  particle_id = 1
  !dummy_particle = start_particles(particle_id)
  !call prepare_particle (dummy_particle, lat)

  do particle_set = 1, n_sets
    ! Prepare set
    print *, 'particle set: ', particle_set
    
    do i = 1, plot_every_n
      omp_particle_id(i) = particle_id
      start = start_particles(particle_id)
      call prepare_particle (start, lat)
      omp_particles(i) = start
      
      ! Check for last set and adjust chunk size
      if (particle_id == size(start_particles)) then
        plot_every_n = particle_id - (n_sets-1)*plot_every_n
        print *, 'setting plot_every_n = ', plot_every_n, ' for particle_id: ', particle_id
        exit
      endif
      particle_id = particle_id + 1
    end do    
    
    ! Track the set
    !$OMP parallel &
    !$OMP default(private), &
    !$OMP private(omp_i), &
    !$OMP shared(omp_n, omp_particles, omp_lat, omp_dc_param, omp_single_track), &
    !$OMP shared(plot_every_n, dc_param, dummy_particle, lat, lat2_name) &
    !$OMP shared(omp_counter)
    !$OMP master
      !print *, 'master'
      !call bmad_parser2(lat2_name, lat) 
    !$OMP end master
    !$OMP do schedule(dynamic)
    do i=1, plot_every_n
      omp_i=1
      !$ omp_i = omp_get_thread_num()+1
      !print *, 'omp_get_thread_num() + 1: ', omp_i
      !call init_saved_orbit (omp_single_track(i), 0)
      omp_single_track(i)%n_pt = -1
      !call old_track_until_dead (omp_particles(i), omp_lat(omp_i), omp_dc_param(omp_i), omp_single_track(i))
      if (dc_param%save_tracks .or. dc_param%plot_on) then
        call track_until_dead (omp_particles(i), omp_lat(omp_i), omp_particles(i), omp_single_track(i))
      else
        call track_until_dead (omp_particles(i), omp_lat(omp_i), omp_particles(i))
      endif
      !call track_until_dead (omp_particles(i), omp_lat(omp_i), omp_dc_param(omp_i))
      ! Worker counter 
      omp_counter(omp_i) = omp_counter(omp_i) + 1
    end do  
    !$OMP end do
    !$OMP end parallel
    
    if (dc_param%plot_on) then
      do i = 1, plot_every_n
        color = omp_single_track(i)%pt(0)%orb%vec(1)/0.01_rp
        !call plcol1(color)
        call plot_particle_track( lat, omp_single_track(i), plot_color_scale, in_global_frame = dc_param%global_frame )
      end do
    endif

  if (dc_param%save_tracks .and. (.not. dc_param%plot_on)) then
  
  end if

  if (.not. dc_param%save_tracks)  then
    do i=1, plot_every_n
      !Write to file: final point
      end_particle = omp_particles(i)
      !Convert to t-coordinates for file output
      if (global_frame) then
        end_particle = particle_in_global_frame (end_particle, branch)
      else
        call convert_particle_coordinates_s_to_t(end_particle, end_particle%s, branch%ele(end_particle%ix_ele)%orientation)
      endif
      call add_to_tally(end_particle)
      write (outfile, '(10es19.10E3, i19)')  end_particle%vec, end_particle%t,  end_particle%charge, &
             end_particle%phase(2), start_particles(omp_particle_id(i))%s,  end_particle%ix_ele
    end do
  ! TODO: FIX!!!
  else if (.not. dc_param%plot_on) then
    do i=1, plot_every_n
      write (outfile, '(a, i8)')   'particle ', omp_particle_id(i)
      call write_particle_track (outfile, lat, track = omp_single_track(i), write_field = dc_param%save_field, global_frame = dc_param%global_frame) 
    enddo
  end if

    !if (dc_param%save_tracks .and. (.not. dc_param%plot_on)) then
    !  do i = 1, omp_n
    !  write (outfile, '(a, i8)'),   'particle ', particle_id
    !  call write_particle_track (outfile, omp_lat(omp_i), track = omp_single_track(i), write_field = dc_param%save_field, global_frame = dc_param%global_frame)   
    !end if
    
  end do
  
  print *, 'Counter: (particles tracked in each thread)', omp_counter

endif

!call system('ls -a > filenames.txt')

!call ran_engine (set ='pseudo', ix_generator = 8)
!call ran_seed_put(seed=0, ix_generator = 8)
!call ran_engine (set ='pseudo', ix_generator = 9)
!call ran_seed_put(seed=0, ix_generator = 9)
!
!random_loop: do particle_id = 1, 1000
!  print *, particle_id
!  call init_coord (start)
!  call ran_uniform_scalar (r1, ix_generator = 8)
!  start%vec(1) = r1*0.01_rp
!  !print *, r1*100
!  
!  call ran_uniform_scalar (r1, ix_generator = 9)
!  !print *, r1*100
!  start%t = r1/1.3e9_rp
!  
!  
!  call track1_dark_current(start) 
!
!enddo random_loop

!Close output file
close(outfile)
print *, "Written: ", outfile_name


!Close monitor file if used
if (element_monitor_list /= '') then
  print *, "Written monitor file: ", monitor_file_name
  close(monitor_file)
end if


if (do_tally) then
 !Prepare output file name
  call file_suffixer (dc_param%particle_file_name, outfile_name, '.tally', .true.)
  outfile = lunget()
  open(outfile, file = outfile_name)
  write(outfile, '(4a10 )') 'name', 'ix_ele', 'charge', 'energy'  
  write(outfile, '(4a10 )') 'string', 'integer', 'C', 'J'  
  do ele_id = 1, branch%n_ele_track
    if (tally(ele_id)%charge == 0 ) cycle
    ele => lat%ele(ele_id)
    write(outfile, '(a, x, i0, 2es19.10E3 )') trim(ele%name), ele%ix_ele, tally(ele_id)%charge, tally(ele_id)%energy
  enddo
  close(outfile)
  print *, "Written: ", outfile_name
endif




!--------------------
! Plotting control
if (plot_on) then
   write (*, '(a)', advance = 'NO') ' Hit any key to end program: '
   read (*, '(a)') lat_name
  call qp_close_page 
endif


contains


subroutine prepare_particle (particle, lat)
! Prepares a particle in time-coordinates imported from a file, 
! Returns a proper s-based particle in-place
! 
!
!
!

type (coord_struct) :: particle
type (lat_struct) :: lat
type (ele_struct), pointer :: ele

real(rp) :: e_tot, pc2
integer ix_ele

!

particle%species = species
!if s < 0, drift particle forward to s = 0
if (particle%s < 0) call drift_orbit_time(particle, 1.0_rp,  delta_s =  -particle%s)
  
  ! Info to screen
  !print *, 'test' 
  !write (*, '(2a)', advance='no') carriage_return, 'def'
  !write (*, '(a, i8)', advance = 'no') carriage_return, particle_id
  !write (*, '( i8)', advance = 'no') particle_id
   
!Find starting element
ix_ele = element_at_s(lat, particle%s, .true., branch%ix_branch)
ele => branch%ele(ix_ele)
  
!Convert to global s-coordinates
!For the starting particle, the reference momentum and time taken at the end of the element. 
particle%state = alive$
particle%location = inside$
particle%ix_ele = ele%ix_ele
particle%p0c = ele%value(p0c$)
particle%direction = nint(sign(1.0_rp, particle%vec(6)))

! beta calc. TODO: fix this entire routine
pc2 = particle%vec(2)**2 + particle%vec(4)**2 +  particle%vec(6)**2
e_tot = sqrt(pc2 + mass_of(particle%species)**2) 
particle%beta = sqrt(pc2)/e_tot


call convert_particle_coordinates_t_to_s(particle, ele, ele%ref_time)




end subroutine prepare_particle


subroutine add_to_tally(particle)
type (coord_struct) :: particle
integer :: ix
if (.not. do_tally) return
ix = particle%ix_ele
tally(ix)%charge = tally(ix)%charge + particle%charge
tally(ix)%energy = tally(ix)%energy + particle%charge * sqrt(particle%vec(2)**2 +particle%vec(4)**2 + particle%vec(6)**2)


end subroutine

end program
