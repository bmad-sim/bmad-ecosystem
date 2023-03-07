program bmad_to_gpt

use gpt_interface_mod
use beam_mod
use time_tracker_mod

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (beam_init_struct) beam_init
type (beam_struct) beam
type (ele_struct), pointer :: ele0
type (coord_struct) :: orb0
type (coord_struct), allocatable :: closed_orb(:)
type (gpt_lat_param_struct) :: gpt_lat_param
integer :: ios, iu
integer :: namelist_file, n_char

character(100) :: bmad_lat_filename, lat_path, base_name, in_file
character(100) :: gpt_particle_filename, time_particle_filename
character(*), parameter :: r_name = 'bmad_to_gpt'

logical :: err
logical :: write_bmad_time_particles, write_gpt_particles

namelist / bmad_to_gpt_params / &
    bmad_lat_filename, write_bmad_time_particles, write_gpt_particles, &
    gpt_lat_param, beam_init

!------------------------------------------
! Defaults for namelist

bmad_lat_filename = 'lat.bmad'
write_bmad_time_particles = .false.
write_gpt_particles = .false.
beam_init%n_particle = 1
beam_init%n_bunch = 1

! Read namelist

in_file = 'bmad_to_gpt.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)

namelist_file = lunget()
print *, 'Opening: ', trim(in_file)
open (namelist_file, file = in_file, status = "old")
read (namelist_file, nml = bmad_to_gpt_params)
close (namelist_file)

! Parse Lattice

call bmad_parser (bmad_lat_filename, lat)

! Prepare gpt file

if (gpt_lat_param%gpt_filename == '') then
  n_char= SplitFileName(bmad_lat_filename, lat_path, base_name) 
  call file_suffixer (base_name, gpt_lat_param%gpt_filename, 'gpt', .true.)
else
  n_char= SplitFileName(gpt_lat_param%gpt_filename, lat_path, base_name) 
endif

! Write file and stop if no particles are to be written.

call write_gpt_lattice_file(lat, gpt_lat_param, err)
if (err) stop

print *, 'Field map dimension: ', gpt_lat_param%fieldmap_dimension
print *, 'Written: ', gpt_lat_param%gpt_filename

if (.not. write_bmad_time_particles .and. .not. write_gpt_particles) stop

!-------------------------------------------
! Write particle file if more than one particle is defined

if (beam_init%n_particle > 1) then

  ! Initialize beam

  ele0 => lat%ele(0)
  call init_beam_distribution (ele0, lat%param, beam_init, beam)
  beam%bunch(1)%particle(:)%p0c = ele0%value(p0c_start$)
  call out_io (s_info$, r_name, 'Initialized bunch with p0c = \es13.3\ ', ele0%value(p0c$) )
  
  if (write_gpt_particles) then
  call file_suffixer (base_name, gpt_particle_filename, 'gpt_particles', .true.)
    iu = lunget()	
    open (iu, file = gpt_particle_filename, iostat = ios)
    if (ios /= 0) then
	    call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(gpt_particle_filename))
  	    stop
    endif

    !Write the first bunch only
    !call write_opal_particle_distribution (iu, beam%bunch(1), mass_of(lat%param%particle),  err)
    call write_time_particle_distribution (iu, beam%bunch(1), ele0, style = 'GPT', branch = lat%branch(0))
    print *, "Written gpt particles: ", gpt_particle_filename
    close(iu)
  endif
  
  if (write_bmad_time_particles) then
    iu = lunget()
    call file_suffixer (base_name, time_particle_filename, 'time_particles', .true.)
    open (iu, file = time_particle_filename, iostat = ios)
    if (ios /= 0) then
	    call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(time_particle_filename))
  	    stop
    endif

    !Write the first bunch only
    call write_time_particle_distribution (iu, beam%bunch(1), ele0, 'BMAD')
    print *, "Written bmad particles: ", time_particle_filename
    close(iu)
    
    ! Track to the end, write to file for reference
    !call track_beam (lat, beam, err=err)
    !call file_suffixer (base_name, time_particle_filename, 'end_particles', .true.)
    !open (iu, file = time_particle_filename, iostat = ios)
    !call write_time_particle_distribution (iu, beam%bunch(1))
    !print *, "Written: ", time_particle_filename
    !close(iu)   
  endif

endif

end program
