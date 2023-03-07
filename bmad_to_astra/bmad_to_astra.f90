program bmad_to_astra

use astra_interface_mod
use time_tracker_mod, only: write_time_particle_distribution
use bmad_struct
use beam_mod

type (lat_struct), target :: lat
type (branch_struct), pointer :: branch
type (beam_init_struct) beam_init
type (beam_struct) beam
type (ele_struct), pointer :: ele0, ele1
type (coord_struct) :: orb0
type (coord_struct), allocatable :: closed_orb(:)
type (astra_lattice_param_struct) :: astra_lattice_param
integer :: ios, iu
integer :: namelist_file, n_char, astra_file

character(100) :: lat_filename, lat2_filename, lat_path, base_name, in_file, astra_filename
character(100) :: astra_particle_filename, time_particle_filename
character(30), parameter :: r_name = 'bmad_to_astra'

logical :: err
logical :: write_time_particles, write_astra_particles

namelist / bmad_to_astra_params / &
    lat_filename, lat2_filename, write_time_particles, write_astra_particles, &
    astra_lattice_param, beam_init, astra_filename

!------------------------------------------
!Defaults for namelist
lat_filename = 'lat.bmad'
lat2_filename = ''
astra_filename = ''
write_time_particles = .false.
write_astra_particles = .false.
astra_lattice_param%fieldmap_dimension = 3
beam_init%n_particle = 1
beam_init%n_bunch = 1

!Read namelist
in_file = 'bmad_to_astra.in'
if (command_argument_count() > 0) call get_command_argument(1, in_file)

namelist_file = lunget()
print *, 'Opening: ', trim(in_file)
open (namelist_file, file = in_file, status = "old")
read (namelist_file, nml = bmad_to_astra_params)
close (namelist_file)

!Trim filename
n_char= SplitFileName(lat_filename, lat_path, base_name) 

!Parse Lattice
call bmad_parser (lat_filename, lat)

!Parse additional settings
if (lat2_filename /= '') then
  print *, 'Parsing: '//trim(lat2_filename)
  call bmad_parser2 (lat2_filename, lat)
endif

!Trim filename
n_char= SplitFileName(lat_filename, lat_path, base_name) 

!Prepare Astra file
if (astra_filename == '') then
  call file_suffixer (base_name, astra_filename, 'astra', .true.)
endif

astra_file = lunget()	
open (astra_file, file = astra_filename, iostat = ios)
if (ios /= 0) then
  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(astra_filename))
  stop
endif



! 
print *, 'Field map dimension: ', astra_lattice_param%fieldmap_dimension

call write_astra_lattice_file(astra_file, lat, astra_lattice_param, err)

print *, 'Written: ', astra_filename

! Stop if no particles are to be written
if ( (.not. write_time_particles) .and. (.not. write_astra_particles) ) then
  stop
endif

!-------------------------------------------
!Particle distribution
!Write particle file if more than one particle is defined
if (beam_init%n_particle > 1 ) then
  !set ele1 to be the init_ele
  ele1 => lat%ele(0)

  !Initialize beam
  call init_beam_distribution (ele1, lat%param, beam_init, beam)
  beam%bunch(1)%particle(:)%p0c = ele1%value(p0c_start$)
  call out_io (s_info$, r_name, 'Initialized bunch with p0c = \es13.3\', ele1%value(p0c$) )
  
  if (write_astra_particles) then
  call file_suffixer (base_name, astra_particle_filename, 'astra_particles', .true.)
    iu = lunget()	
    open (iu, file = astra_particle_filename, iostat = ios)
    if (ios /= 0) then
	    call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(astra_particle_filename))
  	    stop
    endif

    !Write the first bunch only    
    call write_time_particle_distribution  (iu, beam%bunch(1), ele1, 'ASTRA', lat%branch(0))
    
    print *, "Written astra particles: ", astra_particle_filename
    close(iu)
  endif
  
  if (write_time_particles) then
    iu = lunget()
    call file_suffixer (base_name, time_particle_filename, 'time_particles', .true.)
    open (iu, file = time_particle_filename, iostat = ios)
    if (ios /= 0) then
	    call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(time_particle_filename))
  	    stop
    endif

    !Write the first bunch only
    call write_time_particle_distribution  (iu, beam%bunch(1), ele1, 'BMAD', lat%branch(0))
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
