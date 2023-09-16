!+ 
! Program bmad_to_opal_example
!
! Example program to make an OPAL input file from a Bmad lattice
!
! Input (command line)
!   bmad_to_opal_params namelist
!		lat_filename lat.bmad
!    
!
! Output
!   lat_opal.in
!-

program bmad_to_opal_example

use bmad
use beam_mod
use opal_interface_mod
use time_tracker_mod

implicit none

type (lat_struct) lat
character(100) in_filename
character(100) lat_filename
type (beam_init_struct) beam_init
type (beam_struct) beam
type (ele_struct), pointer :: ele1
character(100) opal_filename, opal_particle_filename
character(40)	:: r_name = 'bmad_to_opal_example'
integer :: n_arg
integer :: nmlfile, outfile, iu, ios
logical :: err

namelist / bmad_to_opal_params / lat_filename, beam_init

!------------------------------------------
!Setup

!Get units for files
nmlfile = lunget()
outfile = lunget()

!Get input file from command line						   
n_arg = command_argument_count()
if (n_arg > 1) then
  print *, 'Usage: bmad_to_opal_example <input_file>'
  print *, 'Default: <input_file> = bmad_to_opal_example.in'
  stop
endif

in_filename = 'bmad_to_opal_example.in' !default
if (n_arg == 1) call get_command_argument(1, in_filename)

! Defaults of input file
lat_filename = 'lat.bmad'
beam_init%n_particle = 1
beam_init%n_bunch = 1


!read input file
print *, 'Opening: ', trim(in_filename)
open (nmlfile, file = in_filename, status = "old")
read (nmlfile, nml = bmad_to_opal_params)
close (nmlfile)

!Prepare OPAL file
call file_suffixer (lat_filename, opal_filename, 'opal', .true.)
iu = lunget()	
open (iu, file = opal_filename, iostat = ios)
if (ios /= 0) then
	call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(opal_filename))
  	stop
endif


!Print info to screen
print *, '--------------------------------------'
write (*, '(2a)') 'Bmad lattice:      ', lat_filename
write (*, '(2a)') 'OPAL input file =  ', opal_filename
print *, '--------------------------------------'

!------------------------------------------
!Bmad to OPAL

!Parse Lattice
call bmad_parser (lat_filename, lat)

!Finally make OPAL input file
call write_opal_lattice_file(iu, lat)

!-------------------------------------------
!Cleanup

print *, "Written: ", opal_filename

!Close file
close(iu)


!-------------------------------------------
!Particle distribution
!Write particle file if more than one particle is defined
if (beam_init%n_particle > 1 ) then
  ele1 => lat%ele(0)

  !Initialize beam
  call init_beam_distribution (ele1, lat%param, beam_init, beam)


  !Prepare OPAL file
  call file_suffixer (lat_filename, opal_particle_filename, 'opal_particles', .true.)
  open (iu, file = opal_particle_filename, iostat = ios)
  if (ios /= 0) then
	  call out_io (s_error$, r_name, 'CANNOT OPEN FILE: ' // trim(opal_particle_filename))
  	  stop
  endif

  !Write the first bunch only
  call write_time_particle_distribution (iu, beam%bunch(1), ele1, style='OPAL')
  print *, "Written: ", opal_particle_filename
close(iu)


endif






	
	
end program
