!+
! Program csr_example
!
! Simple program to deminstrate how to track with coherent synchrotron radiation.
!-

program csr_example

use bmad
use beam_mod
use random_mod

implicit none

type (lat_struct) lat
type (beam_init_struct) beam_init
type (beam_struct) beam
type (coord_struct), allocatable :: centroid(:)

real(rp) ave(6)
integer ran_seed, i
character(100) lat_file_name, input_file
logical err_flag

namelist / params / lat_file_name, beam_init, space_charge_com, bmad_com, ran_seed

! Get main input file name from the command line if specified.

if (command_argument_count() > 1) then
  print *, 'Too many command line arguments'
  stop
endif

input_file = 'csr.init'   ! Default
if (command_argument_count() == 1) call get_command_argument (1, input_file)

! Open main file and read namelist.
! Note: space_charge_com and bmad_com are global varibles

ran_seed = 0    ! Default

print *, 'Opening: ', trim(input_file)
open (1, file = input_file, status = 'old')
read (1, nml = params)
close (1)

! Some init

call ran_seed_put(ran_seed)

! Read in lattice

call bmad_parser (lat_file_name, lat)

! init beam

call init_beam_distribution (lat%ele(0), lat%param, beam_init, beam)

! Establish (approximate) centroid trajectory for CSR calculation

do i = 1, 6
  ave(i) = sum(beam%bunch(1)%particle%vec(i)) / size(beam%bunch(1)%particle%vec(i))
enddo
call reallocate_coord (centroid, lat)
call init_coord (centroid(0), ave, lat%ele(0), downstream_end$)
call track_all (lat, centroid)

! and track

call init_beam_distribution (lat%ele(0), lat%param, beam_init, beam)
call track_beam (lat, beam, err = err_flag, centroid = centroid)

! write results for first particle

print *, 'First particle coords at end of lattice:'
print '(5x, 6es15.5)', beam%bunch(1)%particle(1)%vec

open (1, file = 'csr.dat')
do i = 1, size(beam%bunch(1)%particle)
  write (1, '(i8, 6f14.8)') i, beam%bunch(1)%particle(i)%vec
enddo
close(1)

end program
