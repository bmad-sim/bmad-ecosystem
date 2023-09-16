!+ 
! Program parallel_track_example
!
! Example program to track particles in parallel using OpenMP
!
! Input (command line)
!   lat.bmad
!
! Output
!   lat.out
!-

program parallel_track_example

use bmad
use omp_lib

implicit none

type (lat_struct) lat
type (lat_struct),   allocatable :: omp_lat(:)
type (coord_struct), allocatable :: start_particles(:), end_particles(:)
type (ele_struct), pointer :: ele 
character(100) lat_name
character(100) outfile_name
integer :: i, omp_n, omp_i, n_particles
integer :: outfile = 1

type (coord_struct) particle_beg, particle_end

!------------------------------------------
!Get lattice file name
call getarg(1, lat_name)
print *,"Using lattice: ", lat_name

!Prepare output file name
call file_suffixer (lat_name, outfile_name, '.out', .true.)

!Parse Lattice
call bmad_parser (lat_name, lat)


!------------------------------------------
!Track through one element
ele => lat%ele(1) !First real element in lattice

call init_coord(particle_beg, ele=ele, element_end=upstream_end$)
particle_beg%state = alive$
particle_beg%vec = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6] !init
call track1(particle_beg, ele, lat%param, particle_end)

! Get number of threads. Set externally: export OMP_NUM_THREADS=12
omp_n = omp_get_max_threads()
print *, 'omp_get_max_threads(): ', omp_n

allocate(omp_lat(omp_n))
do i=1, omp_n
  omp_lat(i) = lat
end do

n_particles = 10000000

! Set up start/end particle arrays
allocate(start_particles(n_particles))
allocate(end_particles(n_particles))

do i=1, n_particles
  start_particles(i)=particle_beg
end do
print *, size(start_particles), ' particles'
print *, 'Particle Setup Done'

print *, 'Parallel Do'
!$OMP PARALLEL DO &
!$OMP DEFAULT(PRIVATE), &
!$OMP PRIVATE(omp_i), &
!$OMP SHARED(start_particles, n_particles, end_particles, omp_lat)
do i=1, n_particles
  omp_i = omp_get_thread_num()+1
  !print *, 'omp_get_thread_num() + 1: ', omp_i
  call track1(start_particles(i), omp_lat(omp_i)%ele(1), omp_lat(omp_i)%param, end_particles(i))
end do  
!$OMP END PARALLEL DO
print *, 'End Parallel Do'

deallocate(start_particles, end_particles)

!do i=1, n_particles
!	print *, end_particles(i)%vec
!end do


!------------------------------------------
!Write to file
!open(outfile, file = outfile_name)
!write (outfile, '(a)' )  "Start coordinates"
!write (outfile, '(6es18.10)') particle_beg%vec
!write (outfile, '(a)' )  "End coordinates"
!write (outfile, '(6es18.10)') particle_end%vec
!close(outfile)
!print *, "Written: ", outfile_name

end program
